//! MM/ML base-modification tags -> per-cytosine 5mC / 5hmC probabilities.
//!
//! Written from the SAM optional-fields spec (SAMtags §1.7), not from modkit;
//! tests/modkit_parity.rs checks it against `modkit extract` on real dorado
//! output. The points that matter:
//!
//! * Positions count bases of the stated type on the as-sequenced strand, so a
//!   reverse-mapped read (flag 0x10) is walked on the reverse complement of SEQ.
//! * `.` (or no flag) means skipped bases of that type are *unmodified*; `?`
//!   means they are *unknown*. dorado v5 writes `.`, and treating those implied
//!   calls as missing (as methylartist's pysam parse does) overstates 5mC.
//! * `C+hm` interleaves one ML value per code per position; separate `C+h.` and
//!   `C+m.` entries each list their own positions.
//! * Every entry consumes ML values, including ones we ignore (`A+a`, `C-m`),
//!   so they still have to be walked to keep ML aligned.

use std::fmt;

/// One cytosine of a read with a 5mC call.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CytosineCall {
    /// Index into SEQ as stored in the record (i.e. after any reverse complement).
    pub seq_index: usize,
    pub p_m: f32,
    pub p_h: f32,
    /// Summed probability of any other listed C modification (5fC, 5caC, `C+C`, ...).
    pub p_other: f32,
    /// P(unmodified C). For a listed call this is the leftover quality
    /// 255 - sum(ML), decoded like any ML value ((q + 0.5) / 256) -- modkit's
    /// convention, which decides calls sitting on the filter threshold.
    pub p_c: f32,
    /// false when the call is implied by a `.` skip rather than listed in ML.
    pub explicit: bool,
}

#[derive(Debug, PartialEq, Eq)]
pub enum MmError {
    Syntax(String),
    /// ML holds fewer values than MM lists.
    MlTooShort,
    /// MM points past the last base of that type in SEQ (e.g. SEQ was hard-clipped).
    BeyondSequence,
}

impl fmt::Display for MmError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            MmError::Syntax(e) => write!(f, "malformed MM tag: {e}"),
            MmError::MlTooShort => write!(f, "ML has fewer values than MM lists"),
            MmError::BeyondSequence => write!(f, "MM positions run past the end of SEQ"),
        }
    }
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum Code {
    M,
    H,
    Other,
}

fn code_of_letter(c: u8) -> Code {
    match c {
        b'm' => Code::M,
        b'h' => Code::H,
        _ => Code::Other,
    }
}

fn code_of_chebi(n: &[u8]) -> Code {
    match n {
        b"27551" => Code::M,
        b"76792" => Code::H,
        _ => Code::Other,
    }
}

fn complement(b: u8) -> u8 {
    match b.to_ascii_uppercase() {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' => b'A',
        b'U' => b'A',
        _ => b'N',
    }
}

/// One MM entry, e.g. `C+hm?,0,3` -> base C, strand +, codes [h, m], `?`.
struct Entry {
    base: u8,
    strand: u8,
    codes: Vec<Code>,
    /// `.` (or no flag): unlisted bases of this type are unmodified.
    implicit: bool,
    skips: Vec<usize>,
}

/// Parse one MM entry (without the trailing `;`).
fn parse_entry(e: &[u8]) -> Result<Entry, MmError> {
    let bad = || MmError::Syntax(String::from_utf8_lossy(e).into_owned());
    if e.len() < 3 {
        return Err(bad());
    }
    let (base, strand) = (e[0].to_ascii_uppercase(), e[1]);
    if !b"ACGTUN".contains(&base) || !(strand == b'+' || strand == b'-') {
        return Err(bad());
    }
    let rest = &e[2..];
    let code_end = rest
        .iter()
        .position(|&c| c == b',' || c == b'.' || c == b'?')
        .unwrap_or(rest.len());
    let code_str = &rest[..code_end];
    if code_str.is_empty() {
        return Err(bad());
    }
    let codes = if code_str.iter().all(u8::is_ascii_digit) {
        vec![code_of_chebi(code_str)]
    } else if code_str.iter().all(u8::is_ascii_alphabetic) {
        code_str.iter().map(|&c| code_of_letter(c)).collect()
    } else {
        return Err(bad());
    };
    let mut i = code_end;
    let mut implicit = true; // spec: no flag == '.'
    if i < rest.len() && (rest[i] == b'.' || rest[i] == b'?') {
        implicit = rest[i] == b'.';
        i += 1;
    }
    let mut skips = Vec::new();
    if i < rest.len() {
        if rest[i] != b',' {
            return Err(bad());
        }
        for tok in rest[i + 1..].split(|&c| c == b',') {
            let s = std::str::from_utf8(tok).map_err(|_| bad())?;
            skips.push(s.trim().parse::<usize>().map_err(|_| bad())?);
        }
    }
    Ok(Entry {
        base,
        strand,
        codes,
        implicit,
        skips,
    })
}

/// All C (as sequenced) calls carrying 5mC information, ordered by position on
/// the as-sequenced strand. Returns an empty list when MM has no `C+m` data.
///
/// `seq` is SEQ as stored; `reverse` is flag 0x10.
pub fn cytosine_calls(
    mm: &[u8],
    ml: &[u8],
    seq: &[u8],
    reverse: bool,
) -> Result<Vec<CytosineCall>, MmError> {
    let n = seq.len();
    let orig = |i: usize| -> u8 {
        if reverse {
            complement(seq[n - 1 - i])
        } else {
            seq[i].to_ascii_uppercase()
        }
    };
    let c_pos: Vec<usize> = (0..n).filter(|&i| orig(i) == b'C').collect();
    let nc = c_pos.len();

    let mut p_m: Vec<Option<f32>> = vec![None; nc];
    let mut p_h: Vec<Option<f32>> = vec![None; nc];
    let mut p_o: Vec<Option<f32>> = vec![None; nc];
    let mut explicit = vec![false; nc];
    let mut q_sum = vec![0u32; nc];
    let mut have_m = false;
    let mut ml_off = 0usize;

    for entry in mm.split(|&c| c == b';') {
        let entry = entry.trim_ascii();
        if entry.is_empty() {
            continue;
        }
        let Entry {
            base,
            strand,
            codes,
            implicit,
            skips,
        } = parse_entry(entry)?;
        let n_vals = skips.len() * codes.len();
        if ml_off + n_vals > ml.len() {
            return Err(MmError::MlTooShort);
        }
        if base != b'C' || strand != b'+' {
            ml_off += n_vals; // other base / opposite strand: skip its ML values
            continue;
        }
        have_m |= codes.contains(&Code::M);
        if implicit {
            for code in &codes {
                let arr = match code {
                    Code::M => &mut p_m,
                    Code::H => &mut p_h,
                    Code::Other => &mut p_o,
                };
                for v in arr.iter_mut().filter(|v| v.is_none()) {
                    *v = Some(0.0);
                }
            }
        }
        let mut k = 0usize;
        for skip in skips {
            k += skip;
            if k >= nc {
                return Err(MmError::BeyondSequence);
            }
            for code in &codes {
                let p = (f32::from(ml[ml_off]) + 0.5) / 256.0;
                q_sum[k] += u32::from(ml[ml_off]);
                ml_off += 1;
                match code {
                    Code::M => p_m[k] = Some(p),
                    Code::H => p_h[k] = Some(p),
                    Code::Other => p_o[k] = Some(p_o[k].unwrap_or(0.0) + p),
                }
            }
            explicit[k] = true;
            k += 1;
        }
    }
    if !have_m {
        return Ok(Vec::new());
    }
    let mut out = Vec::with_capacity(nc);
    for k in 0..nc {
        if let Some(m) = p_m[k] {
            let o = c_pos[k];
            out.push(CytosineCall {
                seq_index: if reverse { n - 1 - o } else { o },
                p_m: m,
                p_h: p_h[k].unwrap_or(0.0),
                p_other: p_o[k].unwrap_or(0.0),
                p_c: if explicit[k] {
                    (255u32.saturating_sub(q_sum[k]) as f32 + 0.5) / 256.0
                } else {
                    1.0
                },
                explicit: explicit[k],
            });
        }
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn p(v: u8) -> f32 {
        (f32::from(v) + 0.5) / 256.0
    }

    #[test]
    fn spec_example_implicit() {
        // SAMtags: C+m,5,12,0 -> 6th, 19th and 20th C listed, all others unmodified.
        let seq = "C".repeat(25);
        let calls = cytosine_calls(b"C+m,5,12,0;", &[200, 10, 250], seq.as_bytes(), false).unwrap();
        assert_eq!(calls.len(), 25);
        let listed: Vec<_> = calls
            .iter()
            .filter(|c| c.explicit)
            .map(|c| c.seq_index)
            .collect();
        assert_eq!(listed, vec![5, 18, 19]);
        assert_eq!(calls[5].p_m, p(200));
        assert_eq!(calls[0].p_m, 0.0);
        assert!(!calls[0].explicit);
        assert_eq!(calls[0].p_c, 1.0);
        // listed: leftover quality 255 - 200 = 55 -> (55 + 0.5) / 256
        assert_eq!(calls[5].p_c, 55.5 / 256.0);
    }

    #[test]
    fn spec_example_unknown() {
        let seq = "C".repeat(25);
        let calls =
            cytosine_calls(b"C+m?,5,12,0;", &[200, 10, 250], seq.as_bytes(), false).unwrap();
        assert_eq!(
            calls.iter().map(|c| c.seq_index).collect::<Vec<_>>(),
            vec![5, 18, 19]
        );
    }

    #[test]
    fn multi_code_interleaved_equals_separate_entries() {
        let seq = b"ACGTCCGACG";
        // Cs at 1,4,5,8. C+mh,1,1 lists the 2nd and 4th C (indices 4 and 8).
        let a = cytosine_calls(b"C+mh?,1,1;", &[204, 26, 89, 130], seq, false).unwrap();
        let b = cytosine_calls(b"C+m?,1,1;C+h?,1,1;", &[204, 89, 26, 130], seq, false).unwrap();
        assert_eq!(a, b);
        assert_eq!(a.len(), 2);
        assert_eq!((a[0].seq_index, a[0].p_m, a[0].p_h), (4, p(204), p(26)));
        assert_eq!((a[1].seq_index, a[1].p_m, a[1].p_h), (8, p(89), p(130)));
    }

    #[test]
    fn dorado_style_separate_implicit_entries() {
        // h lists only the first C; m lists first and third. Second C is implied
        // unmodified for both codes; third has h implied 0.
        let seq = b"CGCGCG";
        let calls = cytosine_calls(b"C+h.,0;C+m.,0,1;", &[20, 200, 180], seq, false).unwrap();
        assert_eq!(calls.len(), 3);
        assert_eq!(
            (calls[0].p_m, calls[0].p_h, calls[0].explicit),
            (p(200), p(20), true)
        );
        assert_eq!(
            (calls[1].p_m, calls[1].p_h, calls[1].explicit),
            (0.0, 0.0, false)
        );
        assert_eq!(
            (calls[2].p_m, calls[2].p_h, calls[2].explicit),
            (p(180), 0.0, true)
        );
    }

    #[test]
    fn other_bases_consume_ml() {
        // A+a lists two As before the C entry; their ML values must be skipped.
        let seq = b"AACG";
        let calls = cytosine_calls(b"A+a.,0,0;C+m.,0;", &[1, 2, 240], seq, false).unwrap();
        assert_eq!(calls.len(), 1);
        assert_eq!(calls[0].p_m, p(240));
        // ... also when the C entry comes first
        let calls = cytosine_calls(b"C+m.,0;A+a.,0,0;", &[240, 1, 2], seq, false).unwrap();
        assert_eq!(calls[0].p_m, p(240));
    }

    #[test]
    fn reverse_strand_counts_on_original_orientation() {
        // Stored SEQ CGTACG, reverse-mapped: as sequenced it is CGTACG
        // (revcomp of CGTACG). Original Cs at 0 and 4 -> stored indices 5 and 1.
        let seq = b"CGTACG";
        let calls = cytosine_calls(b"C+m?,1;", &[250], seq, true).unwrap();
        assert_eq!(calls.len(), 1);
        assert_eq!(calls[0].seq_index, 1); // original index 4 -> 6-1-4
        assert_eq!(seq[calls[0].seq_index], b'G'); // a C on the other strand
    }

    #[test]
    fn chebi_codes() {
        let seq = b"CG";
        let calls = cytosine_calls(b"C+27551?,0;C+76792?,0;", &[200, 30], seq, false).unwrap();
        assert_eq!((calls[0].p_m, calls[0].p_h), (p(200), p(30)));
    }

    #[test]
    fn explicit_absent_modification() {
        // "C+m;" : no positions, default '.', so every C is unmodified.
        let calls = cytosine_calls(b"C+m;", &[], b"CCGC", false).unwrap();
        assert_eq!(calls.len(), 3);
        assert!(calls.iter().all(|c| c.p_m == 0.0 && !c.explicit));
    }

    #[test]
    fn no_5mc_information() {
        assert!(cytosine_calls(b"A+a.,0;", &[3], b"ACG", false)
            .unwrap()
            .is_empty());
        assert!(cytosine_calls(b"", &[], b"ACG", false).unwrap().is_empty());
        // C+h alone carries no 5mC call
        assert!(cytosine_calls(b"C+h?,0;", &[3], b"ACG", false)
            .unwrap()
            .is_empty());
    }

    #[test]
    fn errors() {
        assert_eq!(
            cytosine_calls(b"C+m?,5;", &[3], b"CG", false),
            Err(MmError::BeyondSequence)
        );
        assert_eq!(
            cytosine_calls(b"C+m?,0,0;", &[3], b"CCG", false),
            Err(MmError::MlTooShort)
        );
        assert!(matches!(
            cytosine_calls(b"C*m,0;", &[3], b"CG", false),
            Err(MmError::Syntax(_))
        ));
        assert!(matches!(
            cytosine_calls(b"C+m?,x;", &[3], b"CG", false),
            Err(MmError::Syntax(_))
        ));
    }
}
