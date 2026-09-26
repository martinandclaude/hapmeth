//! Reads over a locus window -> per-read CpG calls on the reference.
//!
//! BAM (BAI/CSI) and CRAM (CRAI) through noodles. A call is kept when the read's
//! modified C sits on a reference CpG: forward reads at the C, reverse reads at
//! the G, both reported at the CpG's C (strands combined, as modkit
//! --combine-strands). Methylation is never taken from anything but MM/ML.

use std::{
    fs::File,
    io,
    path::{Path, PathBuf},
};

use anyhow::{bail, Context, Result};
use noodles::{
    bam, bgzf,
    core::{Position, Region},
    cram,
    csi::{self, binning_index::ReferenceSequence as _},
    fasta,
    sam::{
        self,
        alignment::record::{
            cigar::op::Kind,
            data::field::{value::Array, Tag, Value},
        },
    },
};

use crate::modbam;

const HP: Tag = Tag::new(b'H', b'P');
const PS: Tag = Tag::new(b'P', b'S');
const MM_LEGACY: Tag = Tag::new(b'M', b'm');
const ML_LEGACY: Tag = Tag::new(b'M', b'l');

#[derive(Clone, Copy, Debug)]
pub struct Filters {
    pub min_mapq: u8,
    pub include_supplementary: bool,
    /// Keep per-read call confidences for the threshold estimate.
    pub collect_conf: bool,
}

/// Reference bases for a window, with one extra base each side so the CpG
/// context of the edge positions can be checked.
#[derive(Clone, Debug)]
pub struct RefWindow {
    /// 0-based position of seq[0].
    pub offset: u64,
    pub seq: Vec<u8>,
}

impl RefWindow {
    pub fn base(&self, pos: u64) -> Option<u8> {
        pos.checked_sub(self.offset)
            .and_then(|i| self.seq.get(i as usize).copied())
    }

    /// CpG C positions in [start, end).
    pub fn cpgs(&self, start: u64, end: u64) -> Vec<u64> {
        (start..end)
            .filter(|&p| self.base(p) == Some(b'C') && self.base(p + 1) == Some(b'G'))
            .collect()
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Call {
    /// 0-based position of the CpG's C.
    pub cpg: u64,
    pub p_m: f32,
    pub p_h: f32,
    pub p_other: f32,
    pub p_c: f32,
    pub explicit: bool,
}

impl Call {
    /// Probability of the most likely state -- what the confidence filter compares.
    pub fn confidence(&self) -> f32 {
        self.p_c.max(self.p_m).max(self.p_h).max(self.p_other)
    }
}

#[derive(Clone, Debug)]
pub struct Read {
    pub name: String,
    /// 1 or 2 from the HP tag; 0 = untagged.
    pub hp: u8,
    pub ps: Option<i64>,
    pub reverse: bool,
    /// Aligned reference span, 0-based half-open.
    pub start: u64,
    pub end: u64,
    /// Sorted by position.
    pub calls: Vec<Call>,
    /// Confidence of every call at a reference CpG along the whole read (not
    /// just the window), implied `.` calls included as modkit does, for the
    /// threshold estimate.
    pub conf: Vec<f32>,
}

#[derive(Clone, Debug, Default)]
pub struct Skipped {
    pub mapq: usize,
    pub flags: usize,
    pub no_mod_tags: usize,
    pub bad_mod_tags: usize,
}

impl Skipped {
    pub fn add(&mut self, o: &Skipped) {
        self.mapq += o.mapq;
        self.flags += o.flags;
        self.no_mod_tags += o.no_mod_tags;
        self.bad_mod_tags += o.bad_mod_tags;
    }
}

enum BamIndex {
    Bai(bam::bai::Index),
    Csi(csi::Index),
}

enum Source {
    Bam(BamIndex),
    Cram {
        index: cram::crai::Index,
        repo: fasta::Repository,
    },
}

/// An indexed alignment file, opened once; `reader()` gives each worker thread
/// its own file handle.
pub struct Alignments {
    path: PathBuf,
    pub header: sam::Header,
    kind: Source,
}

fn with_ext(p: &Path, ext: &str) -> PathBuf {
    let mut s = p.as_os_str().to_owned();
    s.push(ext);
    PathBuf::from(s)
}

impl Alignments {
    pub fn open(path: &Path, reference: &Path) -> Result<Alignments> {
        let is_cram = path
            .extension()
            .is_some_and(|e| e.eq_ignore_ascii_case("cram"));
        if is_cram {
            let crai = with_ext(path, ".crai");
            let index = cram::crai::fs::read(&crai).with_context(|| {
                format!(
                    "reading {} (samtools index {})",
                    crai.display(),
                    path.display()
                )
            })?;
            let fa = fasta::io::indexed_reader::Builder::default()
                .build_from_path(reference)
                .with_context(|| format!("opening reference {}", reference.display()))?;
            let repo = fasta::Repository::new(fasta::repository::adapters::IndexedReader::new(fa));
            let mut r = cram::io::reader::Builder::default()
                .set_reference_sequence_repository(repo.clone())
                .build_from_path(path)
                .with_context(|| format!("opening {}", path.display()))?;
            let header = r
                .read_header()
                .with_context(|| format!("reading header of {}", path.display()))?;
            return Ok(Alignments {
                path: path.into(),
                header,
                kind: Source::Cram { index, repo },
            });
        }
        let mut r = File::open(path)
            .map(bam::io::Reader::new)
            .with_context(|| format!("opening {}", path.display()))?;
        let header = r
            .read_header()
            .with_context(|| format!("reading BAM header of {}", path.display()))?;
        // sample.bam.bai, sample.bai, sample.bam.csi
        let bai = [with_ext(path, ".bai"), path.with_extension("bai")]
            .into_iter()
            .find(|p| p.exists());
        let index = if let Some(bai) = bai {
            BamIndex::Bai(
                bam::bai::fs::read(&bai).with_context(|| format!("reading {}", bai.display()))?,
            )
        } else if with_ext(path, ".csi").exists() {
            BamIndex::Csi(csi::fs::read(with_ext(path, ".csi"))?)
        } else {
            bail!(
                "no index for {} (run: samtools index {})",
                path.display(),
                path.display()
            );
        };
        Ok(Alignments {
            path: path.into(),
            header,
            kind: Source::Bam(index),
        })
    }

    /// Mapped-read counts per contig from the BAM index (what samtools
    /// idxstats reads); None for CRAM or an index without that metadata.
    /// A contig with no reads has no metadata bin: it counts as 0.
    pub fn mapped_counts(&self) -> Option<Vec<(String, usize, u64)>> {
        let Source::Bam(index) = &self.kind else {
            return None;
        };
        let metas: Vec<Option<u64>> = match index {
            BamIndex::Bai(i) => i
                .reference_sequences()
                .iter()
                .map(|r| r.metadata().map(|m| m.mapped_record_count()))
                .collect(),
            BamIndex::Csi(i) => i
                .reference_sequences()
                .iter()
                .map(|r| r.metadata().map(|m| m.mapped_record_count()))
                .collect(),
        };
        if metas.iter().all(Option::is_none) {
            return None;
        }
        Some(
            self.header
                .reference_sequences()
                .iter()
                .zip(metas)
                .map(|((name, rs), m)| (name.to_string(), usize::from(rs.length()), m.unwrap_or(0)))
                .collect(),
        )
    }

    /// Release cached CRAM reference sequences (call between chromosomes).
    pub fn clear_reference_cache(&self) {
        if let Source::Cram { repo, .. } = &self.kind {
            repo.clear();
        }
    }

    pub fn reader(&self) -> Result<Reader<'_>> {
        let inner = match &self.kind {
            Source::Bam(_) => ReaderKind::Bam(File::open(&self.path).map(bam::io::Reader::new)?),
            Source::Cram { repo, .. } => ReaderKind::Cram(
                cram::io::reader::Builder::default()
                    .set_reference_sequence_repository(repo.clone())
                    .build_from_path(&self.path)?,
            ),
        };
        Ok(Reader { src: self, inner })
    }
}

enum ReaderKind {
    Bam(bam::io::Reader<bgzf::io::Reader<File>>),
    Cram(cram::io::Reader<File>),
}

pub struct Reader<'a> {
    src: &'a Alignments,
    inner: ReaderKind,
}

impl Reader<'_> {
    /// Reads overlapping [start, end) with their CpG calls inside that window.
    pub fn fetch(
        &mut self,
        chrom: &str,
        start: u64,
        end: u64,
        refw: &RefWindow,
        filters: Filters,
    ) -> Result<(Vec<Read>, Skipped)> {
        let region = Region::new(
            chrom,
            Position::try_from(start as usize + 1)?..=Position::try_from(end as usize)?,
        );
        let header = &self.src.header;
        let mut reads = Vec::new();
        let mut skipped = Skipped::default();
        match (&mut self.inner, &self.src.kind) {
            (ReaderKind::Bam(r), Source::Bam(index)) => {
                let query = match index {
                    BamIndex::Bai(i) => r.query(header, i, &region),
                    BamIndex::Csi(i) => r.query(header, i, &region),
                }
                .with_context(|| format!("querying {chrom}:{}-{end}", start + 1))?;
                for rec in query.records() {
                    let rec = rec?;
                    if let Some(read) = to_read(&rec, start, end, refw, filters, &mut skipped)? {
                        reads.push(read);
                    }
                }
            }
            (ReaderKind::Cram(r), Source::Cram { index, .. }) => {
                let query = r
                    .query(header, index, &region)
                    .with_context(|| format!("querying {chrom}:{}-{end}", start + 1))?;
                for rec in query.records() {
                    let rec = rec?;
                    if let Some(read) = to_read(&rec, start, end, refw, filters, &mut skipped)? {
                        reads.push(read);
                    }
                }
            }
            _ => unreachable!(),
        }
        Ok((reads, skipped))
    }
}

fn tag_int<R: sam::alignment::Record + ?Sized>(rec: &R, tag: &Tag) -> Option<i64> {
    rec.data()
        .get(tag)
        .and_then(|v| v.ok())
        .and_then(|v| v.as_int())
}

fn mod_tags<R: sam::alignment::Record + ?Sized>(rec: &R) -> io::Result<Option<(Vec<u8>, Vec<u8>)>> {
    let data = rec.data();
    let mm = match data
        .get(&Tag::BASE_MODIFICATIONS)
        .or_else(|| data.get(&MM_LEGACY))
    {
        Some(v) => match v? {
            Value::String(s) => s.to_vec(),
            _ => return Ok(None),
        },
        None => return Ok(None),
    };
    let ml = match data
        .get(&Tag::BASE_MODIFICATION_PROBABILITIES)
        .or_else(|| data.get(&ML_LEGACY))
    {
        Some(v) => match v? {
            Value::Array(Array::UInt8(vals)) => vals.iter().collect::<io::Result<Vec<u8>>>()?,
            _ => return Ok(None),
        },
        // An MM that lists no positions needs no ML ("C+m;").
        None => Vec::new(),
    };
    Ok(Some((mm, ml)))
}

fn to_read<R: sam::alignment::Record + ?Sized>(
    rec: &R,
    win_start: u64,
    win_end: u64,
    refw: &RefWindow,
    f: Filters,
    skipped: &mut Skipped,
) -> Result<Option<Read>> {
    let flags = rec.flags()?;
    if flags.is_unmapped()
        || flags.is_secondary()
        || flags.is_qc_fail()
        || flags.is_duplicate()
        || (flags.is_supplementary() && !f.include_supplementary)
    {
        skipped.flags += 1;
        return Ok(None);
    }
    // missing MAPQ (255) passes, as in samtools
    let mapq = rec.mapping_quality().transpose()?.map_or(255, u8::from);
    if mapq < f.min_mapq {
        skipped.mapq += 1;
        return Ok(None);
    }
    let Some(astart) = rec.alignment_start().transpose()? else {
        skipped.flags += 1;
        return Ok(None);
    };
    let astart = usize::from(astart) as u64 - 1;
    let ops: Vec<(Kind, usize)> = rec
        .cigar()
        .iter()
        .map(|op| op.map(|o| (o.kind(), o.len())))
        .collect::<io::Result<_>>()?;

    let Some((mm, ml)) = mod_tags(rec)? else {
        skipped.no_mod_tags += 1;
        return Ok(None);
    };
    let seq: Vec<u8> = rec.sequence().iter().collect();
    // MN (SEQ length when MM was written) guards against hard-clipped records
    // whose MM was not updated. Without MN, a hard clip is assumed stale.
    let hard_clipped = ops.iter().any(|(k, _)| *k == Kind::HardClip);
    let stale = match tag_int(rec, &Tag::BASE_MODIFICATION_SEQUENCE_LENGTH) {
        Some(mn) => mn != seq.len() as i64,
        None => hard_clipped,
    };
    if seq.is_empty() || stale {
        skipped.bad_mod_tags += 1;
        return Ok(None);
    }
    let reverse = flags.is_reverse_complemented();
    let mut mods = match modbam::cytosine_calls(&mm, &ml, &seq, reverse) {
        Ok(c) => c,
        Err(_) => {
            skipped.bad_mod_tags += 1;
            return Ok(None);
        }
    };
    mods.sort_unstable_by_key(|c| c.seq_index);

    let aend = astart
        + ops
            .iter()
            .filter(|(k, _)| k.consumes_reference())
            .map(|(_, l)| *l as u64)
            .sum::<u64>();
    let mut calls = Vec::new();
    let mut conf = Vec::new();
    let (mut q, mut r, mut i) = (0usize, astart, 0usize);
    for (kind, len) in ops {
        match kind {
            Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                while i < mods.len() && mods[i].seq_index < q + len {
                    let c = &mods[i];
                    if c.seq_index >= q {
                        let rp = r + (c.seq_index - q) as u64;
                        let cpg = if reverse {
                            (rp > 0
                                && refw.base(rp) == Some(b'G')
                                && refw.base(rp - 1) == Some(b'C'))
                            .then(|| rp - 1)
                        } else {
                            (refw.base(rp) == Some(b'C') && refw.base(rp + 1) == Some(b'G'))
                                .then_some(rp)
                        };
                        if let Some(cpg) = cpg {
                            let call = Call {
                                cpg,
                                p_m: c.p_m,
                                p_h: c.p_h,
                                p_other: c.p_other,
                                p_c: c.p_c,
                                explicit: c.explicit,
                            };
                            if f.collect_conf {
                                conf.push(call.confidence());
                            }
                            if (win_start..win_end).contains(&cpg) {
                                calls.push(call);
                            }
                        }
                    }
                    i += 1;
                }
                q += len;
                r += len as u64;
            }
            Kind::Insertion | Kind::SoftClip => {
                while i < mods.len() && mods[i].seq_index < q + len {
                    i += 1;
                }
                q += len;
            }
            Kind::Deletion | Kind::Skip => r += len as u64,
            Kind::HardClip | Kind::Pad => {}
        }
    }
    calls.sort_unstable_by_key(|c| c.cpg);
    let hp = match tag_int(rec, &HP) {
        Some(1) => 1,
        Some(2) => 2,
        _ => 0,
    };
    Ok(Some(Read {
        name: rec.name().map(|n| n.to_string()).unwrap_or_default(),
        hp,
        ps: tag_int(rec, &PS),
        reverse,
        start: astart,
        end: aend.max(astart + 1),
        calls,
        conf,
    }))
}

/// Reference around a window: the window plus `CONTEXT` bp each side, so calls
/// along most of each read can be checked for CpG context (threshold estimate).
pub const CONTEXT: u64 = 100_000;

/// Reference window [start-CONTEXT, end+CONTEXT) clamped to the contig.
pub fn fetch_ref(
    fa: &mut fasta::io::IndexedReader<fasta::io::BufReader<File>>,
    chrom: &str,
    start: u64,
    end: u64,
    chrom_len: u64,
) -> Result<RefWindow> {
    let s = start.saturating_sub(CONTEXT);
    let e = (end + CONTEXT).min(chrom_len);
    let region = Region::new(
        chrom,
        Position::try_from(s as usize + 1)?..=Position::try_from(e as usize)?,
    );
    let rec = fa
        .query(&region)
        .with_context(|| format!("reference {chrom}:{}-{e}", s + 1))?;
    let seq = rec
        .sequence()
        .as_ref()
        .iter()
        .map(u8::to_ascii_uppercase)
        .collect();
    Ok(RefWindow { offset: s, seq })
}
