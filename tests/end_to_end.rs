//! End to end on a synthetic modBAM: every read is built so the expected
//! numbers are known exactly. Covers the cases that decide correctness --
//! dorado-style sparse `.` tags (implied calls count as unmethylated),
//! reverse-strand reads, `?` mode, an indel over a CpG, a CpG inside an
//! insertion, and the MAPQ / hard-clip / MN filters.

use std::{
    fs,
    path::{Path, PathBuf},
    process::Command,
};

use noodles::{
    bam,
    sam::{self, alignment::io::Write as _},
};

const LOCUS: (u64, u64) = (480, 720);

/// 2 kb reference: no C/G except CpGs at 500, 520, ..., 700 and a CA at 510.
fn reference() -> Vec<u8> {
    let mut r: Vec<u8> = b"ATTA".iter().copied().cycle().take(2000).collect();
    for p in (500..=700).step_by(20) {
        r[p] = b'C';
        r[p + 1] = b'G';
    }
    r[510] = b'C';
    r
}

fn revcomp(s: &[u8]) -> Vec<u8> {
    s.iter()
        .rev()
        .map(|b| match b {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            _ => b'N',
        })
        .collect()
}

/// MM skip list + ML values listing `calls` (original-orientation read indices
/// of Cs, with their ML value) among all Cs of `orig`.
fn mm_entry(code: &str, mode: char, orig: &[u8], calls: &[(usize, u8)]) -> (String, Vec<u8>) {
    let cs: Vec<usize> = (0..orig.len()).filter(|&i| orig[i] == b'C').collect();
    let mut s = format!("C+{code}{mode}");
    let mut ml = Vec::new();
    let mut prev = 0usize;
    for (i, (idx, v)) in calls.iter().enumerate() {
        let k = cs
            .iter()
            .position(|c| c == idx)
            .expect("listed position is not a C");
        let skip = if i == 0 { k } else { k - prev - 1 };
        s += &format!(",{skip}");
        prev = k;
        ml.push(*v);
    }
    s.push(';');
    (s, ml)
}

struct R {
    name: &'static str,
    flag: u16,
    pos: usize, // 1-based
    mapq: u8,
    cigar: String,
    seq: Vec<u8>,
    tags: String,
}

fn sam_line(r: &R) -> String {
    format!(
        "{}\t{}\tchrT\t{}\t{}\t{}\t*\t0\t0\t{}\t*\t{}",
        r.name,
        r.flag,
        r.pos,
        r.mapq,
        r.cigar,
        String::from_utf8_lossy(&r.seq),
        r.tags
    )
}

fn ml_tag(v: &[u8]) -> String {
    if v.is_empty() {
        String::new()
    } else {
        format!(
            "\tML:B:C,{}",
            v.iter().map(u8::to_string).collect::<Vec<_>>().join(",")
        )
    }
}

fn build(dir: &Path) -> (PathBuf, PathBuf, PathBuf) {
    let reference = reference();
    let fa = dir.join("ref.fa");
    fs::write(
        &fa,
        format!(">chrT\n{}\n", String::from_utf8_lossy(&reference)),
    )
    .unwrap();
    fs::write(dir.join("ref.fa.fai"), "chrT\t2000\t6\t2000\t2001\n").unwrap();

    let seg = reference[400..800].to_vec();
    let cpg_read_idx: Vec<usize> = (500..=700).step_by(20).map(|p| p - 400).collect();
    let mut reads = Vec::new();

    // HP1: 4 forward reads, every CpG listed methylated; the CA C is implied unmodified.
    let (h, mut ml_h) = mm_entry(
        "h",
        '.',
        &seg,
        &cpg_read_idx.iter().map(|&i| (i, 2)).collect::<Vec<_>>(),
    );
    let (m, ml_m) = mm_entry(
        "m",
        '.',
        &seg,
        &cpg_read_idx.iter().map(|&i| (i, 250)).collect::<Vec<_>>(),
    );
    ml_h.extend(&ml_m);
    for name in ["hp1_a", "hp1_b", "hp1_c", "hp1_d"] {
        reads.push(R {
            name,
            flag: 0,
            pos: 401,
            mapq: 60,
            cigar: "400M".into(),
            seq: seg.clone(),
            tags: format!("MM:Z:{h}{m}{}\tHP:i:1\tPS:i:7\tMN:i:400", ml_tag(&ml_h)),
        });
    }
    // HP1 with MAPQ 5: dropped by the default --min-mapq 10
    reads.push(R {
        name: "hp1_lowmapq",
        flag: 0,
        pos: 401,
        mapq: 5,
        cigar: "400M".into(),
        seq: seg.clone(),
        tags: format!("MM:Z:{h}{m}{}\tHP:i:1\tPS:i:7", ml_tag(&ml_h)),
    });
    // HP1 with the CpG at 600 deleted: 10 methylated calls
    let del_seq: Vec<u8> = [&reference[400..599], &reference[602..800]].concat();
    let del_calls: Vec<(usize, u8)> = (500..=700)
        .step_by(20)
        .filter(|p| *p != 600)
        .map(|p| (if p < 600 { p - 400 } else { p - 403 }, 250))
        .collect();
    let (dm, dml) = mm_entry("m", '.', &del_seq, &del_calls);
    reads.push(R {
        name: "hp1_del",
        flag: 0,
        pos: 401,
        mapq: 60,
        cigar: "199M3D198M".into(),
        seq: del_seq,
        tags: format!("MM:Z:{dm}{}\tHP:i:1\tPS:i:7", ml_tag(&dml)),
    });

    // HP2: 4 reverse reads with NO listed positions ("C+m.;"): every C is implied
    // unmethylated -- the case methylartist's parse gets wrong.
    for name in ["hp2_a", "hp2_b", "hp2_c", "hp2_d"] {
        reads.push(R {
            name,
            flag: 16,
            pos: 401,
            mapq: 60,
            cigar: "400M".into(),
            seq: seg.clone(),
            tags: "MM:Z:C+h.;C+m.;\tHP:i:2\tPS:i:7".into(),
        });
    }
    // HP2 with a CpG inside an insertion, listed as methylated: must be ignored
    let ins_seq: Vec<u8> = [
        &reference[400..650],
        b"ACGT".as_slice(),
        &reference[650..800],
    ]
    .concat();
    let (im, iml) = mm_entry("m", '.', &ins_seq, &[(251, 250)]);
    reads.push(R {
        name: "hp2_ins",
        flag: 0,
        pos: 401,
        mapq: 60,
        cigar: "250M4I150M".into(),
        seq: ins_seq,
        tags: format!("MM:Z:{im}{}\tHP:i:2\tPS:i:7", ml_tag(&iml)),
    });
    // HP2, hard-clipped after tagging (no MN): must be skipped, although it claims methylation
    let clipped = seg[10..].to_vec();
    reads.push(R {
        name: "hp2_clipped",
        flag: 0,
        pos: 411,
        mapq: 60,
        cigar: "10H390M".into(),
        seq: clipped,
        tags: format!("MM:Z:{m}{}\tHP:i:2", ml_tag(&ml_m)),
    });
    // HP2 whose MN disagrees with SEQ: skipped
    reads.push(R {
        name: "hp2_badmn",
        flag: 0,
        pos: 401,
        mapq: 60,
        cigar: "400M".into(),
        seq: seg.clone(),
        tags: format!("MM:Z:{m}{}\tHP:i:2\tMN:i:999", ml_tag(&ml_m)),
    });

    // untagged, '?' mode, alternating 240 / 10 over the 11 CpGs -> 6 of 11 methylated
    let alt: Vec<(usize, u8)> = cpg_read_idx
        .iter()
        .enumerate()
        .map(|(k, &i)| (i, if k % 2 == 0 { 240 } else { 10 }))
        .collect();
    let (um, uml) = mm_entry("m", '?', &seg, &alt);
    reads.push(R {
        name: "untagged",
        flag: 0,
        pos: 401,
        mapq: 60,
        cigar: "400M".into(),
        seq: seg.clone(),
        tags: format!("MM:Z:{um}{}", ml_tag(&uml)),
    });
    // untagged, '?' with nothing listed: no calls at all
    reads.push(R {
        name: "untagged_empty",
        flag: 16,
        pos: 401,
        mapq: 60,
        cigar: "400M".into(),
        seq: seg.clone(),
        tags: "MM:Z:C+m?;".into(),
    });
    // reverse-orientation sanity: stored SEQ is the reference strand, as sequenced it is the revcomp
    assert_eq!(revcomp(&revcomp(&seg)), seg);

    reads.sort_by_key(|r| r.pos);
    let header: sam::Header = "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chrT\tLN:2000\n"
        .parse()
        .unwrap();
    let bam_path = dir.join("syn.bam");
    let mut w = bam::io::Writer::new(fs::File::create(&bam_path).unwrap());
    w.write_header(&header).unwrap();
    for r in &reads {
        let line = sam_line(r);
        let rec = sam::Record::try_from(line.as_bytes()).unwrap_or_else(|e| panic!("{e}: {line}"));
        w.write_alignment_record(&header, &rec).unwrap();
    }
    w.try_finish().unwrap();
    drop(w);
    let index = bam::fs::index(&bam_path).unwrap();
    bam::bai::fs::write(dir.join("syn.bam.bai"), &index).unwrap();

    let bed = dir.join("test.bed");
    fs::write(
        &bed,
        format!("#chrom\tstart\tend\tname\tscore\tstrand\tgene\tdisease\torigin\texpected_lo\texpected_hi\tgroup\n\
                 chrT\t{}\t{}\tTEST_DMR\t0\t.\tTG\tsynthetic\t.\t0\t100\tSynthetic\n", LOCUS.0, LOCUS.1),
    )
    .unwrap();
    (bam_path, fa, bed)
}

fn run(args: &[&str]) -> std::process::Output {
    Command::new(env!("CARGO_BIN_EXE_hapmeth"))
        .args(args)
        .output()
        .unwrap()
}

fn tsv(path: &Path) -> Vec<std::collections::HashMap<String, String>> {
    let text = fs::read_to_string(path).unwrap();
    let mut lines = text.lines().filter(|l| !l.starts_with('#'));
    let cols: Vec<&str> = lines.next().unwrap().split('\t').collect();
    lines
        .map(|l| {
            cols.iter()
                .map(|c| c.to_string())
                .zip(l.split('\t').map(str::to_string))
                .collect()
        })
        .collect()
}

#[test]
fn synthetic_modbam() {
    let dir = std::env::temp_dir().join(format!("hapmeth_e2e_{}", std::process::id()));
    fs::create_dir_all(&dir).unwrap();
    let (bam, fa, bed) = build(&dir);
    let out = dir.join("out");
    let o = run(&[
        "-b",
        bam.to_str().unwrap(),
        "-r",
        fa.to_str().unwrap(),
        "-o",
        out.to_str().unwrap(),
        "-s",
        "syn",
        "--build",
        "hg38",
        "--panels",
        "none",
        "--bed",
        bed.to_str().unwrap(),
        "--filter-threshold",
        "0.6",
        "--sites",
        "--format",
        "both",
    ]);
    assert!(o.status.success(), "{}", String::from_utf8_lossy(&o.stderr));

    let rows = tsv(&out.join("syn.loci.tsv"));
    assert_eq!(rows.len(), 1);
    let r = &rows[0];
    let get = |k: &str| r[k].as_str();
    assert_eq!(get("name"), "TEST_DMR");
    assert_eq!(get("core_cpgs"), "11");
    // HP1: 4 full reads + the deletion read; the MAPQ-5 read is out
    assert_eq!(
        (get("hp1_reads"), get("hp1_cpgs"), get("hp1_pct")),
        ("5", "11", "100.0")
    );
    // HP2: implied '.' calls are unmethylated calls (not missing ones); the
    // insertion CpG and the clipped / bad-MN reads contribute nothing
    assert_eq!(
        (get("hp2_reads"), get("hp2_cpgs"), get("hp2_pct")),
        ("5", "11", "0.0")
    );
    assert_eq!((get("untagged_reads"), get("untagged_pct")), ("1", "54.5"));
    // all: 54 methylated of 4*11 + 10 + 5*11 + 11 = 120 calls
    assert_eq!(
        (get("all_reads"), get("all_pct"), get("hp_diff")),
        ("11", "50.0", "100.0")
    );
    assert_eq!(
        (
            get("expected_lo"),
            get("expected_hi"),
            get("phase_sets"),
            get("status")
        ),
        ("0.0", "100.0", "1", "ok")
    );
    assert_eq!(get("plot"), "plots/Synthetic/TEST_DMR.png");

    // per-CpG: the deleted CpG has one HP1 call fewer; reverse-strand calls land on the C
    let sites = tsv(&out.join("syn.cpgs.tsv"));
    let site = |pos: &str, g: &str| {
        sites
            .iter()
            .find(|s| s["pos"] == pos && s["group"] == g)
            .cloned()
            .unwrap()
    };
    assert_eq!(site("600", "hp1")["n_5mC"], "4");
    assert_eq!(site("620", "hp1")["n_5mC"], "5");
    assert_eq!(
        (
            site("600", "hp2")["n_C"].as_str(),
            site("600", "hp2")["n_5mC"].as_str()
        ),
        ("5", "0")
    );
    assert!(
        sites.iter().all(|s| s["pos"] != "510"),
        "the CA C is not a CpG"
    );

    for f in [
        "syn.overview.png",
        "syn.overview.svg",
        "plots/Synthetic/TEST_DMR.png",
        "plots/Synthetic/TEST_DMR.svg",
    ] {
        let bytes = fs::read(out.join(f)).unwrap_or_else(|_| panic!("{f} missing"));
        if f.ends_with(".png") {
            assert_eq!(&bytes[..8], b"\x89PNG\r\n\x1a\n", "{f} is not a PNG");
        } else {
            assert!(bytes.starts_with(b"<svg"), "{f} is not an SVG");
        }
    }
    fs::remove_dir_all(&dir).ok();
}

#[test]
fn default_threshold_and_min_mapq() {
    let dir = std::env::temp_dir().join(format!("hapmeth_e2e_mapq_{}", std::process::id()));
    fs::create_dir_all(&dir).unwrap();
    let (bam, fa, bed) = build(&dir);
    let out = dir.join("out");
    let o = run(&[
        "-b",
        bam.to_str().unwrap(),
        "-r",
        fa.to_str().unwrap(),
        "-o",
        out.to_str().unwrap(),
        "-s",
        "syn",
        "--build",
        "hg38",
        "--panels",
        "none",
        "--bed",
        bed.to_str().unwrap(),
        "--min-mapq",
        "0",
    ]);
    assert!(o.status.success(), "{}", String::from_utf8_lossy(&o.stderr));
    let r = &tsv(&out.join("syn.loci.tsv"))[0];
    assert_eq!(r["hp1_reads"], "6", "--min-mapq 0 keeps the MAPQ-5 read");
    assert!(String::from_utf8_lossy(&o.stderr).contains("call-confidence threshold"));
    fs::remove_dir_all(&dir).ok();
}

#[test]
fn export_panels_and_errors() {
    let dir = std::env::temp_dir().join(format!("hapmeth_e2e_export_{}", std::process::id()));
    let o = run(&["--export-panels", dir.to_str().unwrap()]);
    assert!(o.status.success());
    let mut names: Vec<String> = fs::read_dir(&dir)
        .unwrap()
        .map(|e| e.unwrap().file_name().to_string_lossy().into_owned())
        .collect();
    names.sort();
    assert_eq!(
        names,
        [
            "imprinted_germline_dmrs.hg38.bed",
            "imprinted_germline_dmrs.hs1.bed",
            "qc_controls.hg38.bed",
            "qc_controls.hs1.bed",
            "xci_panel.hg38.bed",
            "xci_panel.hs1.bed",
        ]
    );
    fs::remove_dir_all(&dir).ok();

    let o = run(&["-b", "/nonexistent.bam", "-r", "/nonexistent.fa"]);
    assert!(!o.status.success());
    assert!(String::from_utf8_lossy(&o.stderr).starts_with("error:"));
}
