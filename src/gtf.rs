//! Optional gene track: transcripts from a GTF (plain or gzipped) that overlap
//! the plotted windows. One transcript per gene is drawn: the one tagged
//! MANE_Select / Ensembl_canonical, else the longest.

use std::{
    collections::HashMap,
    fs::File,
    io::{BufRead, BufReader},
    path::Path,
};

use anyhow::Result;
use flate2::read::MultiGzDecoder;

#[derive(Clone, Debug, PartialEq)]
pub struct Transcript {
    pub gene: String,
    pub chrom: String,
    pub strand: char,
    /// 0-based half-open.
    pub start: u64,
    pub end: u64,
    pub exons: Vec<(u64, u64)>,
    canonical: bool,
}

#[derive(Default)]
pub struct Genes {
    by_chrom: HashMap<String, Vec<Transcript>>,
}

impl Genes {
    /// Transcripts overlapping [start, end) on `chrom` (with or without "chr").
    pub fn overlapping(&self, chrom: &str, start: u64, end: u64) -> Vec<&Transcript> {
        let alt = chrom
            .strip_prefix("chr")
            .map(str::to_string)
            .unwrap_or_else(|| format!("chr{chrom}"));
        [chrom, alt.as_str()]
            .iter()
            .filter_map(|c| self.by_chrom.get(*c))
            .flatten()
            .filter(|t| t.start < end && start < t.end)
            .collect()
    }
}

fn attr<'a>(attrs: &'a str, key: &str) -> Option<&'a str> {
    attrs.split(';').find_map(|kv| {
        let kv = kv.trim();
        let rest = kv.strip_prefix(key)?.trim_start();
        Some(rest.trim_matches('"'))
    })
}

/// Load transcripts overlapping any of `spans` (chrom, start, end).
pub fn load(path: &Path, spans: &[(String, u64, u64)]) -> Result<Genes> {
    let f = File::open(path)?;
    let reader: Box<dyn BufRead> = if path.extension().is_some_and(|e| e == "gz") {
        Box::new(BufReader::new(MultiGzDecoder::new(f)))
    } else {
        Box::new(BufReader::new(f))
    };
    let norm = |c: &str| c.strip_prefix("chr").unwrap_or(c).to_string();
    let mut want: HashMap<String, Vec<(u64, u64)>> = HashMap::new();
    for (c, s, e) in spans {
        want.entry(norm(c)).or_default().push((*s, *e));
    }
    let hit = |c: &str, s: u64, e: u64| {
        want.get(&norm(c))
            .is_some_and(|v| v.iter().any(|(a, b)| s < *b && *a < e))
    };

    let mut tx: HashMap<String, Transcript> = HashMap::new();
    for line in reader.lines() {
        let line = line?;
        if line.starts_with('#') {
            continue;
        }
        let f: Vec<&str> = line.splitn(9, '\t').collect();
        if f.len() < 9 || !(f[2] == "transcript" || f[2] == "exon") {
            continue;
        }
        let (Ok(s), Ok(e)) = (f[3].parse::<u64>(), f[4].parse::<u64>()) else {
            continue;
        };
        let (s, e) = (s.saturating_sub(1), e);
        let Some(tid) = attr(f[8], "transcript_id") else {
            continue;
        };
        // exons outside the window still belong to a transcript crossing it
        // (GTFs list the transcript line before its exons)
        let keep = hit(f[0], s, e) || (f[2] == "exon" && tx.contains_key(tid));
        if !keep {
            continue;
        }
        let t = tx.entry(tid.to_string()).or_insert_with(|| Transcript {
            gene: attr(f[8], "gene_name")
                .or_else(|| attr(f[8], "gene_id"))
                .unwrap_or("?")
                .to_string(),
            chrom: f[0].to_string(),
            strand: f[6].chars().next().unwrap_or('.'),
            start: s,
            end: e,
            exons: Vec::new(),
            canonical: false,
        });
        if f[2] == "transcript" {
            t.start = s;
            t.end = e;
            t.canonical |=
                f[8].contains("\"MANE_Select\"") || f[8].contains("\"Ensembl_canonical\"");
        } else {
            t.exons.push((s, e));
            t.start = t.start.min(s);
            t.end = t.end.max(e);
        }
    }
    // one transcript per gene
    let mut best: HashMap<(String, String), Transcript> = HashMap::new();
    for t in tx.into_values() {
        let key = (t.chrom.clone(), t.gene.clone());
        let better = |a: &Transcript, b: &Transcript| {
            (a.canonical, a.end - a.start) > (b.canonical, b.end - b.start)
        };
        match best.get(&key) {
            Some(b) if !better(&t, b) => {}
            _ => {
                best.insert(key, t);
            }
        }
    }
    let mut g = Genes::default();
    for mut t in best.into_values() {
        t.exons.sort_unstable();
        g.by_chrom.entry(t.chrom.clone()).or_default().push(t);
    }
    for v in g.by_chrom.values_mut() {
        v.sort_by_key(|t| t.start);
    }
    Ok(g)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn picks_canonical_transcript() {
        let dir = std::env::temp_dir().join(format!("hapmeth_gtf_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let p = dir.join("t.gtf");
        std::fs::write(
            &p,
            "chrX\tx\ttranscript\t100\t5000\t.\t-\t.\tgene_id \"G\"; transcript_id \"T1\"; gene_name \"MECP2\"; tag \"basic\";\n\
             chrX\tx\texon\t100\t300\t.\t-\t.\tgene_id \"G\"; transcript_id \"T1\"; gene_name \"MECP2\";\n\
             chrX\tx\ttranscript\t100\t2000\t.\t-\t.\tgene_id \"G\"; transcript_id \"T2\"; gene_name \"MECP2\"; tag \"Ensembl_canonical\";\n\
             chrX\tx\texon\t100\t300\t.\t-\t.\tgene_id \"G\"; transcript_id \"T2\"; gene_name \"MECP2\";\n\
             chrX\tx\texon\t1800\t2000\t.\t-\t.\tgene_id \"G\"; transcript_id \"T2\"; gene_name \"MECP2\";\n\
             chr1\tx\ttranscript\t100\t2000\t.\t+\t.\tgene_id \"H\"; transcript_id \"T3\"; gene_name \"FAR\";\n",
        )
        .unwrap();
        let g = load(&p, &[("chrX".into(), 1500, 1600)]).unwrap();
        let t = g.overlapping("chrX", 0, 10_000);
        assert_eq!(t.len(), 1);
        assert_eq!(
            (t[0].gene.as_str(), t[0].end, t[0].exons.len()),
            ("MECP2", 2000, 2)
        );
        assert!(g.overlapping("chr1", 0, 10_000).is_empty());
        std::fs::remove_dir_all(dir).ok();
    }
}
