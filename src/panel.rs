//! Locus panels: the three built in (imprinting, X-inactivation, QC controls,
//! for hg38 and hs1) and any `--bed` files, all read by the same BED parser.
//!
//! BED3+ is enough (name defaults to chrom_start_end). A `#chrom\tstart\tend...`
//! header line naming extra columns lets a BED carry per-locus metadata:
//! gene, disease, origin, expected_lo, expected_hi (percent 5mC) and group
//! (output subfolder). The files in bed/ are written that way.

use std::{fmt, path::Path};

use anyhow::{bail, Context, Result};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Build {
    Hg38,
    Hs1,
}

impl fmt::Display for Build {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(match self {
            Build::Hg38 => "hg38",
            Build::Hs1 => "hs1",
        })
    }
}

impl Build {
    /// Recognise the build from chr1's length in the BAM header.
    pub fn from_chr1_len(len: usize) -> Option<Build> {
        match len {
            248_956_422 => Some(Build::Hg38),
            248_387_328 => Some(Build::Hs1),
            _ => None,
        }
    }
}

/// Built-in panels, in plotting order.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Builtin {
    Imprinting,
    Xci,
    Qc,
}

pub const BUILTINS: [Builtin; 3] = [Builtin::Imprinting, Builtin::Xci, Builtin::Qc];

impl Builtin {
    pub fn key(self) -> &'static str {
        match self {
            Builtin::Imprinting => "imprinting",
            Builtin::Xci => "xci",
            Builtin::Qc => "qc",
        }
    }

    /// Output subfolder / overview section title.
    pub fn title(self) -> &'static str {
        match self {
            Builtin::Imprinting => "Imprinting",
            Builtin::Xci => "X-inactivation",
            Builtin::Qc => "Controls",
        }
    }

    pub fn parse(s: &str) -> Option<Builtin> {
        BUILTINS.into_iter().find(|b| b.key() == s)
    }

    pub fn bed(self, build: Build) -> (&'static str, &'static str) {
        match (self, build) {
            (Builtin::Imprinting, Build::Hg38) => (
                "imprinted_germline_dmrs.hg38.bed",
                include_str!("../bed/imprinted_germline_dmrs.hg38.bed"),
            ),
            (Builtin::Imprinting, Build::Hs1) => (
                "imprinted_germline_dmrs.hs1.bed",
                include_str!("../bed/imprinted_germline_dmrs.hs1.bed"),
            ),
            (Builtin::Xci, Build::Hg38) => (
                "xci_panel.hg38.bed",
                include_str!("../bed/xci_panel.hg38.bed"),
            ),
            (Builtin::Xci, Build::Hs1) => (
                "xci_panel.hs1.bed",
                include_str!("../bed/xci_panel.hs1.bed"),
            ),
            (Builtin::Qc, Build::Hg38) => (
                "qc_controls.hg38.bed",
                include_str!("../bed/qc_controls.hg38.bed"),
            ),
            (Builtin::Qc, Build::Hs1) => (
                "qc_controls.hs1.bed",
                include_str!("../bed/qc_controls.hs1.bed"),
            ),
        }
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct Locus {
    pub name: String,
    pub chrom: String,
    /// 0-based, half-open.
    pub start: u64,
    pub end: u64,
    pub gene: String,
    pub disease: String,
    /// Which parental allele the atlas reports as methylated (imprinted DMRs).
    pub origin: String,
    /// Expected 5mC percent of the less / more methylated haplotype.
    pub expected: Option<(f64, f64)>,
    /// Output subfolder and overview section.
    pub group: String,
    /// Where the locus came from: built-in panel key or the --bed file name.
    pub source: String,
}

impl Locus {
    pub fn coords(&self) -> String {
        format!("{}:{}-{}", self.chrom, self.start + 1, self.end)
    }
}

fn dot(s: &str) -> String {
    let s = s.trim();
    if s.is_empty() {
        ".".into()
    } else {
        s.into()
    }
}

fn number(s: &str) -> Option<f64> {
    s.trim().parse::<f64>().ok().filter(|v| v.is_finite())
}

/// Parse BED text. `default_group` is used when a row has no group column.
pub fn parse_bed(text: &str, source: &str, default_group: &str) -> Result<Vec<Locus>> {
    let mut cols: Option<Vec<String>> = None;
    let mut out = Vec::new();
    for (i, line) in text.lines().enumerate() {
        let line = line.trim_end_matches('\r');
        if line.trim().is_empty() || line.starts_with("track") || line.starts_with("browser") {
            continue;
        }
        if let Some(h) = line.strip_prefix('#') {
            let f: Vec<String> = h
                .split('\t')
                .map(|s| s.trim().to_ascii_lowercase())
                .collect();
            if f.first().map(String::as_str) == Some("chrom") {
                cols = Some(f);
            }
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 3 {
            bail!(
                "{source} line {}: expected at least 3 tab-separated columns",
                i + 1
            );
        }
        let start: u64 = f[1]
            .trim()
            .parse()
            .with_context(|| format!("{source} line {}: bad start", i + 1))?;
        let end: u64 = f[2]
            .trim()
            .parse()
            .with_context(|| format!("{source} line {}: bad end", i + 1))?;
        if end <= start {
            bail!("{source} line {}: end <= start", i + 1);
        }
        let get = |name: &str| -> Option<&str> {
            let idx = cols.as_ref()?.iter().position(|c| c == name)?;
            f.get(idx)
                .copied()
                .filter(|v| !v.trim().is_empty() && v.trim() != ".")
        };
        let chrom = f[0].trim().to_string();
        let name = f
            .get(3)
            .map(|s| s.trim())
            .filter(|s| !s.is_empty() && *s != ".")
            .map(str::to_string)
            .unwrap_or_else(|| format!("{chrom}_{start}_{end}"));
        let expected = match (
            get("expected_lo").and_then(number),
            get("expected_hi").and_then(number),
        ) {
            (Some(lo), Some(hi)) => Some((lo.min(hi), lo.max(hi))),
            _ => None,
        };
        out.push(Locus {
            gene: dot(get("gene").unwrap_or(".")),
            disease: dot(get("disease").unwrap_or(".")),
            origin: dot(get("origin").unwrap_or(".")),
            expected,
            group: get("group")
                .map(str::to_string)
                .unwrap_or_else(|| default_group.to_string()),
            source: source.to_string(),
            name,
            chrom,
            start,
            end,
        });
    }
    Ok(out)
}

/// Built-in panel rows. Imprinting keeps only disease-anchored DMRs unless
/// `all_imprinted`; X-inactivation only makes sense with two X haplotypes.
pub fn builtin_loci(panel: Builtin, build: Build, all_imprinted: bool) -> Vec<Locus> {
    let (file, text) = panel.bed(build);
    let rows = parse_bed(text, panel.key(), panel.title())
        .unwrap_or_else(|e| panic!("embedded {file}: {e}"));
    rows.into_iter()
        .filter(|l| panel != Builtin::Imprinting || all_imprinted || l.disease != ".")
        .collect()
}

/// A `--bed` file; loci without a group column go under the file's stem
/// (with any .hg38/.hs1 suffix removed).
pub fn bed_file_loci(path: &Path) -> Result<Vec<Locus>> {
    let text = std::fs::read_to_string(path)
        .with_context(|| format!("reading --bed {}", path.display()))?;
    let file = path
        .file_name()
        .map(|s| s.to_string_lossy().into_owned())
        .unwrap_or_default();
    let mut stem = file
        .trim_end_matches(".gz")
        .trim_end_matches(".bed")
        .to_string();
    for sfx in [".hg38", ".hs1", ".grch38", ".chm13", ".t2t"] {
        if let Some(s) = stem.strip_suffix(sfx) {
            stem = s.to_string();
        }
    }
    let rows = parse_bed(&text, &file, &stem)?;
    if rows.is_empty() {
        bail!("--bed {}: no regions", path.display());
    }
    Ok(rows)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn embedded_panels_parse_and_agree_between_builds() {
        for p in BUILTINS {
            let a = builtin_loci(p, Build::Hg38, true);
            let b = builtin_loci(p, Build::Hs1, true);
            assert!(!a.is_empty());
            let names = |v: &[Locus]| {
                let mut n: Vec<_> = v.iter().map(|l| l.name.clone()).collect();
                n.sort();
                n
            };
            assert_eq!(names(&a), names(&b), "{p:?}: hg38 and hs1 panels differ");
        }
        assert_eq!(
            builtin_loci(Builtin::Imprinting, Build::Hg38, true).len(),
            79
        );
        assert_eq!(builtin_loci(Builtin::Xci, Build::Hs1, false).len(), 3);
        assert_eq!(builtin_loci(Builtin::Qc, Build::Hg38, false).len(), 4);
    }

    #[test]
    fn default_imprinting_is_disease_anchored() {
        let d = builtin_loci(Builtin::Imprinting, Build::Hg38, false);
        assert!(d.iter().all(|l| l.disease != "."));
        for must in [
            "H19", "KvDMR1", "SNURF", "GNAS_XL", "MEST", "PLAGL1", "GRB10", "PEG3",
        ] {
            assert!(
                d.iter().any(|l| l.name == must),
                "{must} missing from default imprinting set"
            );
        }
    }

    #[test]
    fn metadata_columns_are_read() {
        let x = builtin_loci(Builtin::Xci, Build::Hg38, false);
        let mecp2 = x.iter().find(|l| l.name == "MECP2_promoter_CGI").unwrap();
        assert_eq!(
            (mecp2.chrom.as_str(), mecp2.start, mecp2.end),
            ("chrX", 154_097_110, 154_098_015)
        );
        assert_eq!(mecp2.gene, "MECP2");
        assert_eq!(mecp2.expected, Some((50.0, 50.0)));
        assert_eq!(mecp2.group, "X-inactivation");
    }

    #[test]
    fn plain_bed() {
        let l = parse_bed("chr1\t10\t20\nchr2\t5\t9\tfoo\t0\t+\n", "x.bed", "x").unwrap();
        assert_eq!(l[0].name, "chr1_10_20");
        assert_eq!(l[1].name, "foo");
        assert_eq!(l[1].group, "x");
        assert_eq!(l[1].expected, None);
        assert!(parse_bed("chr1\t20\t10\n", "x.bed", "x").is_err());
    }
}
