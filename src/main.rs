//! hapmeth -- per-haplotype 5mC locus plots from a haplotagged long-read modBAM.
//!
//! One self-contained binary: reads MM/ML itself (no modkit), draws the figures
//! itself (no methylartist), and carries the imprinting / X-inactivation / QC
//! panels for hg38 and hs1. See README.md.

mod gtf;
mod meth;
mod modbam;
mod panel;
mod plot;
mod reads;
mod render;

use std::{
    collections::{BTreeMap, HashMap},
    fs,
    io::Write,
    path::{Path, PathBuf},
    time::Instant,
};

use anyhow::{bail, Context, Result};
use clap::{Parser, ValueEnum};
use noodles::fasta;
use rayon::prelude::*;

use crate::{
    meth::{Group, LocusMeth},
    panel::{Build, Builtin, Locus, BUILTINS},
    reads::{Alignments, Filters, Read, RefWindow, Skipped},
};

#[derive(Clone, Copy, Debug, PartialEq, Eq, ValueEnum)]
enum BuildArg {
    Auto,
    Hg38,
    Hs1,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, ValueEnum)]
enum SexArg {
    Auto,
    Female,
    Male,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, ValueEnum)]
pub enum Format {
    Png,
    Svg,
    Both,
}

const AFTER_HELP: &str = "\
Outputs (in --out-dir):
  <sample>.overview.png          every locus on one figure: HP1 vs HP2 5mC (for a QC report)
  <sample>.loci.tsv              one row per locus: reads, CpGs and 5mC % per haplotype
  plots/<panel>/<locus>.png      per-locus figure: reads by haplotype + 5mC along the locus

Examples:
  hapmeth -b sample.haplotagged.bam -r hs1.fa
  hapmeth -b s.bam -r GRCh38.fa --all-imprinted --bed focal_dx.hg38.bed
  hapmeth -b s.bam -r GRCh38.fa --panels none --bed my_regions.bed
  hapmeth --export-panels panels/     (write the built-in BEDs and exit)";

/// Per-haplotype 5mC locus plots from a haplotagged long-read modBAM.
///
/// Reads the MM/ML tags directly (dorado's implicit '.' calls count as
/// unmethylated), splits reads by HP tag, and plots each locus of the built-in
/// imprinting, X-inactivation and QC-control panels (hg38 or hs1, picked from
/// the BAM header) plus any --bed regions.
#[derive(Parser, Debug)]
#[command(name = "hapmeth", version, after_help = AFTER_HELP)]
struct Args {
    /// Haplotagged modBAM or CRAM (HP + MM/ML tags), with .bai/.csi/.crai
    #[arg(
        short,
        long,
        required_unless_present = "export_panels",
        help_heading = "Input"
    )]
    bam: Option<PathBuf>,

    /// Reference FASTA the reads are aligned to (.fai; bgzipped also needs .gzi)
    #[arg(
        short,
        long = "ref",
        value_name = "FASTA",
        required_unless_present = "export_panels",
        help_heading = "Input"
    )]
    reference: Option<PathBuf>,

    /// Output directory [default: ./<sample>]
    #[arg(short, long, help_heading = "Output")]
    out_dir: Option<PathBuf>,

    /// Sample label for titles and file names [default: BAM file name]
    #[arg(short, long, help_heading = "Output")]
    sample: Option<String>,

    /// Figure format
    #[arg(long, value_enum, default_value = "png", help_heading = "Output")]
    format: Format,

    /// Built-in panels to plot: imprinting, xci, qc (comma-separated), or none
    #[arg(long, default_value = "imprinting,xci,qc", help_heading = "Loci")]
    panels: String,

    /// All 79 imprinted germline DMRs, not only the disease-anchored ones
    #[arg(long, help_heading = "Loci")]
    all_imprinted: bool,

    /// Extra regions: BED3+ (name in column 4); repeat for several files
    #[arg(long, value_name = "BED", help_heading = "Loci")]
    bed: Vec<PathBuf>,

    /// Only plot loci with this name; repeat for several
    #[arg(long, value_name = "NAME", help_heading = "Loci")]
    locus: Vec<String>,

    /// Genome build of the built-in panels [auto: from chr1's length]
    #[arg(long, value_enum, default_value = "auto", help_heading = "Loci")]
    build: BuildArg,

    /// Sample sex; X-inactivation loci are plotted for females only [auto: chrX/chrY depth from the BAM index]
    #[arg(long, value_enum, default_value = "auto", help_heading = "Loci")]
    sex: SexArg,

    /// Flank (bp) plotted each side of a locus
    #[arg(long, default_value_t = 1500, help_heading = "Tuning")]
    pad: u64,

    /// Minimum mapping quality (lower to 1 for segmental-duplication loci)
    #[arg(long, default_value_t = 10, help_heading = "Tuning")]
    min_mapq: u8,

    /// Also use supplementary alignments
    #[arg(long, help_heading = "Tuning")]
    include_supplementary: bool,

    /// Drop calls whose probability is below this (modkit --filter-threshold)
    #[arg(long, value_name = "P", help_heading = "Tuning")]
    filter_threshold: Option<f32>,

    /// Without --filter-threshold, drop the least confident fraction of calls (modkit's default rule)
    #[arg(long, default_value_t = 0.1, value_name = "Q", help_heading = "Tuning")]
    filter_percentile: f64,

    /// Reads drawn per haplotype in a locus figure (all reads are counted)
    #[arg(long, default_value_t = 60, help_heading = "Tuning")]
    max_reads: usize,

    /// Also write <sample>.cpgs.tsv: per-CpG counts per haplotype in every plotted window
    #[arg(long, help_heading = "Output")]
    sites: bool,

    /// Gene annotation (GTF, plain or gzipped) for a gene track
    #[arg(long, help_heading = "Tuning")]
    gtf: Option<PathBuf>,

    /// Threads [default: all cores]
    #[arg(
        short = 'j',
        long,
        default_value_t = 0,
        hide_default_value = true,
        help_heading = "Tuning"
    )]
    threads: usize,

    /// Write the built-in panel BEDs (hg38 + hs1) to this directory and exit
    #[arg(long, value_name = "DIR", help_heading = "Other")]
    export_panels: Option<PathBuf>,
}

/// A locus with everything needed to plot it.
pub struct Job {
    pub locus: Locus,
    /// Contig name as spelled in the BAM.
    pub contig: String,
    /// Plotted window, 0-based half-open.
    pub win: (u64, u64),
    pub refw: RefWindow,
    pub out_stem: PathBuf,
}

pub struct Done {
    pub job: Job,
    pub meth: LocusMeth,
    /// Reads in the window that were not used, by reason.
    pub skipped: Skipped,
    pub plot: Result<Vec<PathBuf>, String>,
}

/// Why a locus has fewer reads than expected, for the table and the figure.
pub fn skip_note(s: &Skipped, min_mapq: u8) -> String {
    let mut n = Vec::new();
    if s.mapq > 0 {
        n.push(format!("{} reads below --min-mapq {min_mapq}", s.mapq));
    }
    if s.no_mod_tags > 0 {
        n.push(format!("{} reads without MM/ML", s.no_mod_tags));
    }
    if s.bad_mod_tags > 0 {
        n.push(format!("{} reads with unusable MM/ML", s.bad_mod_tags));
    }
    n.join("; ")
}

fn default_sample(bam: &Path) -> String {
    let n = bam
        .file_name()
        .map(|s| s.to_string_lossy().into_owned())
        .unwrap_or_default();
    for sfx in [
        ".haplotagged.bam",
        ".haplotagged.cram",
        ".haplotag.bam",
        ".bam",
        ".cram",
    ] {
        if let Some(s) = n.strip_suffix(sfx) {
            return s.to_string();
        }
    }
    n
}

/// The BAM's spelling of a contig name: exact, else with/without "chr", chrM/MT.
fn resolve_contig(names: &HashMap<String, u64>, chrom: &str) -> Option<String> {
    let alt = match chrom {
        "chrM" => "MT".to_string(),
        "MT" => "chrM".to_string(),
        c => c
            .strip_prefix("chr")
            .map(str::to_string)
            .unwrap_or_else(|| format!("chr{c}")),
    };
    [chrom.to_string(), alt]
        .into_iter()
        .find(|c| names.contains_key(c))
}

fn safe(s: &str) -> String {
    s.chars()
        .map(|c| {
            if c.is_ascii_alphanumeric() || "._-".contains(c) {
                c
            } else {
                '_'
            }
        })
        .collect()
}

/// chrX/chrY depth relative to chr20 from the BAM index (what samtools idxstats reports).
fn infer_sex(aln: &Alignments) -> (Option<plot::Sex>, String) {
    let Some(counts) = aln.mapped_counts() else {
        return (None, "no mapped-read counts in the index (CRAM?)".into());
    };
    let depth = |names: [&str; 2]| {
        counts
            .iter()
            .find(|(n, _, _)| names.contains(&n.as_str()))
            .and_then(|(_, len, m)| (*len > 0 && *m > 0).then(|| *m as f64 / *len as f64))
    };
    let (Some(auto), Some(x)) = (depth(["chr20", "20"]), depth(["chrX", "X"])) else {
        return (None, "chr20/chrX absent or empty in the index".into());
    };
    let y = depth(["chrY", "Y"]).unwrap_or(0.0);
    let (xr, yr) = (x / auto, y / auto);
    let desc = format!("chrX/chr20 {xr:.2}, chrY/chr20 {yr:.2}");
    match sex_from_ratios(xr, yr) {
        Some(s) => (Some(s), desc),
        None => (None, format!("ambiguous ({desc})")),
    }
}

/// chrY present and chrX halved -> male; no chrY and chrX diploid -> female;
/// anything else (XXY, low coverage, contamination) is left undecided.
fn sex_from_ratios(x_ratio: f64, y_ratio: f64) -> Option<plot::Sex> {
    if y_ratio >= 0.15 && x_ratio <= 0.75 {
        Some(plot::Sex::Male)
    } else if y_ratio < 0.05 && x_ratio >= 0.80 {
        Some(plot::Sex::Female)
    } else {
        None
    }
}

fn export_panels(dir: &Path) -> Result<()> {
    fs::create_dir_all(dir)?;
    for p in BUILTINS {
        for b in [Build::Hg38, Build::Hs1] {
            let (name, text) = p.bed(b);
            fs::write(dir.join(name), text)?;
            eprintln!("wrote {}", dir.join(name).display());
        }
    }
    Ok(())
}

fn main() {
    if let Err(e) = run() {
        eprintln!("error: {e:#}");
        std::process::exit(1);
    }
}

fn run() -> Result<()> {
    let a = Args::parse();
    if let Some(dir) = &a.export_panels {
        return export_panels(dir);
    }
    let (bam, reference) = (a.bam.clone().unwrap(), a.reference.clone().unwrap());
    let t0 = Instant::now();
    if a.threads > 0 {
        rayon::ThreadPoolBuilder::new()
            .num_threads(a.threads)
            .build_global()?;
    }
    if !(0.0..1.0).contains(&a.filter_percentile) {
        bail!("--filter-percentile must be in [0, 1)");
    }

    let aln = Alignments::open(&bam, &reference)?;
    let contigs: HashMap<String, u64> = aln
        .header
        .reference_sequences()
        .iter()
        .map(|(n, rs)| (n.to_string(), usize::from(rs.length()) as u64))
        .collect();
    let chr1 = resolve_contig(&contigs, "chr1").and_then(|c| contigs.get(&c).copied());
    let build = match a.build {
        BuildArg::Hg38 => Build::Hg38,
        BuildArg::Hs1 => Build::Hs1,
        BuildArg::Auto => chr1
            .and_then(|l| Build::from_chr1_len(l as usize))
            .with_context(|| {
                format!(
                    "cannot tell hg38 from hs1: chr1 is {} in the BAM header; pass --build",
                    chr1.map_or("absent".into(), |l| format!("{l} bp"))
                )
            })?,
    };
    if let (BuildArg::Hg38 | BuildArg::Hs1, Some(l)) = (a.build, chr1) {
        if Build::from_chr1_len(l as usize).is_some_and(|b| b != build) {
            eprintln!("warning: --build {build} but the BAM's chr1 length says otherwise");
        }
    }

    let sample = a.sample.clone().unwrap_or_else(|| default_sample(&bam));
    let out = a.out_dir.clone().unwrap_or_else(|| PathBuf::from(&sample));
    let (sex, sex_note) = match a.sex {
        SexArg::Female => (Some(plot::Sex::Female), "given".to_string()),
        SexArg::Male => (Some(plot::Sex::Male), "given".to_string()),
        SexArg::Auto => infer_sex(&aln),
    };
    eprintln!(
        "{sample}: build {build}, sex {} ({sex_note})",
        sex.map_or("unknown", plot::Sex::label)
    );

    // ---- loci -------------------------------------------------------------
    let mut loci: Vec<Locus> = Vec::new();
    let wanted: Vec<&str> = a
        .panels
        .split(',')
        .map(str::trim)
        .filter(|s| !s.is_empty() && *s != "none")
        .collect();
    for key in &wanted {
        let p = Builtin::parse(key)
            .with_context(|| format!("unknown panel '{key}' (imprinting, xci, qc, none)"))?;
        let rows = panel::builtin_loci(p, build, a.all_imprinted);
        if p == Builtin::Xci && sex != Some(plot::Sex::Female) {
            let named: Vec<Locus> = rows
                .into_iter()
                .filter(|l| a.locus.contains(&l.name))
                .collect();
            if named.is_empty() {
                eprintln!(
                    "note: X-inactivation loci skipped (sex not female; --sex female to force)"
                );
            } else {
                eprintln!("warning: plotting X-inactivation loci named by --locus although sex is not female");
            }
            loci.extend(named);
            continue;
        }
        loci.extend(rows);
    }
    for b in &a.bed {
        loci.extend(panel::bed_file_loci(b)?);
    }
    if !a.locus.is_empty() {
        loci.retain(|l| a.locus.contains(&l.name));
        for n in &a.locus {
            if !loci.iter().any(|l| &l.name == n) {
                eprintln!("warning: --locus {n} not found in the selected panels");
            }
        }
    }
    if loci.is_empty() {
        bail!("no loci selected (check --panels / --bed / --locus)");
    }

    // windows + reference, in BAM header order
    let order: HashMap<String, usize> = aln
        .header
        .reference_sequences()
        .keys()
        .enumerate()
        .map(|(i, n)| (n.to_string(), i))
        .collect();
    let mut fa = fasta::io::indexed_reader::Builder::default()
        .build_from_path(&reference)
        .with_context(|| {
            format!(
                "opening {} (needs .fai; bgzipped also .gzi)",
                reference.display()
            )
        })?;
    let mut jobs: Vec<Job> = Vec::new();
    for l in loci {
        let Some(contig) = resolve_contig(&contigs, &l.chrom) else {
            eprintln!(
                "warning: {} ({}) skipped: {} not in the BAM header",
                l.name,
                l.coords(),
                l.chrom
            );
            continue;
        };
        let len = contigs[&contig];
        if l.start >= len {
            eprintln!(
                "warning: {} ({}) skipped: beyond the end of {contig}",
                l.name,
                l.coords()
            );
            continue;
        }
        let win = (l.start.saturating_sub(a.pad), (l.end + a.pad).min(len));
        let refw = reads::fetch_ref(&mut fa, &contig, win.0, win.1, len)
            .with_context(|| format!("does --ref match the BAM? ({})", l.name))?;
        let stem = out
            .join("plots")
            .join(safe_path(&l.group))
            .join(safe(&l.name));
        jobs.push(Job {
            locus: l,
            contig,
            win,
            refw,
            out_stem: stem,
        });
    }
    // same name in the same folder (e.g. two IG-DMR rows): add coordinates to all of them
    let mut n_stem: HashMap<PathBuf, usize> = HashMap::new();
    for j in &jobs {
        *n_stem.entry(j.out_stem.clone()).or_default() += 1;
    }
    for j in &mut jobs {
        if n_stem[&j.out_stem] > 1 {
            let l = &j.locus;
            j.out_stem = j.out_stem.with_file_name(format!(
                "{}__{}_{}_{}",
                safe(&l.name),
                l.chrom,
                l.start,
                l.end
            ));
        }
    }
    jobs.sort_by_key(|j| (order.get(&j.contig).copied().unwrap_or(usize::MAX), j.win.0));
    let mut n_by_group: BTreeMap<&str, usize> = BTreeMap::new();
    for j in &jobs {
        *n_by_group.entry(&j.locus.group).or_default() += 1;
    }
    eprintln!(
        "{} loci: {}",
        jobs.len(),
        n_by_group
            .iter()
            .map(|(g, n)| format!("{g} {n}"))
            .collect::<Vec<_>>()
            .join(", ")
    );

    // ---- reads (per chromosome, so a CRAM's reference cache stays small) ----
    let filters = Filters {
        min_mapq: a.min_mapq,
        include_supplementary: a.include_supplementary,
        collect_conf: a.filter_threshold.is_none(),
    };
    let mut fetched: Vec<(Vec<Read>, Skipped)> = Vec::with_capacity(jobs.len());
    let mut i = 0;
    while i < jobs.len() {
        let chrom = jobs[i].contig.clone();
        let n = jobs[i..].iter().take_while(|j| j.contig == chrom).count();
        let part: Vec<Result<(Vec<Read>, Skipped)>> = jobs[i..i + n]
            .par_iter()
            .map_init(
                || aln.reader(),
                |r, j| match r {
                    Ok(r) => r.fetch(&j.contig, j.win.0, j.win.1, &j.refw, filters),
                    Err(e) => Err(anyhow::anyhow!("opening {}: {e}", bam.display())),
                },
            )
            .collect();
        for (j, res) in jobs[i..i + n].iter().zip(part) {
            fetched.push(
                res.with_context(|| format!("reading {} ({})", j.locus.name, j.locus.coords()))?,
            );
        }
        aln.clear_reference_cache();
        i += n;
    }
    let mut skipped = Skipped::default();
    for (_, s) in &fetched {
        skipped.add(s);
    }

    // ---- call-confidence threshold -----------------------------------------
    let threshold = match a.filter_threshold {
        Some(t) => t,
        None => {
            // each read once, even when it spans several loci
            let mut seen = std::collections::HashSet::new();
            let conf: Vec<f32> = fetched
                .iter()
                .flat_map(|(r, _)| r.iter())
                .filter(|r| seen.insert(r.name.as_str()))
                .flat_map(|r| r.conf.iter().copied())
                .collect();
            let n = conf.len();
            match meth::estimate_threshold(conf, a.filter_percentile) {
                Some(t) => {
                    eprintln!(
                        "call-confidence threshold {t:.3} ({:.0}th percentile of {n} CpG calls on {} reads)",
                        a.filter_percentile * 100.0,
                        seen.len()
                    );
                    t
                }
                None => 0.0,
            }
        }
    };

    // ---- aggregate + plot ---------------------------------------------------
    let genes = match &a.gtf {
        Some(p) => {
            let spans: Vec<(String, u64, u64)> = jobs
                .iter()
                .map(|j| (j.locus.chrom.clone(), j.win.0, j.win.1))
                .collect();
            gtf::load(p, &spans).with_context(|| format!("reading --gtf {}", p.display()))?
        }
        None => gtf::Genes::default(),
    };
    let ctx = plot::Context {
        sample: &sample,
        build,
        threshold,
        min_mapq: a.min_mapq,
        max_reads: a.max_reads,
        genes: &genes,
    };
    let format = a.format;
    let done: Vec<Done> = jobs
        .into_par_iter()
        .zip(fetched.into_par_iter())
        .map(|(job, (reads, skipped))| {
            let cpgs = job.refw.cpgs(job.win.0, job.win.1);
            let meth = meth::aggregate(&reads, cpgs, (job.locus.start, job.locus.end), threshold);
            let plot = plot::locus_svg(&ctx, &job, &reads, &meth, &skipped)
                .and_then(|svg| render::write(&svg, &job.out_stem, format))
                .map_err(|e| format!("{e:#}"));
            Done {
                job,
                meth,
                skipped,
                plot,
            }
        })
        .collect();

    fs::create_dir_all(&out).with_context(|| format!("creating {}", out.display()))?;
    let overview = render::write(
        &plot::overview_svg(&ctx, &done, sex),
        &out.join(format!("{}.overview", safe(&sample))),
        format,
    )
    .context("drawing the overview")?;
    let tsv = out.join(format!("{}.loci.tsv", safe(&sample)));
    write_manifest(&tsv, &sample, &out, &done, threshold, a.min_mapq)?;
    if a.sites {
        write_sites(
            &out.join(format!("{}.cpgs.tsv", safe(&sample))),
            &done,
            threshold,
        )?;
    }

    let ok = done.iter().filter(|d| d.plot.is_ok()).count();
    for d in done.iter().filter(|d| d.plot.is_err()) {
        eprintln!(
            "warning: {}: plot failed: {}",
            d.job.locus.name,
            d.plot.as_ref().unwrap_err()
        );
    }
    if skipped.bad_mod_tags > 0 {
        eprintln!(
            "note: {} reads had unusable MM/ML tags (malformed, or SEQ clipped after tagging)",
            skipped.bad_mod_tags
        );
    }
    if skipped.no_mod_tags > 0 {
        eprintln!("note: {} reads had no MM/ML tags", skipped.no_mod_tags);
    }
    eprintln!(
        "{sample}: {ok}/{} locus plots, overview {}, table {} ({:.1}s)",
        done.len(),
        overview
            .first()
            .map_or(String::new(), |p| p.display().to_string()),
        tsv.display(),
        t0.elapsed().as_secs_f64()
    );
    if ok == 0 {
        bail!("no locus plots were written");
    }
    Ok(())
}

/// Group names may contain '/' (MLID/<category>): keep it as a subfolder.
fn safe_path(group: &str) -> PathBuf {
    group
        .split('/')
        .filter(|s| !s.is_empty())
        .map(safe)
        .collect()
}

fn fmt_pct(v: Option<f64>) -> String {
    v.map_or(".".into(), |v| format!("{v:.1}"))
}

/// Per-CpG counts, one row per (locus, CpG, group) with any call; the same
/// quantities as a modkit bedMethyl row (strands combined).
fn write_sites(path: &Path, done: &[Done], threshold: f32) -> Result<()> {
    let mut w = std::io::BufWriter::new(
        fs::File::create(path).with_context(|| format!("creating {}", path.display()))?,
    );
    writeln!(
        w,
        "# hapmeth {}; call-confidence threshold {threshold:.3}; pos = 0-based C of the CpG",
        env!("CARGO_PKG_VERSION")
    )?;
    writeln!(
        w,
        "locus\tchrom\tpos\tgroup\tn_5mC\tn_5hmC\tn_C\tn_other\tn_fail\tpct_5mC"
    )?;
    for d in done {
        for g in meth::GROUPS {
            for (pos, s) in d.meth.cpgs.iter().zip(&d.meth.sites[g.index()]) {
                if s.valid() + s.fail == 0 {
                    continue;
                }
                writeln!(
                    w,
                    "{}\t{}\t{pos}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                    d.job.locus.name,
                    d.job.contig,
                    g.key(),
                    s.m,
                    s.h,
                    s.c,
                    s.other,
                    s.fail,
                    fmt_pct(s.pct())
                )?;
            }
        }
    }
    w.flush()?;
    Ok(())
}

fn write_manifest(
    path: &Path,
    sample: &str,
    out: &Path,
    done: &[Done],
    threshold: f32,
    min_mapq: u8,
) -> Result<()> {
    let mut w = std::io::BufWriter::new(
        fs::File::create(path).with_context(|| format!("creating {}", path.display()))?,
    );
    writeln!(w, "# hapmeth {}; call-confidence threshold {threshold:.3}; 5mC % pooled over the core (un-padded) CpGs", env!("CARGO_PKG_VERSION"))?;
    let cols = [
        "sample",
        "name",
        "group",
        "gene",
        "disease",
        "chrom",
        "start",
        "end",
        "core_cpgs",
        "hp1_reads",
        "hp1_cpgs",
        "hp1_pct",
        "hp2_reads",
        "hp2_cpgs",
        "hp2_pct",
        "untagged_reads",
        "untagged_pct",
        "all_reads",
        "all_pct",
        "hp_diff",
        "expected_lo",
        "expected_hi",
        "phase_sets",
        "plot",
        "status",
        "note",
    ];
    writeln!(w, "{}", cols.join("\t"))?;
    for d in done {
        let l = &d.job.locus;
        let c = &d.meth.core;
        let (h1, h2, un, all) = (
            &c[Group::Hp1.index()],
            &c[Group::Hp2.index()],
            &c[Group::Unphased.index()],
            &c[Group::All.index()],
        );
        let diff = match (h1.pct(), h2.pct()) {
            (Some(x), Some(y)) => format!("{:.1}", (x - y).abs()),
            _ => ".".into(),
        };
        let (plot, status) = match &d.plot {
            Ok(p) => (
                p.first()
                    .map(|p| p.strip_prefix(out).unwrap_or(p).display().to_string())
                    .unwrap_or_default(),
                if all.reads == 0 {
                    "no_reads".to_string()
                } else {
                    "ok".to_string()
                },
            ),
            Err(e) => (
                ".".into(),
                format!("plot_failed: {}", e.replace(['\t', '\n'], " ")),
            ),
        };
        let (elo, ehi) = l.expected.map_or((".".into(), ".".into()), |(lo, hi)| {
            (format!("{lo:.1}"), format!("{hi:.1}"))
        });
        let row = [
            sample.to_string(),
            l.name.clone(),
            l.group.clone(),
            l.gene.clone(),
            l.disease.clone(),
            l.chrom.clone(),
            l.start.to_string(),
            l.end.to_string(),
            d.meth.core_cpgs.to_string(),
            h1.reads.to_string(),
            h1.cpgs_covered.to_string(),
            fmt_pct(h1.pct()),
            h2.reads.to_string(),
            h2.cpgs_covered.to_string(),
            fmt_pct(h2.pct()),
            un.reads.to_string(),
            fmt_pct(un.pct()),
            all.reads.to_string(),
            fmt_pct(all.pct()),
            diff,
            elo,
            ehi,
            d.meth.phase_sets.to_string(),
            plot,
            status,
            skip_note(&d.skipped, min_mapq),
        ];
        writeln!(w, "{}", row.join("\t"))?;
    }
    w.flush()?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sex_rule() {
        assert_eq!(sex_from_ratios(1.0, 0.001), Some(plot::Sex::Female));
        assert_eq!(sex_from_ratios(0.5, 0.45), Some(plot::Sex::Male));
        assert_eq!(sex_from_ratios(0.95, 0.4), None); // XXY
        assert_eq!(sex_from_ratios(0.6, 0.08), None);
    }

    #[test]
    fn contig_names() {
        let names: HashMap<String, u64> = [
            ("1".to_string(), 1),
            ("MT".to_string(), 1),
            ("chrX".to_string(), 1),
        ]
        .into();
        assert_eq!(resolve_contig(&names, "chr1").as_deref(), Some("1"));
        assert_eq!(resolve_contig(&names, "chrM").as_deref(), Some("MT"));
        assert_eq!(resolve_contig(&names, "X").as_deref(), Some("chrX"));
        assert_eq!(resolve_contig(&names, "chr2"), None);
    }

    #[test]
    fn sample_from_file_name() {
        assert_eq!(
            default_sample(Path::new("/x/HG002.haplotagged.bam")),
            "HG002"
        );
        assert_eq!(default_sample(Path::new("s1.cram")), "s1");
    }
}
