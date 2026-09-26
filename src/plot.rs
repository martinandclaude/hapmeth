//! The two figures, built as SVG (rasterised by render.rs).
//!
//! Colour has one job per encoding: haplotype identity is categorical (HP1
//! blue, HP2 orange -- validated colour-blind safe, all pairs); per-call 5mC is
//! a neutral ink ramp (dark = methylated, like a lollipop plot) so it never
//! competes with the haplotype hues; untagged / all reads are muted gray. The
//! loci.tsv table is the figures' table view.

use std::fmt::Write as _;

use anyhow::Result;

use crate::{
    gtf::Genes,
    meth::{self, Group, LocusMeth, Site},
    panel::Build,
    reads::{Read, Skipped},
    Done, Job,
};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Sex {
    Female,
    Male,
}

impl Sex {
    pub fn label(self) -> &'static str {
        match self {
            Sex::Female => "female",
            Sex::Male => "male",
        }
    }
}

pub struct Context<'a> {
    pub sample: &'a str,
    pub build: Build,
    pub threshold: f32,
    pub min_mapq: u8,
    pub max_reads: usize,
    pub genes: &'a Genes,
}

// ---- palette ----------------------------------------------------------------
const SURFACE: &str = "#ffffff";
const INK: &str = "#0b0b0b";
const INK_2: &str = "#52514e";
const MUTED: &str = "#898781";
const GRID: &str = "#e1e0d9";
const AXIS: &str = "#c3c2b7";
const CORE: &str = "#f1f0ec";
const READ_BAR: &str = "#ebeae4";
const HP1: &str = "#2a78d6";
const HP2: &str = "#eb6834";
/// p(5mC) 0 -> 1; validated as an ordinal ramp (light end 2.4:1 on white).
const METH_RAMP: [(u8, u8, u8); 5] = [
    (0xa8, 0xa7, 0x9f),
    (0x7e, 0x7d, 0x77),
    (0x52, 0x51, 0x4e),
    (0x2c, 0x2b, 0x29),
    (0x0b, 0x0b, 0x0b),
];

fn group_color(g: Group) -> &'static str {
    match g {
        Group::Hp1 => HP1,
        Group::Hp2 => HP2,
        Group::Unphased | Group::All => MUTED,
    }
}

fn meth_color(p: f32) -> String {
    let x = p.clamp(0.0, 1.0) * (METH_RAMP.len() - 1) as f32;
    let i = (x.floor() as usize).min(METH_RAMP.len() - 2);
    let t = x - i as f32;
    let (a, b) = (METH_RAMP[i], METH_RAMP[i + 1]);
    let mix = |u: u8, v: u8| (f32::from(u) + (f32::from(v) - f32::from(u)) * t).round() as u8;
    format!(
        "#{:02x}{:02x}{:02x}",
        mix(a.0, b.0),
        mix(a.1, b.1),
        mix(a.2, b.2)
    )
}

// ---- tiny SVG writer ---------------------------------------------------------
fn esc(s: &str) -> String {
    s.replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
        .replace('"', "&quot;")
}

/// Rough DejaVu Sans advance width; used only to keep labels inside their column.
fn text_w(s: &str, size: f64) -> f64 {
    s.chars()
        .map(|c| {
            if c.is_ascii_uppercase() || c == 'm' || c == 'w' {
                0.72
            } else if " .,:;'|il".contains(c) {
                0.32
            } else {
                0.6
            }
        })
        .sum::<f64>()
        * size
}

fn fit(s: &str, size: f64, max_w: f64) -> String {
    if text_w(s, size) <= max_w {
        return s.to_string();
    }
    let mut out = String::new();
    for c in s.chars() {
        if text_w(&format!("{out}{c}…"), size) > max_w {
            break;
        }
        out.push(c);
    }
    format!("{out}…")
}

struct Svg {
    body: String,
}

impl Svg {
    fn new() -> Svg {
        Svg {
            body: String::new(),
        }
    }
    fn rect(&mut self, x: f64, y: f64, w: f64, h: f64, fill: &str) {
        let _ = write!(
            self.body,
            r#"<rect x="{x:.2}" y="{y:.2}" width="{:.2}" height="{:.2}" fill="{fill}"/>"#,
            w.max(0.0),
            h.max(0.0)
        );
    }
    fn line(&mut self, x1: f64, y1: f64, x2: f64, y2: f64, stroke: &str, width: f64) {
        let _ = write!(
            self.body,
            r#"<line x1="{x1:.2}" y1="{y1:.2}" x2="{x2:.2}" y2="{y2:.2}" stroke="{stroke}" stroke-width="{width}"/>"#
        );
    }
    fn circle(&mut self, cx: f64, cy: f64, r: f64, fill: &str, extra: &str) {
        let _ = write!(
            self.body,
            r#"<circle cx="{cx:.2}" cy="{cy:.2}" r="{r}" fill="{fill}" {extra}/>"#
        );
    }
    #[allow(clippy::too_many_arguments)] // mirrors the SVG attributes one to one
    fn text(&mut self, x: f64, y: f64, s: &str, size: f64, fill: &str, anchor: &str, bold: bool) {
        let _ = write!(
            self.body,
            r#"<text x="{x:.2}" y="{y:.2}" font-size="{size}" fill="{fill}" text-anchor="{anchor}"{}>{}</text>"#,
            if bold { r#" font-weight="bold""# } else { "" },
            esc(s)
        );
    }
    fn raw(&mut self, s: &str) {
        self.body.push_str(s);
    }
    fn finish(self, w: f64, h: f64) -> String {
        format!(
            r#"<svg xmlns="http://www.w3.org/2000/svg" width="{w}" height="{h:.0}" viewBox="0 0 {w} {h:.0}" font-family="DejaVu Sans, sans-serif"><rect width="100%" height="100%" fill="{SURFACE}"/>{}</svg>"#,
            self.body
        )
    }
}

fn thousands(n: u64) -> String {
    let s = n.to_string();
    let mut out = String::new();
    for (i, c) in s.chars().enumerate() {
        if i > 0 && (s.len() - i).is_multiple_of(3) {
            out.push(',');
        }
        out.push(c);
    }
    out
}

/// Tick step giving roughly `n` ticks over `span`.
fn nice_step(span: f64, n: f64) -> f64 {
    let raw = span / n;
    let mag = 10f64.powf(raw.log10().floor());
    [1.0, 2.0, 2.5, 5.0, 10.0]
        .iter()
        .map(|m| m * mag)
        .find(|s| *s >= raw)
        .unwrap_or(10.0 * mag)
}

/// Placeholder labels from the source tables ("NA", "none", ...) say nothing.
fn real(v: &str) -> bool {
    !matches!(v, "." | "" | "NA" | "none" | "unknown")
}

fn reads_n(n: usize) -> String {
    if n == 1 {
        "1 read".into()
    } else {
        format!("{n} reads")
    }
}

fn pct(v: Option<f64>) -> String {
    v.map_or("n/a".into(), |v| format!("{v:.0}%"))
}

// ---- per-locus figure ---------------------------------------------------------
const W: f64 = 960.0;
const ML: f64 = 112.0;
const MR: f64 = 28.0;
const ROW: f64 = 4.0;
const BAR: f64 = 3.0;

/// Reads of one group to draw: at most `cap`, evenly spread over the
/// start-sorted list so the subset covers the window like the whole does.
fn pick(mut reads: Vec<&Read>, cap: usize) -> (Vec<&Read>, usize) {
    reads.sort_by_key(|r| (r.start, r.end));
    let n = reads.len();
    if n <= cap || cap == 0 {
        return (reads, n);
    }
    let picked = (0..cap).map(|i| reads[i * n / cap]).collect();
    (picked, n)
}

pub fn locus_svg(
    ctx: &Context,
    job: &Job,
    reads: &[Read],
    meth: &LocusMeth,
    skipped: &Skipped,
) -> Result<String> {
    let l = &job.locus;
    let (w0, w1) = job.win;
    let span = (w1 - w0).max(1) as f64;
    let pw = W - ML - MR;
    let x = |p: f64| ML + (p - w0 as f64) / span * pw;
    let mut s = Svg::new();
    let mut y = 0.0;

    // header
    let mut title = l.name.replace('_', " ");
    if real(&l.gene) && l.gene != l.name {
        title += &format!("  ·  {}", l.gene);
    }
    if real(&l.disease) {
        title += &format!("  ·  {}", l.disease.replace('_', " "));
    }
    s.text(
        16.0,
        24.0,
        &fit(&title, 15.0, W - 32.0),
        15.0,
        INK,
        "start",
        true,
    );
    let mut sub = format!(
        "{}:{}-{} ({})  ·  {}  ·  {}",
        l.chrom,
        thousands(l.start + 1),
        thousands(l.end),
        ctx.build,
        l.group,
        ctx.sample
    );
    if l.origin != "." {
        sub += &format!("  ·  {} allele methylated (atlas)", l.origin.to_lowercase());
    }
    s.text(
        16.0,
        42.0,
        &fit(&sub, 11.0, W - 32.0),
        11.0,
        INK_2,
        "start",
        false,
    );
    // stats line with colour keys (text stays in ink)
    s.text(16.0, 61.0, "Over the locus:", 11.0, INK_2, "start", false);
    let mut sx = 16.0 + text_w("Over the locus:", 11.0) + 10.0;
    let core = &meth.core;
    for g in [Group::Hp1, Group::Hp2, Group::Unphased] {
        let c = &core[g.index()];
        if g == Group::Unphased && c.reads == 0 {
            continue;
        }
        s.rect(sx, 52.0, 10.0, 10.0, group_color(g));
        let t = format!("{} {} · {}", g.label(), pct(c.pct()), reads_n(c.reads));
        s.text(sx + 14.0, 61.0, &t, 11.0, INK, "start", false);
        sx += 14.0 + text_w(&t, 11.0) + 16.0;
    }
    if let Some((lo, hi)) = l.expected {
        let t = if (hi - lo).abs() < 0.5 {
            format!("expected {lo:.0}%")
        } else {
            format!("expected {lo:.0}% / {hi:.0}%")
        };
        s.text(sx, 61.0, &t, 11.0, INK_2, "start", false);
    }
    y += 76.0;

    let top_of_tracks = y;
    // gene track: drawn into its own buffer so the locus shading can go underneath
    let tx: Vec<_> = ctx.genes.overlapping(&l.chrom, w0, w1);
    let mut genes = Svg::new();
    y += draw_genes(&mut genes, &tx, y, w0, w1, &x);

    // CpG rug
    s.text(ML - 8.0, y + 9.0, "CpG", 10.0, MUTED, "end", false);
    let rug_y = y;
    y += 14.0;

    // reads by haplotype
    let mut blocks: Vec<(Group, Vec<&Read>, usize)> = Vec::new();
    for g in [Group::Hp1, Group::Hp2, Group::Unphased] {
        let rs: Vec<&Read> = reads
            .iter()
            .filter(|r| {
                !r.calls.is_empty()
                    && matches!(
                        (g, r.hp),
                        (Group::Hp1, 1) | (Group::Hp2, 2) | (Group::Unphased, 0)
                    )
            })
            .collect();
        if rs.is_empty() {
            continue;
        }
        let cap = if g == Group::Unphased {
            ctx.max_reads.div_ceil(2)
        } else {
            ctx.max_reads
        };
        let (picked, n) = pick(rs, cap);
        blocks.push((g, picked, n));
    }
    // a block is at least tall enough for its two-line label
    let block_h = |n: usize| (n as f64 * ROW).max(26.0);
    let reads_h: f64 = blocks
        .iter()
        .map(|(_, r, _)| block_h(r.len()) + 8.0)
        .sum::<f64>()
        .max(24.0);
    let reads_y = y + 4.0;
    let meth_y = reads_y + reads_h + 16.0;
    let meth_h = 170.0;
    let axis_y = meth_y + meth_h;

    // core shading behind everything below the header
    let (cx0, cx1) = (x(l.start as f64), x(l.end as f64));
    s.rect(cx0, top_of_tracks, cx1 - cx0, axis_y - top_of_tracks, CORE);
    s.raw(&genes.body);

    // rug
    for p in &meth.cpgs {
        let px = x(*p as f64 + 0.5);
        s.line(px, rug_y + 2.0, px, rug_y + 10.0, MUTED, 0.6);
    }

    if blocks.is_empty() {
        let note = crate::skip_note(skipped, ctx.min_mapq);
        let msg = if note.is_empty() {
            "no reads with 5mC calls in this window".to_string()
        } else {
            format!("no usable reads in this window ({note})")
        };
        s.text(
            ML + pw / 2.0,
            reads_y + 16.0,
            &msg,
            11.0,
            MUTED,
            "middle",
            false,
        );
    }
    let mut by = reads_y;
    for (g, rs, n) in &blocks {
        let h = rs.len() as f64 * ROW;
        s.rect(ML - 6.0, by, 3.0, h - 1.0, group_color(*g));
        let count = if rs.len() < *n {
            format!("{} of {n} reads", rs.len())
        } else {
            reads_n(*n)
        };
        s.text(
            ML - 12.0,
            by + h / 2.0 - 1.0,
            g.label(),
            10.5,
            INK,
            "end",
            false,
        );
        s.text(
            ML - 12.0,
            by + h / 2.0 + 11.0,
            &count,
            9.5,
            MUTED,
            "end",
            false,
        );
        for (i, r) in rs.iter().enumerate() {
            let ry = by + i as f64 * ROW;
            let a = x(r.start.max(w0) as f64);
            let b = x(r.end.min(w1) as f64);
            let _ = write!(
                s.body,
                "<g><title>{}  HP {}  PS {}  {} strand  {}:{}-{}</title>",
                esc(&r.name),
                r.hp,
                r.ps.map_or(".".into(), |p| p.to_string()),
                if r.reverse { "-" } else { "+" },
                job.contig,
                r.start + 1,
                r.end
            );
            s.rect(a, ry, b - a, BAR, READ_BAR);
            for c in &r.calls {
                if meth::classify(c, ctx.threshold).is_none() {
                    continue; // below the confidence threshold: not counted, not drawn
                }
                let px = x(c.cpg as f64);
                s.rect(px - 0.8, ry, 1.9, BAR, &meth_color(c.p_m));
            }
            s.raw("</g>");
        }
        by += block_h(rs.len()) + 8.0;
    }

    // methylation panel
    let yv = |v: f64| meth_y + meth_h - v / 100.0 * meth_h;
    for v in [0.0, 25.0, 50.0, 75.0, 100.0] {
        s.line(
            ML,
            yv(v),
            ML + pw,
            yv(v),
            if v == 0.0 { AXIS } else { GRID },
            1.0,
        );
        s.text(
            ML - 8.0,
            yv(v) + 3.5,
            &format!("{v:.0}%"),
            10.0,
            MUTED,
            "end",
            false,
        );
    }
    s.text(
        16.0,
        meth_y + meth_h / 2.0 - 4.0,
        "5mC",
        11.0,
        INK_2,
        "start",
        true,
    );
    s.text(
        16.0,
        meth_y + meth_h / 2.0 + 10.0,
        "per CpG",
        10.0,
        MUTED,
        "start",
        false,
    );
    if let Some((lo, hi)) = l.expected {
        for v in if (hi - lo).abs() < 0.5 {
            vec![lo]
        } else {
            vec![lo, hi]
        } {
            s.line(cx0, yv(v), cx1, yv(v), MUTED, 1.0);
        }
    }
    let have_phased = core[Group::Hp1.index()].reads + core[Group::Hp2.index()].reads > 0;
    let series: Vec<Group> = if have_phased {
        vec![Group::All, Group::Hp1, Group::Hp2]
    } else {
        vec![Group::All]
    };
    for g in &series {
        let sites = &meth.sites[g.index()];
        if *g != Group::All || !have_phased {
            for (p, site) in meth.cpgs.iter().zip(sites) {
                if let Some(v) = site.pct() {
                    s.circle(
                        x(*p as f64),
                        yv(v),
                        1.7,
                        group_color(*g),
                        r#"fill-opacity="0.35""#,
                    );
                }
            }
        }
        let width = if *g == Group::All { 1.25 } else { 2.0 };
        s.raw(&smooth_path(
            &meth.cpgs,
            sites,
            &x,
            &yv,
            group_color(*g),
            width,
        ));
    }
    // legend for the panel (>= 2 series)
    let mut lx = ML + pw;
    if l.expected.is_some() {
        lx -= text_w("expected", 10.0) + 4.0;
        s.text(
            lx + text_w("expected", 10.0),
            meth_y - 5.0,
            "expected",
            10.0,
            INK_2,
            "end",
            false,
        );
        lx -= 16.0;
        s.line(lx, meth_y - 8.5, lx + 12.0, meth_y - 8.5, MUTED, 1.0);
        lx -= 12.0;
    }
    for g in series.iter().rev() {
        let t = g.label();
        lx -= text_w(t, 10.0) + 4.0;
        s.text(
            lx + text_w(t, 10.0),
            meth_y - 5.0,
            t,
            10.0,
            INK_2,
            "end",
            false,
        );
        lx -= 16.0;
        s.line(
            lx,
            meth_y - 8.5,
            lx + 12.0,
            meth_y - 8.5,
            group_color(*g),
            2.0,
        );
        lx -= 12.0;
    }

    // x axis
    let step = nice_step(span, 7.0).max(1.0);
    let mut t = ((w0 as f64) / step).ceil() * step;
    while t <= w1 as f64 {
        let px = x(t);
        s.line(px, axis_y, px, axis_y + 4.0, AXIS, 1.0);
        s.text(
            px,
            axis_y + 16.0,
            &thousands(t as u64),
            10.0,
            MUTED,
            "middle",
            false,
        );
        t += step;
    }
    s.text(
        ML + pw,
        axis_y + 32.0,
        &format!("{} ({})", job.contig, ctx.build),
        10.0,
        MUTED,
        "end",
        false,
    );
    // key for the read colours
    let ky = axis_y + 28.0;
    s.text(
        16.0,
        ky + 4.0,
        "read CpG call:",
        10.0,
        INK_2,
        "start",
        false,
    );
    let mut kx = 16.0 + text_w("read CpG call:", 10.0) + 8.0;
    for (p, lab) in [(0.0, "unmethylated"), (1.0, "methylated")] {
        s.rect(kx, ky - 4.0, 10.0, 10.0, &meth_color(p));
        s.text(kx + 14.0, ky + 4.5, lab, 10.0, INK_2, "start", false);
        kx += 14.0 + text_w(lab, 10.0) + 14.0;
    }
    s.text(
        kx,
        ky + 4.5,
        &format!(
            "(calls below p {:.2} not drawn; shaded = locus)",
            ctx.threshold
        ),
        10.0,
        MUTED,
        "start",
        false,
    );
    Ok(s.finish(W, ky + 16.0))
}

/// Smoothed 5mC line; broken where coverage runs out or CpGs are far apart.
fn smooth_path(
    cpgs: &[u64],
    sites: &[Site],
    x: &dyn Fn(f64) -> f64,
    yv: &dyn Fn(f64) -> f64,
    color: &str,
    width: f64,
) -> String {
    let sm = meth::smooth(sites, 4, 3.0);
    let mut d = String::new();
    let mut prev: Option<u64> = None;
    for (p, v) in cpgs.iter().zip(&sm) {
        match v {
            Some(v) => {
                let cmd = if prev.is_some_and(|q| p - q <= 600) {
                    'L'
                } else {
                    'M'
                };
                let _ = write!(d, "{cmd}{:.2},{:.2}", x(*p as f64), yv(*v));
                prev = Some(*p);
            }
            None => prev = None,
        }
    }
    if d.is_empty() {
        return String::new();
    }
    format!(
        r#"<path d="{d}" fill="none" stroke="{color}" stroke-width="{width}" stroke-linejoin="round" stroke-linecap="round"/>"#
    )
}

const GENE_ROWS: usize = 4;

/// Gene track: one transcript per gene, stacked so labels don't collide (at
/// most GENE_ROWS rows; genes that don't fit are left out and counted).
/// Returns the height used (0 without genes).
fn draw_genes(
    s: &mut Svg,
    tx: &[&crate::gtf::Transcript],
    y: f64,
    w0: u64,
    w1: u64,
    x: &dyn Fn(f64) -> f64,
) -> f64 {
    if tx.is_empty() {
        return 0.0;
    }
    let mut rows: Vec<f64> = Vec::new(); // right edge used per row
    let mut hidden = 0;
    for t in tx {
        let a = x(t.start.max(w0) as f64);
        let b = x(t.end.min(w1) as f64);
        let label = format!("{} {}", t.gene, if t.strand == '-' { "←" } else { "→" });
        let lw = text_w(&label, 10.0);
        let row = match rows.iter().position(|r| *r + 6.0 < a) {
            Some(r) => r,
            None if rows.len() < GENE_ROWS => {
                rows.push(f64::MIN);
                rows.len() - 1
            }
            None => {
                hidden += 1;
                continue;
            }
        };
        rows[row] = b.max(a + lw);
        let gy = y + 4.0 + row as f64 * 24.0;
        s.text(a.max(ML), gy + 9.0, &label, 10.0, INK_2, "start", false);
        s.line(a, gy + 17.0, b, gy + 17.0, INK_2, 1.0);
        for (e0, e1) in &t.exons {
            if *e1 <= w0 || *e0 >= w1 {
                continue;
            }
            let (ea, eb) = (x((*e0).max(w0) as f64), x((*e1).min(w1) as f64));
            s.rect(ea, gy + 13.0, (eb - ea).max(1.0), 8.0, INK_2);
        }
    }
    if hidden > 0 {
        s.text(
            ML - 8.0,
            y + 13.0,
            &format!("+{hidden} genes"),
            9.5,
            MUTED,
            "end",
            false,
        );
    }
    rows.len() as f64 * 24.0 + 6.0
}

// ---- overview -----------------------------------------------------------------
const OW: f64 = 940.0;
const OL: f64 = 300.0; // label column
const OR: f64 = 150.0; // counts column
const OROW: f64 = 20.0;

pub fn overview_svg(ctx: &Context, done: &[Done], sex: Option<Sex>) -> String {
    let pw = OW - OL - OR;
    let xv = |v: f64| OL + v / 100.0 * pw;
    let mut s = Svg::new();
    s.text(
        16.0,
        26.0,
        &format!(
            "{}  ·  5mC per haplotype at {} loci ({})",
            ctx.sample,
            done.len(),
            ctx.build
        ),
        15.0,
        INK,
        "start",
        true,
    );
    s.text(
        16.0,
        44.0,
        "Pooled 5mC over each locus's CpGs. HP1/HP2 are phase-block labels, not parental alleles.",
        11.0,
        INK_2,
        "start",
        false,
    );
    let note = match sex {
        Some(Sex::Female) => String::new(),
        Some(Sex::Male) => "X-inactivation loci not shown (male).".into(),
        None => "X-inactivation loci not shown (sex unknown).".into(),
    };
    // legend
    let ly = 64.0;
    let mut lx = 16.0;
    for (g, lab) in [
        (Group::Hp1, "HP1"),
        (Group::Hp2, "HP2"),
        (Group::All, "all reads (no phased reads)"),
    ] {
        s.circle(
            lx + 5.0,
            ly - 4.0,
            4.5,
            group_color(g),
            &format!(r#"stroke="{SURFACE}" stroke-width="2""#),
        );
        s.text(lx + 14.0, ly, lab, 11.0, INK, "start", false);
        lx += 14.0 + text_w(lab, 11.0) + 18.0;
    }
    s.circle(
        lx + 5.0,
        ly - 4.0,
        3.5,
        SURFACE,
        &format!(r#"stroke="{INK_2}" stroke-width="2""#),
    );
    s.text(
        lx + 14.0,
        ly,
        "fewer than 5 reads",
        11.0,
        INK,
        "start",
        false,
    );
    lx += 14.0 + text_w("fewer than 5 reads", 11.0) + 18.0;
    s.rect(lx + 4.0, ly - 9.0, 2.0, 10.0, INK_2);
    s.text(lx + 12.0, ly, "expected", 11.0, INK, "start", false);
    if !note.is_empty() {
        s.text(16.0, ly + 18.0, &note, 11.0, MUTED, "start", false);
    }

    let mut y = ly + if note.is_empty() { 26.0 } else { 42.0 };
    // axis labels on top
    let axis = |s: &mut Svg, y: f64| {
        for v in [0.0, 25.0, 50.0, 75.0, 100.0] {
            s.text(xv(v), y, &format!("{v:.0}%"), 10.0, MUTED, "middle", false);
        }
    };
    axis(&mut s, y);
    s.text(OW - 16.0, y, "reads HP1 / HP2", 10.0, MUTED, "end", false);
    y += 10.0;

    // sections in first-seen order (panels are sorted by position within)
    let mut groups: Vec<&str> = Vec::new();
    for d in done {
        if !groups.contains(&d.job.locus.group.as_str()) {
            groups.push(&d.job.locus.group);
        }
    }
    let order = |g: &str| match g {
        "Imprinting" => 0,
        "X-inactivation" => 1,
        "Controls" => 2,
        _ => 3,
    };
    groups.sort_by_key(|g| (order(g), g.to_string()));
    for g in groups {
        let rows: Vec<&Done> = done.iter().filter(|d| d.job.locus.group == g).collect();
        y += 16.0;
        s.text(
            16.0,
            y,
            &format!("{g} ({})", rows.len()),
            12.0,
            INK,
            "start",
            true,
        );
        y += 6.0;
        for d in rows {
            y += OROW;
            let cy = y - 6.0;
            let l = &d.job.locus;
            let c = &d.meth.core;
            s.raw("<g>");
            for v in [0.0, 25.0, 50.0, 75.0, 100.0] {
                s.line(xv(v), cy - 8.0, xv(v), cy + 8.0, GRID, 1.0);
            }
            let name = fit(&l.name.replace('_', " "), 11.0, 150.0);
            s.text(16.0, cy + 4.0, &name, 11.0, INK, "start", false);
            let extra = if real(&l.disease) {
                l.disease.replace('_', " ")
            } else if real(&l.gene) {
                l.gene.clone()
            } else {
                String::new()
            };
            if !extra.is_empty() {
                s.text(
                    172.0,
                    cy + 4.0,
                    &fit(&extra, 10.0, OL - 190.0),
                    10.0,
                    MUTED,
                    "start",
                    false,
                );
            }
            if let Some((lo, hi)) = l.expected {
                for v in [lo, hi] {
                    s.rect(xv(v) - 1.0, cy - 6.0, 2.0, 12.0, INK_2);
                }
            }
            let (h1, h2) = (&c[Group::Hp1.index()], &c[Group::Hp2.index()]);
            let mut tip = format!("{}  {}", l.name, l.coords());
            let pts: Vec<(Group, f64, usize)> = if h1.pct().is_some() || h2.pct().is_some() {
                [(Group::Hp1, h1), (Group::Hp2, h2)]
                    .iter()
                    .filter_map(|(g, sm)| sm.pct().map(|v| (*g, v, sm.reads)))
                    .collect()
            } else {
                let a = &c[Group::All.index()];
                a.pct()
                    .map(|v| vec![(Group::All, v, a.reads)])
                    .unwrap_or_default()
            };
            // HP1 a little above, HP2 a little below, so equal values stay visible
            let dy = |g: Group| match g {
                Group::Hp1 => -2.5,
                Group::Hp2 => 2.5,
                _ => 0.0,
            };
            if pts.len() == 2 {
                s.line(
                    xv(pts[0].1),
                    cy + dy(pts[0].0),
                    xv(pts[1].1),
                    cy + dy(pts[1].0),
                    AXIS,
                    2.0,
                );
            }
            for (g, v, n) in &pts {
                let _ = write!(tip, "  {} {v:.1}% ({n} reads)", g.label());
                if *n < 5 {
                    s.circle(
                        xv(*v),
                        cy + dy(*g),
                        3.5,
                        SURFACE,
                        &format!(r#"stroke="{}" stroke-width="2""#, group_color(*g)),
                    );
                } else {
                    s.circle(
                        xv(*v),
                        cy + dy(*g),
                        4.5,
                        group_color(*g),
                        &format!(r#"stroke="{SURFACE}" stroke-width="2""#),
                    );
                }
            }
            if pts.is_empty() {
                s.text(xv(50.0), cy + 4.0, "no reads", 10.0, MUTED, "middle", false);
            }
            let counts = if h1.reads + h2.reads > 0 {
                format!("{} / {}", h1.reads, h2.reads)
            } else {
                format!("{} unphased", c[Group::All.index()].reads)
            };
            s.text(OW - 16.0, cy + 4.0, &counts, 10.5, INK_2, "end", false);
            let _ = write!(s.body, "<title>{}</title></g>", esc(&tip));
        }
    }
    y += 22.0;
    axis(&mut s, y);
    y += 20.0;
    s.text(
        16.0,
        y,
        &format!(
            "hapmeth {}  ·  call threshold p {:.2}  ·  one figure per locus in plots/",
            env!("CARGO_PKG_VERSION"),
            ctx.threshold
        ),
        10.0,
        MUTED,
        "start",
        false,
    );
    s.finish(OW, y + 12.0)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ramp_ends() {
        assert_eq!(meth_color(0.0), "#a8a79f");
        assert_eq!(meth_color(1.0), "#0b0b0b");
        assert_eq!(meth_color(0.5), "#52514e");
    }

    #[test]
    fn formatting() {
        assert_eq!(thousands(2_698_551), "2,698,551");
        assert_eq!(thousands(12), "12");
        assert_eq!(nice_step(5900.0, 7.0), 1000.0);
        assert!(fit("a-very-long-locus-name-indeed", 11.0, 60.0).ends_with('…'));
    }
}
