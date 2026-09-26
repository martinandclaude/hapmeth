# hapmeth — per-haplotype methylation locus plots (long-read WGS)

`hapmeth` (formerly `banana`) plots 5mC per haplotype at a curated set of loci
from a haplotagged Oxford Nanopore modBAM:

- **Imprinting:** germline imprinted DMRs.
- **X-inactivation:** the MECP2, XIST and AR promoters.
- **QC controls:** methylated and unmethylated sanity checks.

It writes one PNG per locus, one overview PNG with every locus (sized for a QC
report), and a TSV with the numbers.

It is a single Rust binary. It reads the MM/ML tags itself (no modkit), draws
the figures itself (no methylartist), and has the panels for **hg38 and hs1**
built in.

```bash
hapmeth -b sample.haplotagged.bam -r GRCh38.fa        # or hs1.fa: build is detected
```

| overview (`<sample>.overview.png`) | one locus (`plots/Imprinting/KvDMR1.png`) |
|---|---|
| ![overview](docs/example_overview.png) | ![KvDMR1](docs/example_KvDMR1.png) |

<sub>HG001, ONT public GIAB data (`giab_2025.01`, dorado v5 SUP 5mC_5hmC).</sub>

## Install

```bash
cargo build --release          # -> target/release/hapmeth
# or: cargo install --path .
```

Rust ≥ 1.89. The binary needs only libc at run time: no conda environment,
htslib or Python. The panels and the font are compiled in.

## Inputs

| input | required | notes |
|---|---|---|
| haplotagged modBAM / CRAM | yes | `HP` tags + `MM`/`ML` tags; `.bai`/`.csi`/`.crai` alongside |
| reference FASTA | yes | the one the reads are aligned to; `.fai` (+ `.gzi` if bgzipped) |
| GTF | no | `--gtf`, plain or gzipped → gene track |

hg38 vs hs1 is taken from chr1's length in the BAM header (`--build` to force).
Contig names with or without `chr` both work. A BAM without `HP` tags still
plots, as a single "untagged" group. CRAM works too, but is slower: decoding
loads whole reference chromosomes.

## Output

```
<out>/<sample>.overview.png        every locus: HP1 vs HP2 5mC, expected values marked
<out>/<sample>.loci.tsv            one row per locus (the numbers behind the figures)
<out>/plots/<panel>/<locus>.png    per-locus figure
<out>/<sample>.cpgs.tsv            with --sites: per-CpG counts per haplotype
```

`--out-dir` defaults to `./<sample>`; `--format svg|both` adds SVGs. Same-named
loci in one folder get their coordinates appended.

**Per-locus figure.** From the top:
- an optional gene track;
- the CpGs in the window;
- the reads, grouped HP1 / HP2 / untagged, with each CpG call inked from light
  (unmethylated) to black (methylated);
- 5mC % per CpG for each haplotype (dots) with a smoothed line, and the expected
  values.

The locus itself is shaded. The header gives each haplotype's pooled 5mC over
the locus.

**`loci.tsv` columns:** `sample, name, group, gene, disease, chrom, start, end,
core_cpgs`, then for `hp1`/`hp2`: `_reads, _cpgs, _pct`, then `untagged_reads,
untagged_pct, all_reads, all_pct, hp_diff` (|HP1−HP2|), `expected_lo,
expected_hi, phase_sets, plot, status` (`ok` / `no_reads` / `plot_failed: …`),
`note` (why reads were dropped, e.g. `11 reads below --min-mapq 10`).
- 5mC % is pooled over the locus's own CpGs (not the padding).
- `_reads` counts reads with at least one call there.
- `phase_sets` counts the distinct `PS` values among the tagged reads. HP1 and
  HP2 are phase-block labels, not parental alleles; there is no trio.

## How the methylation is read

hapmeth parses `MM`/`ML` directly (`src/modbam.rs`), following SAMtags §1.7:
- **Sparse tags.** dorado v5 writes sparse tags (`C+m.`): cytosines that are not
  listed are *unmethylated*, and hapmeth counts them that way. `?` means
  unlisted cytosines are unknown.
- **Strand and codes.** Reverse-strand reads, multi-code entries (`C+hm`), ChEBI
  codes, and the ML values of entries it ignores (`A+a`) are all handled.
- **Stale tags.** Records whose `MN` disagrees with `SEQ`, or that were
  hard-clipped without `MN`, are skipped.

This is why banana needed modkit: methylartist's own parse counts only the
listed calls and overstates 5mC.

Counting matches `modkit pileup --cpg --combine-strands`:
- Each call is the most likely of C / 5mC / 5hmC / other, and calls below the
  confidence threshold are dropped.
- 5mC % = 5mC / (C + 5mC + 5hmC + other).
- The canonical probability uses modkit's convention, `(255 − ΣML + 0.5)/256`.
  This matters for calls sitting on the threshold.
- The threshold follows modkit's default rule: the 10th percentile of call
  confidence (`--filter-percentile`). hapmeth estimates it from the reads at the
  panel loci, so it can differ slightly from what modkit estimates genome-wide.
  Pass modkit's logged value with `--filter-threshold` for exact parity.
- Primary alignments only, and MAPQ ≥ 10 by default (`--min-mapq 1` for
  segmental-duplication loci such as CTAG1B).

**Validation.** The test data was ONT's public HG001 and HG002 runs
(`giab_2025.01`, dorado v5 SUP 5mC_5hmC, sparse `.` tags), sliced to all 151
panel regions ± 1.5 kb. hapmeth's `--sites` counts were compared with
`modkit pileup --cpg --combine-strands --phased` (modkit 0.6.4), both at
threshold 0.7246. That is the value modkit chose itself for these reads.
- **Result:** 106,807 of 106,824 CpG × {HP1, HP2, combined} rows were identical
  (about 2.6 M calls), in 5mC, 5hmC, canonical and filtered counts alike.
- **The other 17 rows,** plus one row only hapmeth reports: they are 3 CpGs
  whose G falls one base outside the region BED. There modkit's
  `--include-bed` drops the reverse-strand calls; it is not a parsing
  difference.
- **Synthetic modBAM:** `tests/end_to_end.rs` checks exact numbers on one that
  covers each case above.

## The panels

Built in, for **hg38 and hs1**:

| panel | loci | default |
|---|---|---|
| `imprinting` | 79 germline imprinted DMRs (atlas: which parental allele is methylated, expected high/low 5mC) | the 15 disease-anchored ones; `--all-imprinted` for all |
| `xci` | MECP2, XIST, AR promoter CpG islands | females only (sex from chrX/chrY depth in the BAM index; `--sex`) |
| `qc` | GAPDH, B2M, RPL13A (unmethylated on both haplotypes); CTAG1B (methylated on both) | always |

`--panels imprinting,qc` picks a subset, and `--panels none` plots only `--bed`
regions.

**Any other regions via `--bed`** (repeatable). A plain tab-separated BED3/BED4
works. A `#chrom start end name …` header naming extra columns (`gene disease
origin expected_lo expected_hi group`) adds labels, expected values and the
output subfolder; without `group`, the subfolder is the file name. For example,
`focal_dx.hg38.bed`:

```
#chrom	start	end	name	score	strand	gene	disease
chr3	36992737	36993865	MLH1	0	.	MLH1	Lynch constitutional epimutation (MLH1)
chr9	27572968	27573989	C9orf72	0	.	C9orf72	ALS/FTD (C9orf72 G4C2 5'CpGI hyper)
chr19	45767441	45771706	DMPK	0	.	DMPK	DM1 (DMPK CTG 3'CpGI hyper; congenital)
chrX	147911573	147912682	FMR1	0	.	FMR1	Fragile-X (FMR1 CGG, promoter hyper)
```

`--bed focal_dx.hg38.bed` adds these loci under `plots/focal_dx/`. The built-in
BEDs in `bed/` use the same format; `hapmeth --export-panels DIR` writes them
out.

**hs1 coordinates** were lifted from hg38 with the T2T consortium's GRCh38 →
CHM13v2.0 chain
([`grch38-chm13v2.chain`](https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/chain/v1_nflo/grch38-chm13v2.chain),
the same alignment UCSC ships as `hg38ToHs1.over.chain.gz`):
- All 151 regions lift to the same chromosome with ≥ 96 % of their bases
  aligned.
- Nearly all have ≥ 98 % sequence identity. The larger differences are real
  repeat-allele differences between the assemblies (FMR1, C9orf72, DMPK, XYLT1,
  and one imprinted DMR with a 2.2 kb insertion).
- The promoter CpG islands land exactly on hs1's own CpG-island track, next to
  the same genes (CTAG1B, not its CTAG1A paralog).
- On real reads realigned to hs1, the per-haplotype values match hg38 (KvDMR1
  1 % / 86 % vs 1 % / 87 %).

**Panel fixes in this version:**
- **MECP2** pointed at a CpG island 39 kb upstream of the gene, at the 5′ end of
  one long Ensembl transcript. It now uses the promoter CGI that contains the
  MANE TSS (chrX:154,097,110-154,098,015).
- **Default imprinting set.** H19 (IC1), IGF2 DMR0 and GNAS-XL now carry their
  disease labels (a join by gene name had missed them), so they are in the
  default set.

## Options (`--help` for all)

| option | purpose |
|---|---|
| `-o, --out-dir`, `-s, --sample` | output directory (default `./<sample>`), label (default: BAM file name) |
| `--panels`, `--all-imprinted`, `--bed`, `--locus` | choose loci |
| `--build auto\|hg38\|hs1`, `--sex auto\|female\|male` | override detection |
| `--gtf` | gene track |
| `--pad N` | flank plotted each side (default 1500 bp) |
| `--min-mapq`, `--include-supplementary` | read filters |
| `--filter-threshold`, `--filter-percentile` | call-confidence filter (modkit semantics) |
| `--max-reads N` | reads drawn per haplotype (default 60; all are counted) |
| `--sites`, `--format` | per-CpG table; png / svg / both |
| `-j, --threads` | default: all cores |

## Using it from a pipeline (e.g. pika)

hapmeth is a standalone command with fixed output names, so a pipeline only
needs to run it and pick up the files:
- `<out>/<sample>.overview.png` for a report;
- `<out>/<sample>.loci.tsv` for a table.

It exits non-zero on any error, including when no plot could be written. Pass
`--build` and `--sex` when the pipeline already knows them. A Snakemake rule
could look like this:

```python
rule hapmeth:
    input:  bam = "{sample}.haplotagged.bam", bai = "{sample}.haplotagged.bam.bai", ref = REF
    output: png = "qc/hapmeth/{sample}/{sample}.overview.png",
            tsv = "qc/hapmeth/{sample}/{sample}.loci.tsv"
    params: build = "hs1",                         # or hg38
            sex   = lambda wc: SEX.get(wc.sample, "auto")
    threads: 4
    log: "logs/hapmeth/{sample}.log"
    shell: "hapmeth -b {input.bam} -r {input.ref} -o qc/hapmeth/{wildcards.sample} "
           "-s {wildcards.sample} --build {params.build} --sex {params.sex} -j {threads} 2> {log}"
```

A whole WGS sample takes seconds with the default panels, because only the
panel windows are read.

## Adapting the panels

- **Add loci:** put them in a BED and pass it with `--bed`; no rebuild needed.
- **Edit a built-in panel:** change its BED in `bed/`, for both builds, then
  rebuild the binary; the built-ins are compiled in.
- **Lift a BED from hg38 to hs1:** use UCSC `liftOver` with the T2T chain linked
  above. Check the lifted lengths at repeats (FMR1, C9orf72, DMPK, XYLT1): the
  repeat alleles differ between the assemblies.

## Not carried over from banana

- `--control BAM`
- methylartist's variant ticks (`--vcf`, `--splitvar`)
- `--no-highlight`

Say if any of these are needed.

## Layout

```
src/            main.rs (CLI), modbam.rs (MM/ML), reads.rs (BAM/CRAM), meth.rs (counting),
                panel.rs (panels), plot.rs (figures), render.rs (SVG→PNG), gtf.rs
bed/            built-in panels (imprinting, xci, qc) as .hg38.bed / .hs1.bed, compiled in
assets/fonts/   DejaVu Sans (subset), embedded in the binary
tests/          end_to_end.rs (synthetic modBAM)
docs/           example figures
LICENSE         MIT
```

## Acknowledgements

The per-locus figure is modelled on the `locus` plots of
[methylartist](https://github.com/adamewing/methylartist), which banana used to
draw its figures:

> Cheetham SW, Kindlova M, Ewing AD. Methylartist: tools for visualizing
> modified bases from nanopore sequence data. *Bioinformatics*
> 38(11):3109–3112 (2022). https://doi.org/10.1093/bioinformatics/btac292

## Licence

hapmeth is MIT-licensed (`LICENSE`).

**Third-party code:**
- Rust crates: noodles (MIT); resvg, usvg and tiny-skia (Apache-2.0/MIT;
  BSD-3); clap, rayon, anyhow and flate2 (MIT/Apache-2.0).
- DejaVu Sans: Bitstream Vera licence and public domain; see
  `assets/fonts/LICENSE_DEJAVU.txt`.
- No modkit code (ONT Public License, research use only).
- No methylartist code (MIT); only the figure's layout idea is borrowed from
  it.
