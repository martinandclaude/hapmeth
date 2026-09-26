//! Calls -> per-CpG, per-haplotype 5mC and locus summaries.
//!
//! Counting follows modkit pileup (--cpg --combine-strands): a call is the most
//! likely of {C, 5mC, 5hmC, other}; calls whose probability is below the
//! confidence threshold are dropped (modkit's Nfail); 5mC % = 5mC / (C + 5mC +
//! 5hmC + other) over the remaining calls. The threshold is by default modkit's
//! rule too: the 10th percentile of call confidence (--filter-percentile).

use std::collections::BTreeSet;

use crate::reads::{Call, Read};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum State {
    Canonical,
    Methyl,
    Hydroxy,
    Other,
}

/// Most likely state, or None if its probability is below `threshold`.
pub fn classify(c: &Call, threshold: f32) -> Option<State> {
    let mut best = (c.p_c, State::Canonical);
    for (p, s) in [
        (c.p_m, State::Methyl),
        (c.p_h, State::Hydroxy),
        (c.p_other, State::Other),
    ] {
        if p > best.0 {
            best = (p, s);
        }
    }
    (best.0 >= threshold).then_some(best.1)
}

/// `percentile` (0..1) of the call confidences (modkit's --filter-percentile).
pub fn estimate_threshold(mut conf: Vec<f32>, percentile: f64) -> Option<f32> {
    if conf.is_empty() {
        return None;
    }
    conf.sort_unstable_by(f32::total_cmp);
    let i = ((conf.len() as f64 * percentile).floor() as usize).min(conf.len() - 1);
    Some(conf[i])
}

/// The groups summarised per locus.
pub const GROUPS: [Group; 4] = [Group::Hp1, Group::Hp2, Group::Unphased, Group::All];

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Group {
    Hp1,
    Hp2,
    Unphased,
    All,
}

impl Group {
    pub fn index(self) -> usize {
        self as usize
    }
    pub fn key(self) -> &'static str {
        match self {
            Group::Hp1 => "hp1",
            Group::Hp2 => "hp2",
            Group::Unphased => "unphased",
            Group::All => "all",
        }
    }
    pub fn label(self) -> &'static str {
        match self {
            Group::Hp1 => "HP1",
            Group::Hp2 => "HP2",
            Group::Unphased => "untagged",
            Group::All => "all reads",
        }
    }
    fn of(read: &Read) -> Group {
        match read.hp {
            1 => Group::Hp1,
            2 => Group::Hp2,
            _ => Group::Unphased,
        }
    }
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Site {
    pub m: u32,
    pub h: u32,
    pub c: u32,
    pub other: u32,
    pub fail: u32,
}

impl Site {
    pub fn valid(&self) -> u32 {
        self.m + self.h + self.c + self.other
    }
    pub fn pct(&self) -> Option<f64> {
        (self.valid() > 0).then(|| 100.0 * f64::from(self.m) / f64::from(self.valid()))
    }
    fn add(&mut self, s: Option<State>) {
        match s {
            Some(State::Methyl) => self.m += 1,
            Some(State::Hydroxy) => self.h += 1,
            Some(State::Canonical) => self.c += 1,
            Some(State::Other) => self.other += 1,
            None => self.fail += 1,
        }
    }
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct Summary {
    /// Reads with at least one passing call in the core region.
    pub reads: usize,
    /// Core CpGs with at least one passing call.
    pub cpgs_covered: usize,
    pub valid: u64,
    pub m: u64,
    pub h: u64,
}

impl Summary {
    /// Pooled 5mC % over the core CpGs.
    pub fn pct(&self) -> Option<f64> {
        (self.valid > 0).then(|| 100.0 * self.m as f64 / self.valid as f64)
    }
}

#[derive(Clone, Debug)]
pub struct LocusMeth {
    /// Reference CpG positions in the plotted window.
    pub cpgs: Vec<u64>,
    /// Per group, aligned with `cpgs`.
    pub sites: [Vec<Site>; 4],
    /// Per group, over the core (un-padded) region.
    pub core: [Summary; 4],
    /// Number of core CpGs in the reference.
    pub core_cpgs: usize,
    /// Distinct PS values among HP-tagged reads covering the core.
    pub phase_sets: usize,
}

pub fn aggregate(reads: &[Read], cpgs: Vec<u64>, core: (u64, u64), threshold: f32) -> LocusMeth {
    let n = cpgs.len();
    let mut sites: [Vec<Site>; 4] = std::array::from_fn(|_| vec![Site::default(); n]);
    let mut summ: [Summary; 4] = Default::default();
    let mut core_cov: [Vec<bool>; 4] = std::array::from_fn(|_| vec![false; n]);
    let mut ps = BTreeSet::new();
    let in_core = |p: u64| p >= core.0 && p < core.1;

    for read in reads {
        let g = Group::of(read).index();
        let mut covers_core = false;
        for c in &read.calls {
            let Ok(k) = cpgs.binary_search(&c.cpg) else {
                continue;
            };
            let s = classify(c, threshold);
            sites[g][k].add(s);
            sites[Group::All.index()][k].add(s);
            if let (Some(st), true) = (s, in_core(c.cpg)) {
                covers_core = true;
                for gi in [g, Group::All.index()] {
                    summ[gi].valid += 1;
                    summ[gi].m += u64::from(st == State::Methyl);
                    summ[gi].h += u64::from(st == State::Hydroxy);
                    core_cov[gi][k] = true;
                }
            }
        }
        if covers_core {
            summ[g].reads += 1;
            summ[Group::All.index()].reads += 1;
            if read.hp != 0 {
                if let Some(p) = read.ps {
                    ps.insert(p);
                }
            }
        }
    }
    for g in 0..4 {
        summ[g].cpgs_covered = core_cov[g].iter().filter(|&&b| b).count();
    }
    let core_cpgs = cpgs.iter().filter(|&&p| in_core(p)).count();
    LocusMeth {
        cpgs,
        sites,
        core: summ,
        core_cpgs,
        phase_sets: ps.len(),
    }
}

/// Coverage-weighted Hann smoothing over +/- `k` neighbouring CpGs; None where
/// the window holds fewer than `min_calls` weighted calls.
pub fn smooth(sites: &[Site], k: usize, min_calls: f64) -> Vec<Option<f64>> {
    let w: Vec<f64> = (0..=k)
        .map(|d| 0.5 * (1.0 + (std::f64::consts::PI * d as f64 / (k as f64 + 1.0)).cos()))
        .collect();
    (0..sites.len())
        .map(|i| {
            let (mut num, mut den) = (0.0, 0.0);
            for j in i.saturating_sub(k)..=(i + k).min(sites.len() - 1) {
                let wj = w[i.abs_diff(j)];
                num += wj * f64::from(sites[j].m);
                den += wj * f64::from(sites[j].valid());
            }
            (den >= min_calls).then(|| 100.0 * num / den)
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn call(cpg: u64, p_m: f32, p_h: f32) -> Call {
        Call {
            cpg,
            p_m,
            p_h,
            p_other: 0.0,
            p_c: 1.0 - p_m - p_h,
            explicit: true,
        }
    }

    #[test]
    fn classification_and_filter() {
        assert_eq!(classify(&call(0, 0.9, 0.05), 0.7), Some(State::Methyl));
        assert_eq!(classify(&call(0, 0.05, 0.9), 0.7), Some(State::Hydroxy));
        assert_eq!(classify(&call(0, 0.1, 0.1), 0.7), Some(State::Canonical));
        assert_eq!(classify(&call(0, 0.5, 0.2), 0.7), None); // 0.5 < 0.7
        assert_eq!(classify(&call(0, 0.0, 0.0), 0.99), Some(State::Canonical)); // implied call
    }

    #[test]
    fn threshold_percentile() {
        let v: Vec<f32> = (0..100).map(|i| i as f32 / 100.0).collect();
        assert_eq!(estimate_threshold(v, 0.1), Some(0.1));
        assert_eq!(estimate_threshold(vec![], 0.1), None);
    }

    #[test]
    fn aggregate_groups_and_core() {
        let read = |hp: u8, ps: i64, pm: f32| Read {
            name: String::new(),
            hp,
            ps: Some(ps),
            reverse: false,
            start: 0,
            end: 100,
            calls: vec![call(10, pm, 0.0), call(50, pm, 0.0), call(90, pm, 0.0)],
            conf: vec![],
        };
        let reads = vec![
            read(1, 7, 0.95),
            read(1, 7, 0.95),
            read(2, 7, 0.02),
            read(0, 7, 0.6),
        ];
        let m = aggregate(&reads, vec![10, 50, 90], (40, 95), 0.7);
        let hp1 = &m.core[Group::Hp1.index()];
        assert_eq!(
            (hp1.reads, hp1.cpgs_covered, hp1.valid, hp1.m),
            (2, 2, 4, 4)
        );
        assert_eq!(hp1.pct(), Some(100.0));
        assert_eq!(m.core[Group::Hp2.index()].pct(), Some(0.0));
        // untagged read's calls are all below threshold -> fail, no coverage
        assert_eq!(m.core[Group::Unphased.index()].reads, 0);
        assert_eq!(m.sites[Group::Unphased.index()][0].fail, 1);
        assert_eq!(m.core[Group::All.index()].pct(), Some(100.0 * 4.0 / 6.0));
        assert_eq!(m.core_cpgs, 2);
        assert_eq!(m.phase_sets, 1);
    }

    #[test]
    fn smoothing_weights_by_coverage() {
        let s = vec![
            Site {
                m: 10,
                c: 0,
                ..Default::default()
            },
            Site {
                m: 0,
                c: 1,
                ..Default::default()
            },
            Site::default(),
        ];
        let v = smooth(&s, 1, 1.0);
        // centre: (0.5*10 + 1*0 + 0.5*0) / (0.5*10 + 1*1 + 0) = 5/6
        assert!((v[1].unwrap() - 100.0 * 5.0 / 6.0).abs() < 1e-9);
        assert_eq!(smooth(&[Site::default()], 2, 1.0), vec![None]);
    }
}
