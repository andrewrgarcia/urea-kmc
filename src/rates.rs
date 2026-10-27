//! Rate lookup from the decoded SI tables, plus the modeling switches the paper leaves open.

use crate::rates_table::*;

#[derive(Clone, Copy, PartialEq, Debug)]
pub enum Solvent {
    Water,
    Methanol,
}

/// How to use the three-valued (pseudo-second-neighbor) cells.
#[derive(Clone, Copy, PartialEq, Debug)]
pub enum Split {
    /// Classify by second neighbors (lattice.rs, "Second neighbors") and use i / ii / iii.
    Nnn,
    /// Nearest-neighbor model: arithmetic mean of i, ii, iii.
    Mean,
}

/// Cells that were never observed in MD (0 in both the dissolution and the crystallization table).
#[derive(Clone, Copy, PartialEq, Debug)]
pub enum Unobserved {
    /// Geometric mean of the observed cells with the same neighbor count. The six-neighbor
    /// (bulk) cell keeps zero dissolution; a buried vacancy fills like a five-neighbor site.
    Fill,
    /// Keep the published zeros: such a molecule can never leave, such a hole never fills.
    Published,
}

pub fn col(n1: usize, n2: usize) -> usize {
    match (n1, n2) {
        (0, 0) => 0,
        (1, 0) => 1,
        (2, 0) => 2,
        (0, 1) => 3,
        (1, 1) => 4,
        (2, 1) => 5,
        (0, 2) => 6,
        (1, 2) => 7,
        (2, 2) => 8,
        _ => unreachable!("at most two partners of each in-plane type"),
    }
}

pub fn row(a: bool, b: bool) -> usize {
    a as usize + 2 * b as usize
}

pub const COL_NAMES: [&str; 9] = ["-", "1", "11", "2", "12", "112", "22", "122", "1122"];
pub const ROW_NAMES: [&str; 4] = ["-", "A", "B", "AB"];

pub fn neighbors(r: usize, c: usize) -> usize {
    let inplane = [0, 1, 2, 1, 2, 3, 2, 3, 4][c];
    let caxis = [0, 1, 1, 2][r];
    inplane + caxis
}

pub struct Rates {
    pub kd: Table,
    pub kc: Table,
    pub split_on: [[bool; 9]; 4],
    pub filled: Vec<String>,
}

fn geomean(v: &[f64]) -> f64 {
    (v.iter().map(|x| x.ln()).sum::<f64>() / v.len() as f64).exp()
}

impl Rates {
    pub fn new(solvent: Solvent, split: Split, unobserved: Unobserved) -> Rates {
        let (mut kd, mut kc) = match solvent {
            Solvent::Water => (WATER_DISSOLUTION, WATER_CRYSTALLIZATION),
            Solvent::Methanol => (METHANOL_DISSOLUTION, METHANOL_CRYSTALLIZATION),
        };
        let mut split_on = SPLIT;
        if split == Split::Mean {
            for r in 0..4 {
                for c in 0..9 {
                    if split_on[r][c] {
                        let md = kd[r][c].iter().sum::<f64>() / 3.0;
                        let mc = kc[r][c].iter().sum::<f64>() / 3.0;
                        kd[r][c] = [md; 3];
                        kc[r][c] = [mc; 3];
                        split_on[r][c] = false;
                    }
                }
            }
        }
        let mut filled = Vec::new();
        if unobserved == Unobserved::Fill {
            let observed = |n: usize, t: &Table| -> Vec<f64> {
                let mut v = Vec::new();
                for r in 0..4 {
                    for c in 0..9 {
                        if (r, c) != (0, 0) && neighbors(r, c) == n {
                            v.extend(t[r][c].iter().copied().filter(|x| *x > 0.0));
                        }
                    }
                }
                v
            };
            let (kd0, kc0) = (kd, kc);
            for r in 0..4 {
                for c in 0..9 {
                    let n = neighbors(r, c);
                    if kd0[r][c][0] == 0.0 && kc0[r][c][0] == 0.0 && n > 0 {
                        if n < 6 {
                            kd[r][c] = [geomean(&observed(n, &kd0)); 3];
                            kc[r][c] = [geomean(&observed(n, &kc0)); 3];
                        } else {
                            kc[r][c] = [geomean(&observed(5, &kc0)); 3];
                        }
                        filled.push(format!(
                            "{}/{} (n={n}): kd={:.3e} kc={:.3e}",
                            COL_NAMES[c], ROW_NAMES[r], kd[r][c][0], kc[r][c][0]
                        ));
                    }
                }
            }
        }
        Rates { kd, kc, split_on, filled }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn pinned_values_from_the_supplement() {
        let w = Rates::new(Solvent::Water, Split::Nnn, Unobserved::Published);
        // "12/B" water dissolution: 0.4868E+09 0.1353E+10 0.3241E+10
        assert_eq!(w.kd[row(false, true)][col(1, 1)], [4.868e8, 1.353e9, 3.241e9]);
        // "1122/B" water crystallization: 0.1953E+12
        assert_eq!(w.kc[row(false, true)][col(2, 2)][0], 1.953e11);
        // "112/AB" water dissolution: 0.1905E+08
        assert_eq!(w.kd[row(true, true)][col(2, 1)][0], 1.905e7);
        let m = Rates::new(Solvent::Methanol, Split::Nnn, Unobserved::Published);
        // "1/AB" methanol crystallization: 0.3358E+11
        assert_eq!(m.kc[row(true, true)][col(1, 0)][0], 3.358e10);
        // "122/A" methanol dissolution: 0.8228E+10 0.3055E+10 0.3163E+09
        assert_eq!(m.kd[row(true, false)][col(1, 2)], [8.228e9, 3.055e9, 3.163e8]);
    }

    #[test]
    fn site_type_count_matches_the_paper() {
        // 34 surface types: 10 kinks (3), 8 steps (4), 4 terraces (5), 12 with one or two neighbors.
        let mut by_n = [0usize; 7];
        for r in 0..4 {
            for c in 0..9 {
                by_n[neighbors(r, c)] += 1;
            }
        }
        assert_eq!(by_n, [1, 4, 8, 10, 8, 4, 1]);
        assert_eq!(by_n[1] + by_n[2], 12);
    }

    #[test]
    fn bulk_never_dissolves() {
        for s in [Solvent::Water, Solvent::Methanol] {
            let r = Rates::new(s, Split::Nnn, Unobserved::Fill);
            assert_eq!(r.kd[3][8][0], 0.0);
            assert!(r.kc[3][8][0] > 0.0);
        }
    }
}
