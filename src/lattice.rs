//! Urea crystal lattice: tetragonal P-42_1m, a = b = 5.565 A, c = 4.684 A, Z = 2
//! (Swaminathan, Craven, Spackman & Stewart, Acta Cryst. B40, 398 (1984)).
//!
//! One kMC site per molecule. Two sublattices related by the -4 axis along c:
//!   alpha (s = 0): C at (0, 1/2, +zC), C=O points +c, molecular plane contains [110]
//!   beta  (s = 1): C at (1/2, 0, -zC), C=O points -c, molecular plane contains [1-10]
//! so the "dipole up/down" of the paper is the sublattice (paper: "the orientation of
//! the dipole is equivalent to defining the plane in which the urea molecule is situated").
//!
//! First neighbors (Fig. 2 of the paper; six hydrogen-bonded partners):
//!   A    chain neighbor along c on the molecule's NH2 side (it donates nothing; its O takes our N-H), as Fig. 2
//!        draws it. The SI rate rows are read the other way round by default (sim.rs, --ab co; README).
//!   B    chain neighbor along c on the molecule's C=O side
//!   1,1  the two in-plane partners along the molecule's OWN plane diagonal: our N-H donate to their O
//!   2,2  the two in-plane partners along the perpendicular diagonal: their N-H donate to our O
//! All four in-plane partners sit 0.348 c toward the molecule's C=O side (checked from the atomic
//! coordinates: every N...O distance is 2.99 A).
//!
//! ASSUMPTION (README, "1 vs 2"): which pair the paper calls 1 and which 2 is read from the projection
//! in Fig. 2 (the 2s sit next to the carbonyl, the 1s further out along the molecular plane).
//! `swap12` in the rate lookup tests the other reading.
//!
//! Second neighbors (paper's pseudo-second-neighbor split, SI):
//!   Z  group: the four non-bonded images of the in-plane partners, displaced 0.652 c the other way (4.98 A)
//!   XY group: the four same-sublattice molecules at +-a, +-b (5.565 A)
//! ASSUMPTION (README, "split"): this is our reading of "second nearest neighbors on the Z direction"
//! versus "on the X or Y direction".

pub const A_LEN: f64 = 5.565; // Angstrom
pub const C_LEN: f64 = 4.684; // Angstrom
pub const Z_C: f64 = 0.3260; // fractional z of the carbon atom

/// Molar concentration of crystalline urea: 2 molecules per a*a*c cell.
pub fn crystal_molarity() -> f64 {
    let v_mol_a3 = A_LEN * A_LEN * C_LEN / 2.0;
    1e27 / v_mol_a3 / 6.022_140_76e23
}

/// (target sublattice, di, dj, dk)
pub type Off = (usize, i64, i64, i64);

/// First neighbors: [type1, type1, type2, type2, A, B]
pub const NN: [[Off; 6]; 2] = [
    // alpha (i, j+1/2, k+zC); C=O points +c
    [(1, 0, 1, 0), (1, -1, 0, 0), (1, 0, 0, 0), (1, -1, 1, 0), (0, 0, 0, -1), (0, 0, 0, 1)],
    // beta (i+1/2, j, k+1-zC); C=O points -c
    [(0, 0, 0, 0), (0, 1, -1, 0), (0, 0, -1, 0), (0, 1, 0, 0), (1, 0, 0, 1), (1, 0, 0, -1)],
];

/// Second neighbors: [Z, Z, Z, Z, XY, XY, XY, XY]
pub const NNN: [[Off; 8]; 2] = [
    [
        (1, 0, 1, -1), (1, -1, 0, -1), (1, 0, 0, -1), (1, -1, 1, -1),
        (0, 1, 0, 0), (0, -1, 0, 0), (0, 0, 1, 0), (0, 0, -1, 0),
    ],
    [
        (0, 0, 0, 1), (0, 1, -1, 1), (0, 0, -1, 1), (0, 1, 0, 1),
        (1, 1, 0, 0), (1, -1, 0, 0), (1, 0, 1, 0), (1, 0, -1, 0),
    ],
];

/// Two geometries share one index layout (s, i, j, k):
///
/// * box: (i, j, k) are unit cells along a, b, c; the outer cell layer is inactive.
/// * slab: a [110] face, periodic along [1-10] and c. Index i is the position n along the layer
///   (period nx), j is the (110) layer number m (ny layers, the outer two inactive), k the cell along c
///   (period nz). Crystal coordinates are a = n, b = m - n. Both molecules of a cell with a + b = m sit on
///   the same plane x + y = (m + 1/2) a, so layer m is one complete (110) molecular layer, 3.935 A thick.
///   Wrapping n by nx is the lattice translation (nx, -nx, 0), so the face is a true periodic crystal face.
#[derive(Clone, Copy, Debug)]
pub struct Lattice {
    pub nx: usize,
    pub ny: usize,
    pub nz: usize,
    pub slab: bool,
}

/// (110) interlayer spacing in nm.
pub fn d110_nm() -> f64 {
    A_LEN / std::f64::consts::SQRT_2 / 10.0
}

impl Lattice {
    pub fn boxed(nx: usize, ny: usize, nz: usize) -> Lattice {
        Lattice { nx, ny, nz, slab: false }
    }

    /// `period` molecules-pairs along [1-10] (0.787 nm each), `layers` (110) layers, `nz` cells along c.
    pub fn slab(period: usize, layers: usize, nz: usize) -> Lattice {
        Lattice { nx: period, ny: layers, nz, slab: true }
    }

    pub fn n_sites(&self) -> usize {
        2 * self.nx * self.ny * self.nz
    }

    #[inline]
    pub fn idx(&self, s: usize, i: usize, j: usize, k: usize) -> usize {
        (((k * self.ny + j) * self.nx + i) << 1) | s
    }

    #[inline]
    pub fn coords(&self, idx: usize) -> (usize, usize, usize, usize) {
        let s = idx & 1;
        let c = idx >> 1;
        let i = c % self.nx;
        let j = (c / self.nx) % self.ny;
        let k = c / (self.nx * self.ny);
        (s, i, j, k)
    }

    /// Sites on the outer cell layer are never active, so every active site has all
    /// first and second neighbors inside the box (all offsets are within +-1 cell).
    #[inline]
    pub fn interior(&self, i: usize, j: usize, k: usize) -> bool {
        if self.slab {
            return j > 0 && j + 1 < self.ny; // every offset changes the layer by at most one
        }
        i > 0 && j > 0 && k > 0 && i + 1 < self.nx && j + 1 < self.ny && k + 1 < self.nz
    }

    /// Within two cells of the box wall: growth there ends the run.
    #[inline]
    pub fn near_wall(&self, i: usize, j: usize, k: usize) -> bool {
        if self.slab {
            return j < 2 || j + 2 >= self.ny;
        }
        i < 2 || j < 2 || k < 2 || i + 2 >= self.nx || j + 2 >= self.ny || k + 2 >= self.nz
    }

    /// Offset neighbor; caller guarantees (i, j, k) is interior.
    #[inline]
    pub fn shift(&self, i: usize, j: usize, k: usize, o: Off) -> usize {
        if self.slab {
            // crystal offset (da, db) moves the layer by da + db and the in-layer position by da
            let n = (i as i64 + o.1).rem_euclid(self.nx as i64) as usize;
            let m = (j as i64 + o.1 + o.2) as usize;
            let kk = (k as i64 + o.3).rem_euclid(self.nz as i64) as usize;
            return self.idx(o.0, n, m, kk);
        }
        let ii = (i as i64 + o.1) as usize;
        let jj = (j as i64 + o.2) as usize;
        let kk = (k as i64 + o.3) as usize;
        self.idx(o.0, ii, jj, kk)
    }

    /// Cartesian position of the carbon atom, in nanometers.
    pub fn position_nm(&self, idx: usize) -> (f64, f64, f64) {
        let (s, i, j, k) = self.coords(idx);
        let (i, j) = if self.slab { (i as f64, j as f64 - i as f64) } else { (i as f64, j as f64) };
        let (fx, fy, fz) = if s == 0 {
            (i, j + 0.5, k as f64 + Z_C)
        } else {
            (i + 0.5, j, k as f64 + 1.0 - Z_C)
        };
        (fx * A_LEN / 10.0, fy * A_LEN / 10.0, fz * C_LEN / 10.0)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Sites to test: the box center, and for the slab also sites on the wrap-around edges.
    fn cases() -> Vec<(Lattice, usize, usize, usize)> {
        let b = Lattice::boxed(8, 8, 8);
        let sl = Lattice::slab(8, 8, 8);
        vec![(b, 4, 4, 4), (sl, 4, 4, 4), (sl, 0, 4, 0), (sl, 7, 3, 7), (sl, 0, 2, 7)]
    }

    #[test]
    fn first_neighbors_are_mutual_with_swapped_types() {
        // If Y is X's type-1 partner (X donates to Y), X is Y's type-2 partner, and vice versa;
        // if Y is X's A, X is Y's B.
        let kind = |n: usize| match n { 0 | 1 => 1, 2 | 3 => 2, 4 => 10, _ => 20 };
        for (l, i, j, k) in cases() {
            for s in 0..2 {
                let x = l.idx(s, i, j, k);
                for (n, o) in NN[s].iter().enumerate() {
                    let y = l.shift(i, j, k, *o);
                    let (sy, iy, jy, ky) = l.coords(y);
                    let back: Vec<usize> = (0..6).filter(|&m| l.shift(iy, jy, ky, NN[sy][m]) == x).collect();
                    assert_eq!(back.len(), 1, "{l:?} ({i},{j},{k}) s={s} n={n}");
                    let want = match kind(n) { 1 => 2, 2 => 1, 10 => 20, _ => 10 };
                    assert_eq!(kind(back[0]), want, "{l:?} s={s} n={n}");
                }
            }
        }
    }

    #[test]
    fn second_neighbors_are_mutual() {
        for (l, i, j, k) in cases() {
            for s in 0..2 {
                let x = l.idx(s, i, j, k);
                for (n, o) in NNN[s].iter().enumerate() {
                    let y = l.shift(i, j, k, *o);
                    let (sy, iy, jy, ky) = l.coords(y);
                    let back: Vec<usize> = (0..8).filter(|&m| l.shift(iy, jy, ky, NNN[sy][m]) == x).collect();
                    assert_eq!(back.len(), 1, "{l:?} s={s} n={n}");
                    assert_eq!(back[0] < 4, n < 4, "Z stays Z, XY stays XY");
                }
            }
        }
    }

    #[test]
    fn distances_match_the_crystal_structure() {
        // center sites only: across a periodic edge the stored positions are a period apart
        for l in [Lattice::boxed(8, 8, 8), Lattice::slab(8, 8, 8)] {
            let d = |a: usize, b: usize| {
                let (x1, y1, z1) = l.position_nm(a);
                let (x2, y2, z2) = l.position_nm(b);
                ((x1 - x2).powi(2) + (y1 - y2).powi(2) + (z1 - z2).powi(2)).sqrt() * 10.0
            };
            for s in 0..2 {
                let x = l.idx(s, 4, 4, 4);
                for (n, o) in NN[s].iter().enumerate() {
                    let r = d(x, l.shift(4, 4, 4, *o));
                    let want = if n < 4 { 4.259 } else { 4.684 };
                    assert!((r - want).abs() < 0.01, "NN s={s} n={n} r={r}");
                }
                for (n, o) in NNN[s].iter().enumerate() {
                    let r = d(x, l.shift(4, 4, 4, *o));
                    let want = if n < 4 { 4.98 } else { 5.565 };
                    assert!((r - want).abs() < 0.01, "NNN s={s} n={n} r={r}");
                }
                // in-plane partners sit toward the molecule's own C=O side
                let zc = l.position_nm(x).2;
                for o in &NN[s][..4] {
                    let dz = l.position_nm(l.shift(4, 4, 4, *o)).2 - zc;
                    assert_eq!(dz > 0.0, s == 0, "alpha C=O points +c, beta -c");
                }
            }
        }
    }

    #[test]
    fn slab_layers_are_110_planes() {
        // every molecule of layer m sits at x + y = (m + 1/2) a, i.e. (m + 1/2) * d110 along [110]
        let l = Lattice::slab(8, 8, 8);
        for idx in 0..l.n_sites() {
            let (_, _, m, _) = l.coords(idx);
            let (x, y, _) = l.position_nm(idx);
            let u = (x + y) / std::f64::consts::SQRT_2;
            assert!((u - (m as f64 + 0.5) * d110_nm()).abs() < 1e-9);
        }
    }

    #[test]
    fn molarity() {
        assert!((crystal_molarity() - 22.89).abs() < 0.01);
    }
}
