//! Exact continuous-time kMC (Bortz-Kalos-Lebowitz n-fold way) on the urea lattice.
//!
//! Every site is either crystal (occupied) or solution. A crystal site dissolves at kd(config);
//! a solution site with at least one crystal neighbor crystallizes at kc(config) * C, where C is
//! the urea concentration as a fraction of lattice sites (ASSUMPTION, README "units"), held constant
//! (an infinite, well-mixed reservoir; no diffusion). Unlike the paper's 50-ps synchronous sweep this
//! is serial-exact, so there is no "probability above 1" to clip.
//!
//! Sites with equal rate share a class: 2 kinds x 4 rows x 9 columns x 3 splits = 216 classes.

use crate::lattice::{Lattice, NN, NNN};
use crate::rates::{col, row, Rates};
use crate::rng::Rng;

pub const NONE: u8 = 255;
const N_CLASSES: usize = 216;
const HALF: usize = 108; // classes 0..108 dissolve, 108..216 crystallize

pub struct Sim {
    pub lat: Lattice,
    pub occ: Vec<u8>,
    cls: Vec<u8>,
    pos: Vec<u32>,
    members: Vec<Vec<u32>>,
    krate: [f64; N_CLASSES],
    sum_d: f64,
    sum_c: f64, // without the factor C
    pub c: f64,
    pub t: f64, // seconds
    pub events: u64,
    pub n_crystal: u64,
    pub hit_wall: bool,
    /// None: exact continuous-time kMC. Some((dt, random_order)): the paper's synchronous sweep.
    pub sync: Option<(f64, bool)>,
    swap12: bool,
    /// 0: A on the molecule's NH2 side (Fig. 2 as drawn); 1: A on its C=O side (default); 2: A = +c for every molecule
    ab: u8,
    split_on: [[bool; 9]; 4],
    rng: Rng,
}

pub struct Stats {
    pub n: u64,
    pub coord: [u64; 7],
    pub lc_nm: f64,
    pub w110_nm: f64,
    pub w1m10_nm: f64,
}

impl Sim {
    pub fn new(lat: Lattice, rates: &Rates, c: f64, swap12: bool, ab: u8, seed: u64) -> Sim {
        let n = lat.n_sites();
        assert!(n < u32::MAX as usize, "box too large for u32 site ids");
        let mut krate = [0.0; N_CLASSES];
        for kind in 0..2 {
            for r in 0..4 {
                for cc in 0..9 {
                    for s in 0..3 {
                        let id = ((kind * 4 + r) * 9 + cc) * 3 + s;
                        krate[id] = if kind == 0 { rates.kd[r][cc][s] } else { rates.kc[r][cc][s] };
                    }
                }
            }
        }
        Sim {
            lat,
            occ: vec![0; n],
            cls: vec![NONE; n],
            pos: vec![0; n],
            members: vec![Vec::new(); N_CLASSES],
            krate,
            sum_d: 0.0,
            sum_c: 0.0,
            c,
            t: 0.0,
            events: 0,
            n_crystal: 0,
            hit_wall: false,
            sync: None,
            swap12,
            ab,
            split_on: rates.split_on,
            rng: Rng::new(seed),
        }
    }

    /// Fill a block of w x w x h unit cells (both sublattices) centered in the box, then classify all sites.
    pub fn seed_block(&mut self, w: usize, h: usize) {
        let l = self.lat;
        assert!(w + 4 < l.nx && w + 4 < l.ny && h + 4 < l.nz, "nucleus does not fit in the box");
        let (i0, j0, k0) = ((l.nx - w) / 2, (l.ny - w) / 2, (l.nz - h) / 2);
        for k in k0..k0 + h {
            for j in j0..j0 + w {
                for i in i0..i0 + w {
                    for s in 0..2 {
                        self.occ[l.idx(s, i, j, k)] = 1;
                        self.n_crystal += 1;
                    }
                }
            }
        }
        for site in 0..l.n_sites() {
            let c = self.classify(site);
            self.set_class(site, c);
        }
    }

    /// Slab geometry: fill (110) layers m0 .. m0 + w completely, then classify all sites.
    pub fn seed_layers(&mut self, m0: usize, w: usize) {
        let l = self.lat;
        assert!(l.slab && m0 >= 2 && m0 + w + 2 <= l.ny, "layers do not fit in the slab");
        for k in 0..l.nz {
            for m in m0..m0 + w {
                for n in 0..l.nx {
                    for s in 0..2 {
                        self.occ[l.idx(s, n, m, k)] = 1;
                        self.n_crystal += 1;
                    }
                }
            }
        }
        for site in 0..l.n_sites() {
            let c = self.classify(site);
            self.set_class(site, c);
        }
    }

    /// Slab geometry: for every column (n, k), the highest and lowest occupied layer (-1 if empty).
    pub fn face_heights(&self) -> (Vec<i32>, Vec<i32>) {
        let l = self.lat;
        let cols = l.nx * l.nz;
        let mut top = vec![-1i32; cols];
        let mut bot = vec![i32::MAX; cols];
        for site in 0..self.occ.len() {
            if self.occ[site] == 0 {
                continue;
            }
            let (_, n, m, k) = l.coords(site);
            let c = k * l.nx + n;
            top[c] = top[c].max(m as i32);
            bot[c] = bot[c].min(m as i32);
        }
        for b in bot.iter_mut() {
            if *b == i32::MAX {
                *b = -1;
            }
        }
        (top, bot)
    }

    fn classify(&self, site: usize) -> u8 {
        let l = &self.lat;
        let (s, i, j, k) = l.coords(site);
        if !l.interior(i, j, k) {
            return NONE;
        }
        let (mut n1, mut n2, mut a, mut b) = (0usize, 0usize, false, false);
        for (m, o) in NN[s].iter().enumerate() {
            if self.occ[l.shift(i, j, k, *o)] != 0 {
                match m {
                    0 | 1 => n1 += 1,
                    2 | 3 => n2 += 1,
                    4 => a = true,
                    _ => b = true,
                }
            }
        }
        if self.swap12 {
            std::mem::swap(&mut n1, &mut n2);
        }
        // NN lists A on the molecule's NH2 side. alpha's NH2 side is -c, beta's is +c, so "A = +c for
        // every molecule" swaps alpha only.
        if self.ab == 1 || (self.ab == 2 && s == 0) {
            std::mem::swap(&mut a, &mut b);
        }
        let occupied = self.occ[site] != 0;
        if !occupied && n1 + n2 + a as usize + b as usize == 0 {
            return NONE; // not a surface site: no homogeneous nucleation
        }
        let (r, c) = (row(a, b), col(n1, n2));
        let sp = if self.split_on[r][c] {
            let mut nz = 0;
            let mut nxy = 0;
            for (m, o) in NNN[s].iter().enumerate() {
                if self.occ[l.shift(i, j, k, *o)] != 0 {
                    if m < 4 { nz += 1 } else { nxy += 1 }
                }
            }
            if nz > nxy { 0 } else if nz == nxy { 1 } else { 2 }
        } else {
            0
        };
        let kind = if occupied { 0 } else { 1 };
        let id = ((kind * 4 + r) * 9 + c) * 3 + sp;
        if self.krate[id] > 0.0 { id as u8 } else { NONE }
    }

    fn set_class(&mut self, site: usize, new: u8) {
        let old = self.cls[site];
        if old == new {
            return;
        }
        if old != NONE {
            let o = old as usize;
            let p = self.pos[site] as usize;
            let v = &mut self.members[o];
            let last = v.pop().unwrap();
            if p < v.len() {
                v[p] = last;
                self.pos[last as usize] = p as u32;
            }
            if o < HALF { self.sum_d -= self.krate[o] } else { self.sum_c -= self.krate[o] }
        }
        if new != NONE {
            let nn = new as usize;
            self.pos[site] = self.members[nn].len() as u32;
            self.members[nn].push(site as u32);
            if nn < HALF { self.sum_d += self.krate[nn] } else { self.sum_c += self.krate[nn] }
        }
        self.cls[site] = new;
    }

    fn resum(&mut self) {
        self.sum_d = (0..HALF).map(|q| self.members[q].len() as f64 * self.krate[q]).sum();
        self.sum_c = (HALF..N_CLASSES).map(|q| self.members[q].len() as f64 * self.krate[q]).sum();
    }

    pub fn total_rate(&self) -> f64 {
        self.sum_d + self.c * self.sum_c
    }

    /// One event. Returns false if nothing can happen (or growth reached the box wall).
    pub fn step(&mut self) -> bool {
        if self.events & 0xF_FFFF == 0 {
            self.resum();
        }
        let total = self.total_rate();
        if total <= 0.0 || self.hit_wall {
            return false;
        }
        let mut u = self.rng.f64() * total;
        let (range, scale) = if u < self.sum_d {
            (0..HALF, 1.0)
        } else {
            u -= self.sum_d;
            (HALF..N_CLASSES, self.c)
        };
        let mut pick = NONE as usize;
        for q in range {
            let len = self.members[q].len();
            if len == 0 {
                continue;
            }
            pick = q;
            let w = len as f64 * self.krate[q] * scale;
            if u < w {
                break;
            }
            u -= w;
        }
        if pick == NONE as usize {
            return false;
        }
        let m = &self.members[pick];
        let site = m[self.rng.below(m.len())] as usize;
        self.apply(site);
        self.events += 1;
        self.t += -(1.0 - self.rng.f64()).ln() / total;
        true
    }

    /// Flip one site (crystal <-> solution) and reclassify it and every site whose rate depends on it.
    pub fn apply(&mut self, site: usize) {
        let l = self.lat;
        let (s, i, j, k) = l.coords(site);
        if self.occ[site] == 0 {
            self.occ[site] = 1;
            self.n_crystal += 1;
            if l.near_wall(i, j, k) {
                self.hit_wall = true;
            }
        } else {
            self.occ[site] = 0;
            self.n_crystal -= 1;
        }
        let c = self.classify(site);
        self.set_class(site, c);
        for o in NN[s].iter() {
            let y = l.shift(i, j, k, *o);
            let c = self.classify(y);
            self.set_class(y, c);
        }
        for o in NNN[s].iter() {
            let y = l.shift(i, j, k, *o);
            let c = self.classify(y);
            self.set_class(y, c);
        }
    }

    /// The paper's scheme (SI): every surface site, occupied or not, is visited once per interval dt and
    /// changes with probability rate * dt, capped at 1; then the clock advances by dt. Sites are visited
    /// in place ("sequentially"), so a change is seen by sites visited later in the same sweep. Sites that
    /// become surface sites during a sweep wait for the next one.
    /// `random_order`: shuffle the visiting order each sweep; otherwise lattice order (bottom layer first).
    /// Returns false if nothing is left to do.
    pub fn sweep(&mut self, dt: f64, random_order: bool) -> bool {
        if self.hit_wall || self.total_rate() <= 0.0 {
            return false;
        }
        let mut list: Vec<u32> = self.members.iter().flatten().copied().collect();
        if random_order {
            for i in (1..list.len()).rev() {
                let j = self.rng.below(i + 1);
                list.swap(i, j);
            }
        } else {
            list.sort_unstable();
        }
        for &site in &list {
            let site = site as usize;
            let q = self.cls[site];
            if q == NONE {
                continue; // no longer a surface site, or its rate became zero
            }
            let q = q as usize;
            let rate = if q < HALF { self.krate[q] } else { self.krate[q] * self.c };
            let p = (rate * dt).min(1.0);
            if self.rng.f64() < p {
                self.apply(site);
                self.events += 1;
                if self.hit_wall {
                    break;
                }
            }
        }
        self.t += dt;
        true
    }

    /// Run until simulated time t_end (seconds) or until stopped.
    pub fn run_until(&mut self, t_end: f64) {
        while self.t < t_end {
            let more = match self.sync {
                None => self.step(),
                Some((dt, random_order)) => self.sweep(dt, random_order),
            };
            if !more {
                break;
            }
        }
    }

    pub fn coordination(&self, site: usize) -> usize {
        let (s, i, j, k) = self.lat.coords(site);
        if !self.lat.interior(i, j, k) {
            return 0;
        }
        NN[s].iter().filter(|o| self.occ[self.lat.shift(i, j, k, **o)] != 0).count()
    }

    pub fn stats(&self) -> Stats {
        let mut coord = [0u64; 7];
        let (mut zmin, mut zmax) = (f64::MAX, f64::MIN);
        let (mut umin, mut umax, mut vmin, mut vmax) = (f64::MAX, f64::MIN, f64::MAX, f64::MIN);
        let r2 = std::f64::consts::FRAC_1_SQRT_2;
        for site in 0..self.occ.len() {
            if self.occ[site] == 0 {
                continue;
            }
            coord[self.coordination(site)] += 1;
            let (x, y, z) = self.lat.position_nm(site);
            let (u, v) = ((x + y) * r2, (x - y) * r2);
            zmin = zmin.min(z);
            zmax = zmax.max(z);
            umin = umin.min(u);
            umax = umax.max(u);
            vmin = vmin.min(v);
            vmax = vmax.max(v);
        }
        let n: u64 = coord.iter().sum();
        if n == 0 {
            return Stats { n, coord, lc_nm: 0.0, w110_nm: 0.0, w1m10_nm: 0.0 };
        }
        Stats { n, coord, lc_nm: zmax - zmin, w110_nm: umax - umin, w1m10_nm: vmax - vmin }
    }

    /// Crystal molecules with fewer than six neighbors: x, y, z (nm), sublattice, coordination.
    pub fn surface_csv(&self) -> String {
        let mut out = String::from("x_nm,y_nm,z_nm,sublattice,n\n");
        for site in 0..self.occ.len() {
            if self.occ[site] == 0 {
                continue;
            }
            let n = self.coordination(site);
            if n < 6 {
                let (x, y, z) = self.lat.position_nm(site);
                out.push_str(&format!("{x:.4},{y:.4},{z:.4},{},{n}\n", site & 1));
            }
        }
        out
    }

    /// Debug check: incremental bookkeeping equals a full reclassification.
    #[cfg(test)]
    pub fn consistent(&self) -> bool {
        (0..self.occ.len()).all(|s| self.classify(s) == self.cls[s])
            && self.members.iter().enumerate().all(|(q, v)| {
                v.iter().enumerate().all(|(p, &s)| self.cls[s as usize] as usize == q && self.pos[s as usize] as usize == p)
            })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::rates::{Rates, Solvent, Split, Unobserved};

    #[test]
    fn slab_bookkeeping_survives_many_events() {
        let rates = Rates::new(Solvent::Water, Split::Nnn, Unobserved::Fill);
        let mut sim = Sim::new(Lattice::slab(12, 16, 10), &rates, 0.3, false, 1, 3);
        sim.seed_layers(6, 4);
        assert_eq!(sim.n_crystal, 2 * 12 * 4 * 10);
        for _ in 0..20_000 {
            if !sim.step() {
                break;
            }
        }
        assert!(sim.consistent());
        let n = sim.occ.iter().filter(|&&o| o != 0).count() as u64;
        assert_eq!(n, sim.n_crystal);
    }

    #[test]
    fn bookkeeping_survives_many_events() {
        for swap in [false, true] {
            let rates = Rates::new(Solvent::Methanol, Split::Nnn, Unobserved::Fill);
            let lat = Lattice::boxed(16, 16, 16);
            let mut sim = Sim::new(lat, &rates, 0.2, swap, if swap { 2 } else { 1 }, 7);
            sim.seed_block(6, 6);
            for _ in 0..20_000 {
                if !sim.step() {
                    break;
                }
            }
            assert!(sim.consistent());
            let n = sim.occ.iter().filter(|&&o| o != 0).count() as u64;
            assert_eq!(n, sim.n_crystal);
            for _ in 0..50 {
                if !sim.sweep(50e-12, swap) {
                    break;
                }
            }
            assert!(sim.consistent());
            let n = sim.occ.iter().filter(|&&o| o != 0).count() as u64;
            assert_eq!(n, sim.n_crystal);
        }
    }
}
