//! Rate lookup from the decoded SI tables (or a rate file in the same JSON schema, `--rates`), plus the
//! modeling switches the paper leaves open.

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
    /// The compiled SI tables.
    pub fn new(solvent: Solvent, split: Split, unobserved: Unobserved) -> Rates {
        let (kd, kc) = match solvent {
            Solvent::Water => (WATER_DISSOLUTION, WATER_CRYSTALLIZATION),
            Solvent::Methanol => (METHANOL_DISSOLUTION, METHANOL_CRYSTALLIZATION),
        };
        Rates::from_tables(kd, kc, SPLIT, split, unobserved)
    }

    /// A rate file in the schema of data/piana_gale_2005_rates.json (e.g. a predicted table). The same
    /// split and unobserved-cell rules apply as for the SI tables.
    pub fn from_file(path: &str, solvent: Solvent, split: Split, unobserved: Unobserved) -> Result<Rates, String> {
        let text = std::fs::read_to_string(path).map_err(|e| format!("{path}: {e}"))?;
        let (kd, kc, mask) = parse_tables(&text, solvent).map_err(|e| format!("{path}: {e}"))?;
        Ok(Rates::from_tables(kd, kc, mask, split, unobserved))
    }

    fn from_tables(mut kd: Table, mut kc: Table, mut split_on: [[bool; 9]; 4], split: Split, unobserved: Unobserved) -> Rates {
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

// ---- reading a rate file: just enough JSON for {"<table>": {"<col>/<row>": [numbers], ...}, ...} ----

enum Json {
    Num(f64),
    Str(String),
    Arr(Vec<Json>),
    Obj(Vec<(String, Json)>),
}

struct Parser<'a> {
    b: &'a [u8],
    i: usize,
}

impl Parser<'_> {
    fn ws(&mut self) {
        while self.i < self.b.len() && self.b[self.i].is_ascii_whitespace() {
            self.i += 1;
        }
    }
    fn eat(&mut self, c: u8) -> Result<(), String> {
        self.ws();
        if self.b.get(self.i) == Some(&c) {
            self.i += 1;
            Ok(())
        } else {
            Err(format!("expected '{}' at byte {}", c as char, self.i))
        }
    }
    fn value(&mut self) -> Result<Json, String> {
        self.ws();
        match self.b.get(self.i) {
            Some(b'{') => {
                self.i += 1;
                let mut v = Vec::new();
                self.ws();
                if self.b.get(self.i) == Some(&b'}') {
                    self.i += 1;
                    return Ok(Json::Obj(v));
                }
                loop {
                    let k = match self.value()? {
                        Json::Str(s) => s,
                        _ => return Err(format!("object key must be a string, byte {}", self.i)),
                    };
                    self.eat(b':')?;
                    v.push((k, self.value()?));
                    self.ws();
                    match self.b.get(self.i) {
                        Some(b',') => self.i += 1,
                        Some(b'}') => {
                            self.i += 1;
                            return Ok(Json::Obj(v));
                        }
                        _ => return Err(format!("expected ',' or '}}' at byte {}", self.i)),
                    }
                }
            }
            Some(b'[') => {
                self.i += 1;
                let mut v = Vec::new();
                self.ws();
                if self.b.get(self.i) == Some(&b']') {
                    self.i += 1;
                    return Ok(Json::Arr(v));
                }
                loop {
                    v.push(self.value()?);
                    self.ws();
                    match self.b.get(self.i) {
                        Some(b',') => self.i += 1,
                        Some(b']') => {
                            self.i += 1;
                            return Ok(Json::Arr(v));
                        }
                        _ => return Err(format!("expected ',' or ']' at byte {}", self.i)),
                    }
                }
            }
            Some(b'"') => {
                self.i += 1;
                let start = self.i;
                while self.i < self.b.len() && self.b[self.i] != b'"' {
                    if self.b[self.i] == b'\\' {
                        return Err(format!("escapes in strings are not supported, byte {}", self.i));
                    }
                    self.i += 1;
                }
                let s = std::str::from_utf8(&self.b[start..self.i]).map_err(|e| e.to_string())?.to_string();
                self.eat(b'"')?;
                Ok(Json::Str(s))
            }
            Some(_) => {
                let start = self.i;
                while self.i < self.b.len() && matches!(self.b[self.i], b'0'..=b'9' | b'-' | b'+' | b'.' | b'e' | b'E') {
                    self.i += 1;
                }
                let t = std::str::from_utf8(&self.b[start..self.i]).unwrap();
                t.parse::<f64>().map(Json::Num).map_err(|_| format!("bad number '{t}' at byte {start}"))
            }
            None => Err("unexpected end of file".into()),
        }
    }
}

/// Dissolution table, crystallization table and split mask (cells with three values) for one solvent.
fn parse_tables(text: &str, solvent: Solvent) -> Result<(Table, Table, [[bool; 9]; 4]), String> {
    let mut p = Parser { b: text.as_bytes(), i: 0 };
    let root = match p.value()? {
        Json::Obj(v) => v,
        _ => return Err("top level must be an object".into()),
    };
    let s = match solvent {
        Solvent::Water => "water",
        Solvent::Methanol => "methanol",
    };
    let mut mask = [[false; 9]; 4];
    let read = |kind: &str, mask: &mut [[bool; 9]; 4]| -> Result<Table, String> {
        let name = format!("{s}_{kind}");
        let cells = match root.iter().find(|(k, _)| *k == name) {
            Some((_, Json::Obj(v))) => v,
            Some(_) => return Err(format!("{name} must be an object")),
            None => return Err(format!("no table {name}")),
        };
        let mut t: Table = [[[0.0; 3]; 9]; 4];
        for r in 0..4 {
            for c in 0..9 {
                let key = format!("{}/{}", COL_NAMES[c], ROW_NAMES[r]);
                let vals: Vec<f64> = match cells.iter().find(|(k, _)| *k == key) {
                    Some((_, Json::Arr(a))) => a
                        .iter()
                        .map(|x| match x {
                            Json::Num(n) if *n >= 0.0 => Ok(*n),
                            _ => Err(format!("{name} {key}: rates must be non-negative numbers")),
                        })
                        .collect::<Result<_, _>>()?,
                    _ => return Err(format!("{name}: missing cell {key}")),
                };
                t[r][c] = match vals.len() {
                    1 => [vals[0]; 3],
                    3 => {
                        mask[r][c] = true;
                        [vals[0], vals[1], vals[2]]
                    }
                    n => return Err(format!("{name} {key}: {n} values (1 or 3 expected)")),
                };
            }
        }
        Ok(t)
    };
    let kd = read("dissolution", &mut mask)?;
    let kc = read("crystallization", &mut mask)?;
    Ok((kd, kc, mask))
}

#[cfg(test)]
mod tests {
    use super::*;

    const SI_JSON: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/data/piana_gale_2005_rates.json");

    #[test]
    fn rate_file_reproduces_the_compiled_tables() {
        // --rates data/piana_gale_2005_rates.json must be the default model, value for value.
        let text = std::fs::read_to_string(SI_JSON).unwrap();
        for (s, kd0, kc0) in [
            (Solvent::Water, WATER_DISSOLUTION, WATER_CRYSTALLIZATION),
            (Solvent::Methanol, METHANOL_DISSOLUTION, METHANOL_CRYSTALLIZATION),
        ] {
            let (kd, kc, mask) = parse_tables(&text, s).unwrap();
            for r in 0..4 {
                for c in 0..9 {
                    for k in 0..3 {
                        assert_eq!(kd[r][c][k].to_bits(), kd0[r][c][k].to_bits(), "{s:?} kd {r} {c} {k}");
                        assert_eq!(kc[r][c][k].to_bits(), kc0[r][c][k].to_bits(), "{s:?} kc {r} {c} {k}");
                    }
                }
            }
            assert_eq!(mask, SPLIT);
            for (sp, un) in [(Split::Nnn, Unobserved::Fill), (Split::Mean, Unobserved::Published)] {
                let a = Rates::new(s, sp, un);
                let b = Rates::from_file(SI_JSON, s, sp, un).unwrap();
                assert_eq!(a.kd, b.kd);
                assert_eq!(a.kc, b.kc);
                assert_eq!(a.split_on, b.split_on);
                assert_eq!(a.filled, b.filled);
            }
        }
    }

    #[test]
    fn rate_file_errors_name_the_problem() {
        let text = std::fs::read_to_string(SI_JSON).unwrap();
        let missing = text.replacen("\"12/B\"", "\"12/X\"", 1);
        let e = parse_tables(&missing, Solvent::Water).err().unwrap();
        assert!(e.contains("missing cell 12/B"), "{e}");
        let negative = text.replacen("18790000000000.0", "-1.0", 1);
        assert!(parse_tables(&negative, Solvent::Water).is_err());
        assert!(parse_tables("{}", Solvent::Methanol).err().unwrap().contains("no table methanol_dissolution"));
    }

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
