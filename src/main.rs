//! urea-kmc: 3D kinetic Monte Carlo of urea crystal growth from solution, using the MD-derived
//! site rates of Piana, Reyhani & Gale, Nature 438, 70 (2005).
//!
//!   urea-kmc scan --solvent water             net growth vs concentration -> saturation estimate
//!   urea-kmc grow --solvent water --cs 0.22   one growth run -> morphology log and surface dump
//!
//! Run `urea-kmc help` for every option.

mod lattice;
mod rates;
mod rates_table;
mod rng;
mod sim;

use lattice::{crystal_molarity, Lattice};
use rates::{Rates, Solvent, Split, Unobserved};
use sim::Sim;
use std::collections::HashMap;
use std::io::Write;
use std::time::Instant;

const HELP: &str = "\
urea-kmc scan|grow [options]

common
  --solvent water|methanol      (water)
  --split nnn|mean              second-neighbour split as in the SI, or its mean (nnn)
  --unobserved fill|published   cells never seen in MD: fill from same-coordination mean, or keep zeros (fill)
  --swap12                      swap which in-plane pair is called 1 and which 2
  --ab co|nh2|crystal           which chain neighbour the SI calls A: the molecule's C=O side (co, default: the
                                only reading that reproduces the paper's solubilities), its NH2 side (nh2, Fig. 2
                                as drawn), or the one above (+c) for every molecule (crystal: breaks the symmetry)
  --swapAB                      same as --ab co (kept for old commands)
  --nucleus W | WxH             starting block in unit cells (16)
  --box N | NXxNYxNZ            box in unit cells (scan 48, grow 48x48x400)
  --seed N                      (1)
  --scheme exact|sync           exact continuous-time kMC, or the paper's sweep: every surface site gets a
                                chance rate*dt (capped at 1) once per dt (exact)
  --dt-ps T                     sweep interval for --scheme sync (50, as in the paper)
  --order random|lattice        visiting order within a sweep (random)

scan
  --c-list a,b,c                concentrations as a fraction of lattice sites (0.05,0.1,...,0.5)
  --t-ns T                      simulated time per concentration (50)

grow
  --c X                         concentration (fraction of sites), or:
  --cs X --sigma S              C = cs * (1 + S); S is the paper's supersaturation C* (0.03)
  --t-ns T                      stop after T ns of simulated time (10000)
  --report-ns R                 log interval (t-ns / 40)
  --out DIR                     writes DIR/log.csv and DIR/surface.csv (runs/<solvent>)

face (one [110] side face, periodic along [1-10] and c, grown from whole molecular layers)
  --c X | --cs X --sigma S      concentration, as for grow
  --face-nm L                   face size, L x L nm (200)
  --layers W                    starting slab thickness in (110) layers (20)
  --margin M                    empty layers on each side; reaching them ends the run (20)
  --t-ns T, --report-ns R       as for grow (100, T/20)
  --out DIR                     DIR/log.csv and DIR/top.csv, the top face's height map (runs/face-<solvent>)

concentration: 1.0 = every site full = crystal density (22.89 M)";

fn parse_args() -> (String, HashMap<String, String>) {
    let mut args = std::env::args().skip(1);
    let cmd = args.next().unwrap_or_else(|| "help".into());
    let mut map = HashMap::new();
    let mut key: Option<String> = None;
    for a in args {
        if a.trim().is_empty() {
            continue;
        }
        if let Some(k) = a.strip_prefix("--") {
            if let Some(prev) = key.take() {
                map.insert(prev, "true".into());
            }
            key = Some(k.to_string());
        } else if let Some(k) = key.take() {
            map.insert(k, a);
        } else {
            die(&format!("unexpected argument {a}"));
        }
    }
    if let Some(k) = key {
        map.insert(k, "true".into());
    }
    (cmd, map)
}

fn die(msg: &str) -> ! {
    eprintln!("error: {msg}\n\n{HELP}");
    std::process::exit(2)
}

fn get<'a>(m: &'a HashMap<String, String>, k: &str, default: &'a str) -> &'a str {
    m.get(k).map(|s| s.as_str()).unwrap_or(default)
}

fn num(m: &HashMap<String, String>, k: &str, default: &str) -> f64 {
    get(m, k, default).parse().unwrap_or_else(|_| die(&format!("--{k} needs a number")))
}

fn dims3(s: &str) -> (usize, usize, usize) {
    let v: Vec<usize> = s.split('x').map(|p| p.parse().unwrap_or_else(|_| die(&format!("bad size {s}")))).collect();
    match v.as_slice() {
        [n] => (*n, *n, *n),
        [a, b, c] => (*a, *b, *c),
        _ => die(&format!("bad size {s}")),
    }
}

fn dims2(s: &str) -> (usize, usize) {
    let v: Vec<usize> = s.split('x').map(|p| p.parse().unwrap_or_else(|_| die(&format!("bad size {s}")))).collect();
    match v.as_slice() {
        [n] => (*n, *n),
        [w, h] => (*w, *h),
        _ => die(&format!("bad nucleus {s}")),
    }
}

struct Model {
    solvent: Solvent,
    rates: Rates,
    swap12: bool,
    ab: u8,
    sync: Option<(f64, bool)>,
    label: String,
}

fn model(m: &HashMap<String, String>) -> Model {
    let solvent = match get(m, "solvent", "water") {
        "water" => Solvent::Water,
        "methanol" => Solvent::Methanol,
        s => die(&format!("unknown solvent {s}")),
    };
    let split = match get(m, "split", "nnn") {
        "nnn" => Split::Nnn,
        "mean" => Split::Mean,
        s => die(&format!("unknown split {s}")),
    };
    let unobserved = match get(m, "unobserved", "fill") {
        "fill" => Unobserved::Fill,
        "published" => Unobserved::Published,
        s => die(&format!("unknown unobserved {s}")),
    };
    let swap12 = m.contains_key("swap12");
    let ab_name = if m.contains_key("swapAB") { "co" } else { get(m, "ab", "co") };
    let ab = match ab_name {
        "nh2" => 0,
        "co" => 1,
        "crystal" => 2,
        s => die(&format!("unknown --ab {s}")),
    };
    let rates = Rates::new(solvent, split, unobserved);
    let sync = match get(m, "scheme", "exact") {
        "exact" => None,
        "sync" => {
            let dt = num(m, "dt-ps", "50") * 1e-12;
            let random = match get(m, "order", "random") {
                "random" => true,
                "lattice" => false,
                s => die(&format!("unknown --order {s}")),
            };
            Some((dt, random))
        }
        s => die(&format!("unknown --scheme {s}")),
    };
    let scheme_label = match sync {
        None => "exact".to_string(),
        Some((dt, r)) => format!("sync dt={}ps order={}", dt * 1e12, if r { "random" } else { "lattice" }),
    };
    let label = format!(
        "solvent={solvent:?} split={split:?} unobserved={unobserved:?} swap12={swap12} ab={ab_name} scheme={scheme_label}"
    );
    Model { solvent, rates, swap12, ab, sync, label }
}

fn header(md: &Model) {
    println!("# urea-kmc · {}", md.label);
    for f in &md.rates.filled {
        println!("# filled unobserved cell {f}");
    }
}

fn scan(m: &HashMap<String, String>) {
    let md = model(m);
    header(&md);
    let (nx, ny, nz) = dims3(get(m, "box", "48"));
    let (w, h) = dims2(get(m, "nucleus", "16"));
    let t_ns = num(m, "t-ns", "50");
    let seed = num(m, "seed", "1") as u64;
    let list: Vec<f64> = match m.get("c-list") {
        Some(s) => s.split(',').map(|p| p.parse().unwrap_or_else(|_| die("bad --c-list"))).collect(),
        None => (1..=10).map(|i| 0.05 * i as f64).collect(),
    };
    let molar = crystal_molarity();
    println!("# box {nx}x{ny}x{nz} cells, nucleus {w}x{w}x{h} cells, {t_ns} ns per point");
    println!("c,c_molar,t_ns,events,n_start,n_end,net_per_ns,hit_wall,seconds");
    let mut pts = Vec::new();
    for (q, &c) in list.iter().enumerate() {
        let wall = Instant::now();
        let mut sim = Sim::new(Lattice::boxed(nx, ny, nz), &md.rates, c, md.swap12, md.ab, seed + q as u64);
        sim.sync = md.sync;
        sim.seed_block(w, h);
        let n0 = sim.n_crystal;
        sim.run_until(t_ns * 1e-9);
        let t = sim.t * 1e9;
        let net = (sim.n_crystal as f64 - n0 as f64) / t.max(1e-12);
        println!(
            "{c:.4},{:.2},{t:.2},{},{n0},{},{net:.3},{},{:.1}",
            c * molar, sim.events, sim.n_crystal, sim.hit_wall, wall.elapsed().as_secs_f64()
        );
        pts.push((c, net));
    }
    for p in pts.windows(2) {
        let ((c0, r0), (c1, r1)) = (p[0], p[1]);
        if r0 < 0.0 && r1 >= 0.0 {
            let cs = c0 + (c1 - c0) * (-r0) / (r1 - r0);
            println!(
                "# net growth changes sign near C = {cs:.4} ({:.2} M); paper (2D screw dislocation on (110)): {} M",
                cs * molar,
                if md.solvent == Solvent::Water { "5" } else { "4" }
            );
        }
    }
}

fn grow(m: &HashMap<String, String>) {
    let md = model(m);
    header(&md);
    let (nx, ny, nz) = dims3(get(m, "box", "48x48x400"));
    let (w, h) = dims2(get(m, "nucleus", "16"));
    let c = if let Some(c) = m.get("c") {
        c.parse().unwrap_or_else(|_| die("--c needs a number"))
    } else if m.contains_key("cs") {
        num(m, "cs", "0") * (1.0 + num(m, "sigma", "0.03"))
    } else {
        die("grow needs --c, or --cs (with optional --sigma)")
    };
    let t_ns = num(m, "t-ns", "10000");
    let report = num(m, "report-ns", &format!("{}", t_ns / 40.0));
    let seed = num(m, "seed", "1") as u64;
    let solvent = if md.solvent == Solvent::Water { "water" } else { "methanol" };
    let out = get(m, "out", &format!("runs/{solvent}")).to_string();
    std::fs::create_dir_all(&out).unwrap_or_else(|e| die(&format!("cannot create {out}: {e}")));
    let mut log = std::fs::File::create(format!("{out}/log.csv")).unwrap();
    writeln!(log, "# {}", md.label).unwrap();
    writeln!(log, "# C = {c:.5} ({:.3} M), box {nx}x{ny}x{nz}, nucleus {w}x{w}x{h}", c * crystal_molarity()).unwrap();
    let cols = "t_ns,events,n,length_c_nm,width_110_nm,width_1m10_nm,aspect,n3,n4,n5,seconds";
    writeln!(log, "{cols}").unwrap();
    println!("# C = {c:.5} ({:.3} M), box {nx}x{ny}x{nz} cells, nucleus {w}x{w}x{h}", c * crystal_molarity());
    println!("{cols}");

    let wall = Instant::now();
    let mut sim = Sim::new(Lattice::boxed(nx, ny, nz), &md.rates, c, md.swap12, md.ab, seed);
    sim.sync = md.sync;
    sim.seed_block(w, h);
    let mut next = 0.0;
    loop {
        sim.run_until(next * 1e-9);
        let st = sim.stats();
        let width = 0.5 * (st.w110_nm + st.w1m10_nm);
        let line = format!(
            "{:.1},{},{},{:.2},{:.2},{:.2},{:.2},{},{},{},{:.1}",
            sim.t * 1e9, sim.events, st.n, st.lc_nm, st.w110_nm, st.w1m10_nm,
            if width > 0.0 { st.lc_nm / width } else { 0.0 },
            st.coord[3], st.coord[4], st.coord[5], wall.elapsed().as_secs_f64()
        );
        println!("{line}");
        writeln!(log, "{line}").unwrap();
        if sim.hit_wall || st.n == 0 || sim.t * 1e9 >= t_ns || sim.total_rate() == 0.0 {
            break;
        }
        next = (next + report).min(t_ns);
    }
    let why = if sim.hit_wall {
        "crystal reached the box wall"
    } else if sim.n_crystal == 0 {
        "crystal dissolved"
    } else {
        "time limit"
    };
    println!("# stopped: {why}");
    std::fs::write(format!("{out}/surface.csv"), sim.surface_csv()).unwrap();
    println!("# wrote {out}/log.csv and {out}/surface.csv");
}

fn face(m: &HashMap<String, String>) {
    let md = model(m);
    header(&md);
    let c = if let Some(c) = m.get("c") {
        c.parse().unwrap_or_else(|_| die("--c needs a number"))
    } else if m.contains_key("cs") {
        num(m, "cs", "0") * (1.0 + num(m, "sigma", "0"))
    } else {
        die("face needs --c, or --cs (with optional --sigma)")
    };
    let side = num(m, "face-nm", "200");
    let period = (side / (2.0 * lattice::d110_nm())).round().max(4.0) as usize; // 0.787 nm per step along [1-10]
    let nz = (side / (lattice::C_LEN / 10.0)).round().max(4.0) as usize;
    let w = num(m, "layers", "20") as usize;
    let margin = num(m, "margin", "20") as usize;
    let t_ns = num(m, "t-ns", "100");
    let report = num(m, "report-ns", &format!("{}", t_ns / 20.0));
    let seed = num(m, "seed", "1") as u64;
    let solvent = if md.solvent == Solvent::Water { "water" } else { "methanol" };
    let out = get(m, "out", &format!("runs/face-{solvent}")).to_string();
    std::fs::create_dir_all(&out).unwrap_or_else(|e| die(&format!("cannot create {out}: {e}")));

    let lat = Lattice::slab(period, w + 2 * margin, nz);
    let along = period as f64 * 2.0 * lattice::d110_nm();
    let up = nz as f64 * lattice::C_LEN / 10.0;
    let per_layer = (2 * period * nz) as f64; // molecules in one (110) layer
    let mut log = std::fs::File::create(format!("{out}/log.csv")).unwrap();
    let head = format!(
        "# C = {c:.5} ({:.3} M), face {along:.0} nm along [1-10] x {up:.0} nm along c, {w} layers + {margin} margin each side, {} sites",
        c * crystal_molarity(), lat.n_sites()
    );
    writeln!(log, "# {}", md.label).unwrap();
    writeln!(log, "{head}").unwrap();
    println!("{head}");
    let cols = "t_ns,events,n,layers_grown_per_face,top_mean,top_rms,top_min,top_max,n3,n4,n5,seconds";
    writeln!(log, "{cols}").unwrap();
    println!("{cols}");

    let wall = Instant::now();
    let mut sim = Sim::new(lat, &md.rates, c, md.swap12, md.ab, seed);
    sim.sync = md.sync;
    sim.seed_layers(margin, w);
    let n0 = sim.n_crystal as f64;
    let top0 = (margin + w - 1) as f64;
    let mut next = 0.0;
    loop {
        sim.run_until(next * 1e-9);
        let (top, _) = sim.face_heights();
        let h: Vec<f64> = top.iter().filter(|&&t| t >= 0).map(|&t| t as f64 - top0).collect();
        let mean = h.iter().sum::<f64>() / h.len().max(1) as f64;
        let rms = (h.iter().map(|x| (x - mean).powi(2)).sum::<f64>() / h.len().max(1) as f64).sqrt();
        let (lo, hi) = h.iter().fold((f64::MAX, f64::MIN), |(a, b), &x| (a.min(x), b.max(x)));
        let mut coord = [0u64; 7];
        for site in 0..sim.occ.len() {
            if sim.occ[site] != 0 {
                coord[sim.coordination(site)] += 1;
            }
        }
        let grown = (sim.n_crystal as f64 - n0) / (2.0 * per_layer);
        let line = format!(
            "{:.1},{},{},{grown:.3},{mean:.2},{rms:.2},{lo:.0},{hi:.0},{},{},{},{:.1}",
            sim.t * 1e9, sim.events, sim.n_crystal, coord[3], coord[4], coord[5], wall.elapsed().as_secs_f64()
        );
        println!("{line}");
        writeln!(log, "{line}").unwrap();
        if sim.hit_wall || sim.n_crystal == 0 || sim.t * 1e9 >= t_ns || sim.total_rate() == 0.0 {
            break;
        }
        next = (next + report).min(t_ns);
    }
    let grown = (sim.n_crystal as f64 - n0) / (2.0 * per_layer);
    let t_us = sim.t * 1e6;
    println!(
        "# each face grew {grown:.3} layers = {:.3} nm in {:.1} ns: {:.2} nm/us",
        grown * lattice::d110_nm(), sim.t * 1e9, if t_us > 0.0 { grown * lattice::d110_nm() / t_us } else { 0.0 }
    );
    if sim.hit_wall {
        println!("# stopped: growth reached the margin (raise --margin)");
    }
    let (top, _) = sim.face_heights();
    let mut csv = String::new();
    for k in 0..nz {
        let row: Vec<String> = (0..period).map(|n| (top[k * period + n] - top0 as i32).to_string()).collect();
        csv.push_str(&row.join(","));
        csv.push('\n');
    }
    std::fs::write(format!("{out}/top.csv"), csv).unwrap();
    println!("# wrote {out}/log.csv and {out}/top.csv (rows = c, columns = [1-10]; layers above the start)");
}

fn main() {
    let (cmd, m) = parse_args();
    match cmd.as_str() {
        "scan" => scan(&m),
        "grow" => grow(&m),
        "face" => face(&m),
        _ => println!("{HELP}"),
    }
}
