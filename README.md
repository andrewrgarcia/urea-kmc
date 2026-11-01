# urea-kmc

Three-dimensional kinetic Monte Carlo of urea crystals growing from water and methanol, built to
replicate the morphology results of

> S. Piana, M. Reyhani & J. D. Gale, "Simulating micrometre-scale crystal growth from solution",
> *Nature* **438**, 70–73 (2005), doi:10.1038/nature04173, and its Supplementary Information.

The site rates are the paper's own, measured by counting transitions in molecular dynamics
(`data/piana_gale_2005_rates.json`, decoded from the SI). Everything else (the lattice, the
neighbor definitions, the kMC algorithm) is written from the paper and the crystal structure, and every
choice the paper leaves open is a command-line switch, listed below.

## Run

```
makers test                      # geometry, pinned SI values, bookkeeping
makers scan                      # net growth vs concentration (seconds)
makers scan-big                  # larger nucleus, finer grid (a few minutes)
makers grow                      # the water needle run: C = 0.24, reading --ab co
CS=0.30 makers grow              # growth at another concentration (C = CS * (1 + SIGMA), SIGMA default 0)
AB=nh2 makers scan-big           # any task with Fig. 2 as drawn (AB=co is the default)
SOLVENT=methanol CS=0.21 makers face && SOLVENT=methanol makers plot-face
                                 # one [110] side face, 200 x 200 nm, periodic: does it grow?
SOLVENT=methanol CS=0.19 ISLAND_NM=20 FACE_TNS=100 makers face
                                 # one seeded island per face: how fast do its edges move toward +c and -c?
makers plot                      # runs/water/morphology.png
SOLVENT=methanol makers scan     # any task, other solvent
```

or directly: `cargo run --release -- help`.

## The model

**Lattice.** One site per molecule on the real urea structure (P-42₁m, a = 5.565 Å, c = 4.684 Å,
two molecules per cell). The two molecules of the cell point their C=O in opposite directions along c and
lie in perpendicular planes; that is the paper's "dipole up or down." Every molecule has six hydrogen-bonded
neighbors (the paper's Fig. 2):

| label | where | count |
|---|---|---|
| A | next molecule along c, on this molecule's NH₂ side | 1 |
| B | next molecule along c, on this molecule's C=O side | 1 |
| 1 | in-plane, along this molecule's own plane: its N–H bond to their O | 2 |
| 2 | in-plane, across: their N–H bond to this molecule's O | 2 |

A site's type is the set of occupied neighbors, e.g. `12/B` = one type-1, one type-2 and B. That gives the
paper's 34 surface types (10 kinks, 8 steps, 4 terraces, 12 with one or two neighbors).

**Reading the SI tables.** Rows are the c-axis neighbors (`-`, `A`, `B`, `AB`), columns the in-plane ones
(`-`, `1`, `11`, `2`, `12`, `112`, `22`, `122`, `1122`). Cells with three numbers are the second-neighbor
split: (i) more second neighbors "along Z", (ii) equal, (iii) more "along X or Y."

**Dynamics.** Exact continuous-time kMC (the Bortz–Kalos–Lebowitz "n-fold way"): a crystal molecule
dissolves at `kd(type)`; an empty site touching the crystal fills at `kc(type) × C`. C is the urea
concentration as a fraction of lattice sites (1.0 = crystal density = 22.89 M), held constant (an infinite
reservoir). The paper instead swept all sites every 50 ps with probability `k × 50 ps`; that scheme is
approximate, and at these rates some probabilities exceed 1 (`kc × C × Δt` up to 2.2). This code has no
time step.

## Choices the paper leaves open (switches)

| ⚠ | what | default | switch |
|---|---|---|---|
| 1 vs 2 | which in-plane pair is "1": read from the projection in Fig. 2 (shape barely depends on it) | 1 = along the molecule's own plane | `--swap12` |
| A vs B | which chain neighbor is "A": Fig. 2 shows A on the NH₂ side, but the SI rates match the paper's solubilities only the other way round | A = C=O side (since 2026-10-26) | `--ab nh2\|co\|crystal` (`--swapAB` = `--ab co`) |
| split | what "second neighbors along Z / along X or Y" are: here the four non-bonded images 0.652 c away (4.98 Å) vs the four same-orientation molecules at ±a, ±b (5.565 Å) | use i / ii / iii | `--split mean` ignores the split |
| unobserved | cells never seen in MD (0 in both tables, e.g. `11/-`, `22/AB`): a molecule there could never leave | fill from the mean of observed cells with the same neighbor count; bulk never dissolves | `--unobserved published` |
| units | C as a fraction of lattice sites; molarity = C × 22.89 M | | |
| placeholders | some cells are 1.000E+10 in both tables (e.g. `112/-`, `122/-`), probably filled in by hand | used as published | |
| solution | well-mixed, constant C, no diffusion, no nucleation away from the crystal | | (finite-solution depletion: not yet) |

## Results so far (2026-10-26, water and methanol, defaults)

**Concentration scans** (`makers scan`, `makers scan-big`): net growth over 100 ns from a cube seed,
and where it changes sign.

| seed | water | methanol | paper's kMC | experiment |
|---|---|---|---|---|
| 16 cells (~9 nm) | 0.47 of sites (10.8 M) | | | |
| 40 cells (~22 nm) | 0.43 (9.9 M) | 0.37 (8.5 M) | 5 M / 4 M | 11 M / 4 M |

The order is right (urea less soluble in methanol) and the bigger seed gives a lower value, as it should.
Both numbers are almost exactly twice the paper's (9.9/5 = 1.98, 8.5/4 = 2.1); a factor-of-two
difference in how concentration is counted would explain both, but nothing yet shows which.
Same seed, same output: the methanol scan run twice was identical.

**But the scan overstates growth.** A 1,000-ns run at C = 0.444 (water, 3% above the 40-cell crossover,
40-cell seed, 100x100x500 box) grows for about 150 ns while the cube fills its corners into a square prism
with {110} sides, then dissolves steadily: 128,000 → 208,000 molecules at 100 ns → 4,000 at 1,000 ns.
Length/width stays near 0.9: no needle. So the 100-ns scan was mostly measuring that corner-filling
burst, and the break-even for a squared-off crystal is above 0.444.

**Which reading of Fig. 2? (2026-10-26)** Four readings, water, 40-cell seed, 140x140x600 box, 60 ns;
speeds in nm per microsecond measured after the first 20 ns (each c end vs each {110} side):

| reading | C | end speed | side speed | end/side |
|---|---|---|---|---|
| default | 0.46 / 0.50 / 0.55 | 35 / 111 / 187 | 7 / 49 / 101 | 4.8 / 2.3 / 1.9 |
| `--swap12` | same | 47 / 100 / 187 | 20 / 57 / 98 | 2.4 / 1.8 / 1.9 |
| `--swapAB` | same | 945 / 1064 / 1250 | 175 / 220 / 285 | 5.4 / 4.8 / 4.4 |
| both | same | 939 / 1045 / 1224 | 175 / 223 / 283 | 5.4 / 4.7 / 4.3 |

- Which in-plane pair is "1" barely matters: `--swap12` changes little.
- Which chain neighbor is "A" matters a lot. With `--swapAB` (A on the molecule's C=O side, B on its
  NH₂ side) the 40-cell scans break even at **0.20 of sites = 4.7 M in water** and **0.15 = 3.4 M in
  methanol**, against the paper's 5 M and 4 M. The factor of two in the default reading is gone without
  touching the concentration units. So the SI's rows A and B most likely mean the opposite of the
  NH₂-up / C=O-down picture in Fig. 2.
- In every reading the ends outgrow the sides, and more so closer to break-even; no needle yet at these
  supersaturations. The paper's needle is at 3% above saturation.

**The water needle (2026-10-26, `--swapAB`, 40-cell seed, 100x100x500 box, Andrew's machine).**

| C | what happened |
|---|---|
| 0.21 (4.8 M) | ends grow while sides dissolve: 18 → 81 nm long, 31 → 12 nm wide, then the thin needle dissolves (250 ns) |
| 0.24 (5.5 M) | sides frozen at 29 nm from 25 ns on; each end grows ~0.27 nm/ns; 18 → 231 nm long when it reached the box end (395 ns, 21 s of compute) |

At 0.24 the sides did not advance by one molecular layer while the ends grew 200 nm: end speed over side
speed above ~170. The plot shows a square rod with flat {110} sides (five-neighbor molecules) and rough
ends, the paper's Fig. 4a ("large [110] faces", no growth on them at this supersaturation). The aspect
ratio of 7.7 is a box limit, not the crystal's.

**Methanol, and a third reading (2026-10-26).** End shape measured as the width 1–30 nm from each end,
along [110] and [1-10] separately.

| reading | break-even water / methanol | water | methanol |
|---|---|---|---|
| `--ab nh2` (Fig. 2 as drawn) | 9.9 M / 8.5 M | no needle; stubby | stubby (length ≈ 1.5 × width), chisel ends twisted 90°: top 5 nm thin along [110] and 33 nm along [1-10], bottom the reverse (C = 0.42) |
| `--ab co` (default since 2026-10-26) | 4.7 M / 3.4 M | needle, sides frozen (C = 0.24) | rod with pointed ends, only a weak 90° twist (top 10.6 vs 7.9 nm, bottom 9.1 vs 11.8 nm, 3 nm from the end); sides grow slowly, ends ~90× faster (C = 0.22) |
| `--ab crystal` (A = the molecule above, for every molecule) | 6.4 M / 4.9 M | crystal walks: grows at the top, dissolves at the bottom | same |
| paper | 5 M / 4 M | needle, aspect > 1000 | aspect 20–100× smaller than water, [001] replaced by [111], ends twisted 90° |

- `--ab crystal` is ruled out on physics alone: the top and bottom (001) faces of urea are related by the
  crystal's -4 symmetry, so they cannot grow and dissolve in the same solution. A walking crystal means
  the rate table, read that way, breaks the symmetry.
- `--ab co` reproduces the numbers the paper reports (both solubilities, the water needle) but only
  partly the methanol shape; `--ab nh2` reproduces the methanol ends but neither number. Choosing a reading
  because it matches the target is one fitted bit; any write-up says so and shows both.

**Default reading and methanol near break-even (2026-10-26).** `--ab co` is now the default (Andrew's decision).
Methanol, 60-cell seed, 140x140x500 box:

| C | shape |
|---|---|
| 0.19 | sides shrink 46 → 37 nm, ends grow 28 → 232 nm (box end at 248 ns) |
| 0.20 | sides 46 → 43 nm, ends to 231 nm (208 ns) |
| 0.21 | sides frozen at 46 nm, ends to 232 nm (187 ns) |

Ends are pointed, tapering to ~5–20 nm over the last 15 nm, and the taper is twisted 90° between the two
ends in every run (8 nm from the end at C = 0.19: top 11 nm wide along [110] and 7 along [1-10], bottom 6
and 10). So in this reading methanol gives the paper's pointed, twisted ends, but not its stubby shape:
the crystal still grows mostly along c. In water the sides are frozen and the ends nearly flat; in
methanol the sides creep (end/side ≈ 90 at 0.22) and the ends taper. Whether methanol ends up 20–100×
less elongated than water needs a true aspect ratio for each, which these boxes are too short to give.

**The paper's own time stepping (2026-10-26).** `--scheme sync` visits every surface site once per 50 ps
and flips it with probability rate × 50 ps, capped at 1 (`--order random|lattice` for the visiting order).
Random order, 40-cell seed, 200 ns: water at 0.24 is a needle with sides frozen at 28 nm (length 18 → 150 nm);
methanol at 0.21 is also a needle (18 → 223 nm, sides 28 nm). So the cap alone does not make methanol stubby.

Lattice order (Andrew's machine, 200 ns): same shapes as random order, and no end runs ahead (water: top end
68 nm, bottom 67 nm; methanol 95 / 97 nm). The time-stepping method does not change the morphology here,
which agrees with the paper's SI ("the two approaches give very similar results").

**Size, the likely reason methanol is still a needle (2026-10-26).** The paper's simulated methanol crystal
(Fig. 4b) is over 1 µm wide and its [110] side face is covered with growth islands ("birth and spread").
Ours are 30–60 nm wide. A flat face grows only when an island forms on it, and islands form in proportion to
the face's area, so small faces stay frozen; the paper itself calls the methanol aspect ratio
"size-dependent". In our runs the methanol plot is spindle-shaped with rough sides and pointed, twisted
ends; the water plot is a square rod with smooth faces and flat ends. Doubling the seed (80 cells, 62 nm
wide, C = 0.21, 100 ns) still leaves the sides frozen. A full 3D crystal 1 µm wide is out of reach on a dense
lattice; a periodic [110] slab (the paper's own 2D test, 0.6 × 0.6 µm) is the way to measure side growth at
that scale.

**The flat-face test (2026-10-26).** `urea-kmc face` grows one [110] side face, periodic along [1-10] and c
so it has no edges, starting from whole molecular layers (each layer 0.39 nm). It reports layers gained per
face and writes the top face's height map (`makers plot-face`). First check, 50 x 50 nm, 50 ns: water at
C = 0.24 gains 0.007 layers (frozen: single molecules come and go); methanol at 0.21 gains 0.49 layers and
speeds up (0.05 → 0.49 over five 10-ns steps), and the height map shows second- and third-layer islands
spreading, stretched along c. That is the paper's claim seen directly: water's side face does not grow,
methanol's grows by islands.

**Side-face speeds at 200 x 200 nm and the predicted shapes (2026-10-26, Andrew's machine, 200 ns each).**

| run | layers per face in 200 ns | side speed | height map |
|---|---|---|---|
| water, C = 0.24 | 0.076, still speeding up (0.042 in the last 50 ns) | 0.15 nm/µs average, 0.33 in the last 50 ns | two or three late islands, long streaks along c |
| methanol, C = 0.21 | 3.13, steady after ~80 ns | 6.8 nm/µs (100-200 ns) | many islands up to six layers high, stretched along c ("birth and spread", like the paper's Fig. 1c) |
| methanol, C = 0.19 | 0.021, flat | ~0 | no islands, single molecules only |

Predicted length/width at full size = end speed (3D runs: water 0.24 → 266 nm/µs per end; methanol 0.21 →
532 nm/µs per end) / side speed: **water ~800–1,800** (paper: > 1,000), **methanol ~80** (paper predicts
10–40). Water over methanol ≈ 10–20× (paper: 20–100×). Side growth on methanol switches on between C = 0.19
and 0.21, so these numbers move strongly with concentration; the two solvents here are not at matched
supersaturation (water 0.24 is ~18% above its scan crossover, methanol 0.21 ~40%).

**Are the measured rates bond-additive? (2026-10-26)** Least-squares fit of ln(kd/kc), ln kd and ln kc over
all measured cells (zeros and the 1.000E+10 placeholders excluded) to four per-neighbor terms plus a constant:

| | ln(kd/kc) | ln kd | ln kc | per-neighbor ΔG for the ratio, kT (1, 2, A, B) |
|---|---|---|---|---|
| water (44 values) | R² 0.93, scatter ×2.5 | R² 0.83 | R² 0.83 | −2.5, −2.6, −1.8, −4.7 |
| methanol (48) | R² 0.85, scatter ×3.2 | R² 0.73 | **R² 0.49** | −2.3, −2.0, −1.1, −4.1 |

Four bond energies reproduce the equilibrium ratio of every site type to within the MD counting noise; the two
in-plane bonds are equal, the two chain bonds are not (that asymmetry is the polar growth). The attachment
rates alone in methanol are not additive: that is the kinetic information a bond-counting model (CrystalGrower's
energy ladder) cannot carry, and what measuring the rates buys. Script: the fit is twelve lines of numpy over
`data/piana_gale_2005_rates.json` (tools/additivity.py).

**Polar step speeds: the paper's "twice as fast" (2026-10-26, cloud).** `--island-nm 20` puts one 20 x 20 nm
island, one layer high, at the center of each face of a 200 nm slab and tracks how far its edges move (median
over the island's central rows and columns; island.csv). Speeds from 10 to 100 ns, nm/µs:

| run | face | toward +c | toward −c | +c / −c | along [1-10] (each side) |
|---|---|---|---|---|---|
| methanol 0.19 | top (+[110]) | 343 | 36 | 9.4 | ~0, slightly shrinking |
| | bottom (−[110]) | 354 | 94 | 3.8 | ~0 |
| methanol 0.20 | top | 437 | 198 | 2.2 | 26–61 |
| | bottom | 494 | 286 | 1.7 | 61–70 |
| water 0.24 | top | 567 | 338 | 1.7 | 61–79 |
| | bottom | 510 | 396 | 1.3 | 79 |

- In every run, on both faces, steps move faster toward +c than toward −c. The paper: "on the [110] and
  [-1-10] faces, deposition on steps in the +c direction is twice that in the −c one", reversed on [1-10] and
  [-110] (which our model gets automatically: the -4 symmetry that maps those faces also flips c).
- At C = 0.20 the ratio is ≈ 2 (2.2 and 1.7): the paper's number. Closer to break-even (0.19) the −c step almost
  stalls and the ratio grows (4–9). Water is less polar (1.3–1.7); the paper claims the asymmetry for methanol.
- The 0.19 height map shows the paper's other sentence directly: the island grows into a strip with its
  starting square near one end ("the initial nucleation site remains close to one edge of the island").
- Steps run 6–7× faster along c than along [1-10], which is why every island in these maps is a streak along c.
- Caveats: one run per condition; at 0.20 islands also form on their own and some may touch the seeded one;
  which end the paper calls +c is a sign convention we cannot check, only the factor.

## What differs from the paper on purpose

- exact continuous-time kMC by default instead of 50-ps sweeps (`--scheme sync` runs theirs; the comparison is itself worth a figure);
- explicit, switchable choices instead of unstated ones;
- no solid-on-solid restriction: overhangs and vacancies are allowed.

## Layout

- `src/lattice.rs` geometry, derivation in the comments, tests on distances and symmetry
- `src/rates.rs`, `src/rates_table.rs` (generated) the SI tables and the switches
- `src/sim.rs` the kMC engine
- `src/main.rs` `scan` and `grow`
- `tools/gen_rates.py` JSON → Rust table; `tools/plot.py` projections shaded like Fig. 4; `tools/plot_face.py` face height map

## License

The code is MIT-licensed (`LICENSE`). The rate tables in `data/piana_gale_2005_rates.json` (and
`src/rates_table.rs`, generated from them) are a transcription of the Supplementary Information of
S. Piana, M. Reyhani and J. D. Gale, *Nature* **438**, 70–73 (2005), doi:10.1038/nature04173. They are
not covered by the MIT license; cite the original paper when you use them.
