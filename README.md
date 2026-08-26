# 2-UFLP

Solver for the bi-objective 0/1 Uncapacitated Facility Location Problem (**2-UFLP**).

This code implements the three-phase algorithm described in:

> Xavier Gandibleux, Anthony Przybylski — *Paving and computing the set of
> nondominated points for the bi-objective 0/1 uncapacitated facility
> location problem*.
> Preprint available on Optimization Online:
> <https://optimization-online.org/?p=34431>

## Algorithm overview

The solver computes the complete set of nondominated points $Y_N$ in three phases:

1. **Paving** — a branch-and-bound builds a set of "boxes" that cover $Y_N$,
   pruning subsets of candidate facilities with three dominance-based tests.
2. **Refinement** — boxes are shrunk or eliminated pairwise using supported
   points computed within each box.
3. **Generation** — a labeling algorithm (dynamic programming) solves the
   allocation subproblem associated with each remaining box and extracts the
   locally nondominated points.

A generic ε-constraint method (`JuMP` + `Gurobi`) is available for cross-validation.

## Requirements

- **Julia** ≥ 1.10 (developed and tested with Julia 1.12.5)
- **Gurobi** with a valid license (academic or commercial).
- Julia packages to install beforehand (`using Pkg; Pkg.add("...")` in the REPL, or via the package manager):
  - `JuMP`, `Gurobi`
  - `DataStructures` 
  - `PyPlot`  
  - `Printf`, `Test` — part of the standard library, nothing to install

## Project structure

```
biUFLP2026/
├── data/
│   ├── dataDidactic/     # didactic instance (8 customers, 5 facilities)
│   ├── dataFernandez/    # 28 instances F, 30 customers × 90 facilities
│   └── dataHarris/       # 5 instances H, 10 customers × {2000..10000} facilities
├── output/               # generated figures (created automatically)
├── test/
│   ├── runtests.jl               # automated regression tests (specific solver)
│   └── test_epsilonconstraint.jl # optional cross-check (requires Gurobi)
└── src/
    ├── main.jl           # main entry point
    ├── datastru.jl       # data structures (Instance, Box, Label)
    ├── parser.jl         # instance file reading
    ├── biUFLP.jl         # orchestration of the 3 phases (biUFLPsolver)
    ├── pavingBnB.jl      # phase 1 — branch-and-bound paving
    ├── reduce.jl         # phase 2 — local reduction routines
    ├── reducePaving.jl   # phase 2 — paving refinement
    ├── labeling.jl       # phase 3 — labeling algorithm    
    ├── mopRoutines.jl    # common multi-objective routines
    ├── voptMOA.jl        # generic ε-constraint method (JuMP + Gurobi)
    └── graphicpyplot.jl  # optional visualization (PyPlot backend)
```

## Usage

```bash
cd src
julia main.jl
```

The working directory is automatically fixed to the location of the file
(`cd(@__DIR__)`), so this command behaves identically whether run from a
REPL, VS Code, or the command line.

### Configuration

Everything is set at the top of `src/main.jl`:

| Constant | Role |
|---|---|
| `experiment` | `false`: a single instance; `true`: every instance of a collection |
| `epsilonCst` | `false`: specific solver (`biUFLPsolver`); `true`: generic solver (`ε-constraint`) |
| `graphics` | enables figure generation (PyPlot backend) — single-instance mode only |
| `verboseProd` / `verboseDev` / `verboseDev2` | on-screen trace levels |
| `dnameRun` / `fnameRun` | folder / file of the instance to solve (single-instance mode) |
| `dnameExperiment` | collection to iterate over in experiment mode |
| `figDir` | output folder for figures (`../output` by default) |

To change instance or collection, comment/uncomment the corresponding lines
in the `# instance selection` block.

### Output

**Specific solver** (`epsilonCst = false`):

```
      fnames         tPav         tRed         tTot    #BoxPav    #BoxRed         tLab        #YN         tTOT
      F50-56     0.004700     0.017400     0.022100         23         17     0.003700        729     0.025800
```

`tPav`/`tRed`/`tLab`: time (s) of phases 1/2/3; `#BoxPav`/`#BoxRed`: number
of boxes after paving/refinement; `#YN`: cardinality of the resulting
nondominated set.

**Generic solver** (`epsilonCst = true`):

```
      fnames         tOpt        #YN
      F50-56   327.079400        729
```

### Figures

If `graphics = true`, a figure `<instance>.png` is saved to `output/` at
the end of the resolution (single-instance mode only).

## Datasets

- **`dataDidactic`** — toy instance (8 customers, 5 facilities) used in the
  paper as a step-by-step illustration of the algorithm.
- **`dataFernandez`** — 28 instances derived from Fernández & Puerto (2003),
  30 customers × 90 facilities, objective correlation ≈ 0.
- **`dataHarris`** — 5 instances from Harris et al. (2009, 2011),
  10 customers × {2000, 4000, 6000, 8000, 10000} facilities, correlation ≈
  0.99 (green logistics context: cost and CO₂ emissions).

Full details (cost ranges, provenance): see the paper.

## Tests

Automated regression tests are provided in `test/`, checking that the
solver returns the correct number of nondominated points and the correct
box counts on a set of reference instances, plus a few structural
invariants (nondomination of the returned set, consistency of each box's
remarkable points) that do not depend on any hardcoded reference value.

Run the default (fast) test tier:

```bash
julia test/runtests.jl
```

Run the full tier as well (every instance of datasets F and H, including
the two largest H instances, whose labeling phase alone can take 10-20
seconds):

```bash
BIUFLP_FULL_TESTS=true julia test/runtests.jl
```

An optional, separate cross-check against the generic ε-constraint method
is available in `test/test_epsilonconstraint.jl` (requires a licensed
Gurobi installation; not run by `runtests.jl`):

```bash
julia test/test_epsilonconstraint.jl
```

## License

MIT.
