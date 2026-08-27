# ==============================================================================
# fname : test/runtests.jl
# ==============================================================================
#
# Automated regression tests for the dedicated three-phase 2-UFLP solver
# (biUFLPsolver). Two kinds of checks are performed:
#
#   1. Reference-value regression tests: for a set of benchmark instances,
#      check that the number of boxes after paving (#BoxPav), after
#      refinement (#BoxRed), and the number of nondominated points (#YN)
#      match known reference values. These references come from Table 4 of
#      the paper and were cross-validated multiple times, across several
#      code revisions, during this project's manual test campaigns —
#      including one correction versus the paper itself: H10-10000.txt
#      yields #YN = 493100, not 4931 as printed in the paper's table (this
#      was confirmed to be a typo in the paper, not a solver bug).
#
#   2. Structural invariant tests, which do not depend on knowing the
#      "correct" answer for an instance: the returned nondominated set is
#      checked to be mutually nondominated, and every box's remarkable
#      points are checked to be consistently ordered.
#
# This file intentionally avoids `include`-ing voptMOA.jl (and therefore
# does not `using JuMP, Gurobi`): the dedicated solver does not need them,
# and requiring a licensed Gurobi installation just to run these tests would
# make them fail in environments that only care about the dedicated solver.
# A separate, optional cross-validation test against the generic
# epsilon-constraint method is provided in test/test_epsilonconstraint.jl.
#
# Usage:
#   julia --project=. test/runtests.jl              # fast tier only (default)
#   BIUFLP_FULL_TESTS=true julia --project=. test/runtests.jl   # full tier too
#
# Run as a fresh `julia test/runtests.jl` process (not repeated `include`d
# in a persistent REPL session) — the source files declare `const`s that
# Julia refuses to redefine on a second load in the same session.
# ==============================================================================

using Test

# ------------------------------------------------------------------------------
# load the solver's algorithmic core directly — this project is not packaged
# as an installable module (no src/biUFLP2026.jl entry point), so we mirror
# what main.jl does, minus the JuMP/Gurobi-dependent files.
# ------------------------------------------------------------------------------

const srcDir  = joinpath(@__DIR__, "..", "src")
const dataDir = joinpath(@__DIR__, "..", "data")

# flags expected as globals by the algorithmic files (see main.jl); kept
# quiet here so test output only shows @test results.
const verboseProd = false
const verboseDev  = false
const verboseDev2 = false

include(joinpath(srcDir, "datastru.jl"))
include(joinpath(srcDir, "parser.jl"))
include(joinpath(srcDir, "biUFLP.jl"))
include(joinpath(srcDir, "reduce.jl"))
include(joinpath(srcDir, "mopRoutines.jl"))
include(joinpath(srcDir, "pavingBnB.jl"))
include(joinpath(srcDir, "reducePaving.jl"))
include(joinpath(srcDir, "labeling.jl"))


# ==============================================================================
# reference values: (nBoxPaving, nBoxReducing, #YN)
# ==============================================================================

# fast tier: representative subset, expected to run in well under a minute
# on any machine — safe to run on every change.
const referenceFast = [
    ("dataFernandez", "F56-57.txt", (3, 3, 3)),
    ("dataFernandez", "F55-57.txt", (4, 4, 4)),
    ("dataFernandez", "F55-56.txt", (5, 5, 5)),
    ("dataFernandez", "F54-57.txt", (8, 8, 20)),
    ("dataFernandez", "F53-57.txt", (9, 9, 173)),
    ("dataFernandez", "F52-54.txt", (8, 8, 47)),
    ("dataFernandez", "F51-56.txt", (14, 12, 755)),
    ("dataFernandez", "F50-56.txt", (23, 17, 729)),
    ("dataFernandez", "F50-51.txt", (152, 10, 1229)),  # worst-case #BoxPav in dataset F
    ("dataHarris",    "H10-2000.txt", (11, 11, 13412)),
]

# full tier: every instance of datasets F and H, including the two largest H
# instances (H10-8000 / H10-10000), whose labeling phase alone can take
# 10-20 seconds depending on the machine. Opt-in only.
const referenceFull = [
    ("dataFernandez", "F28-29.txt", (10, 6, 244)),
    ("dataFernandez", "F50-52.txt", (11, 5, 408)),
    ("dataFernandez", "F50-53.txt", (21, 7, 771)),
    ("dataFernandez", "F50-54.txt", (16, 13, 700)),
    ("dataFernandez", "F50-55.txt", (13, 9, 513)),
    ("dataFernandez", "F50-57.txt", (17, 14, 616)),
    ("dataFernandez", "F51-52.txt", (38, 15, 635)),
    ("dataFernandez", "F51-53.txt", (32, 9, 1047)),
    ("dataFernandez", "F51-54.txt", (24, 19, 1013)),
    ("dataFernandez", "F51-55.txt", (31, 22, 1111)),
    ("dataFernandez", "F51-57.txt", (18, 12, 796)),
    ("dataFernandez", "F52-53.txt", (14, 8, 435)),
    ("dataFernandez", "F52-55.txt", (3, 3, 20)),
    ("dataFernandez", "F52-56.txt", (6, 6, 15)),
    ("dataFernandez", "F52-57.txt", (6, 6, 16)),
    ("dataFernandez", "F53-54.txt", (7, 7, 333)),
    ("dataFernandez", "F53-55.txt", (6, 6, 306)),
    ("dataFernandez", "F53-56.txt", (5, 4, 318)),
    ("dataFernandez", "F54-55.txt", (5, 5, 37)),
    ("dataFernandez", "F54-56.txt", (4, 4, 22)),
    ("dataHarris",    "H10-4000.txt",  (11, 11, 69167)),
    ("dataHarris",    "H10-6000.txt",  (13, 13, 172278)),
    ("dataHarris",    "H10-8000.txt",  (12, 12, 731385)),
    ("dataHarris",    "H10-10000.txt", (12, 12, 493100)),  # see header note
]


# ==============================================================================
"""
    isMutuallyNondominated(points::Vector{Tuple{Int64,Int64}})

    Structural check: true if no point in the list weakly dominates another
    distinct point. Used to validate the output of biUFLPsolver independently
    of any instance-specific reference value.

    Implemented in O(n log n) (sort by f1, then a single pass checking that f2
    is strictly decreasing) rather than a naive O(n²) pairwise comparison —
    the latter would be far too slow on the largest H instances, where #YN
    reaches several hundred thousand points.
"""
function isMutuallyNondominated(points::Vector{Tuple{Int64,Int64}})
    length(points) <= 1 && return true
    sorted = sort(points, by = p -> (p[1], p[2]))
    for i in 2:length(sorted)
        sorted[i] == sorted[i-1] && continue  # identical duplicate: not a violation
        if sorted[i][2] >= sorted[i-1][2]
            return false
        end
    end
    return true
end


# ==============================================================================
"""
    runInstanceTests(dname, fname, expected)

    Solve one instance with the dedicated solver and check both the
    reference values and the structural invariants.
"""
function runInstanceTests(dname::String, fname::String, expected::Tuple{Int64,Int64,Int64})
    data = load2UFLP(joinpath(dataDir, dname), fname)
    paving, tPav, tRed, nBoxPav, nBoxRed, tLab, ND_YN = biUFLPsolver(data)

    @testset "$fname" begin
        @test nBoxPav        == expected[1]
        @test nBoxRed        == expected[2]
        @test length(ND_YN)  == expected[3]

        @test isMutuallyNondominated(ND_YN)

        # every remaining box must have consistently ordered remarkable points:
        # yI <= y12 <= yN and yI <= y21 <= yN, componentwise
        for b in paving
            @test all(b.yI .<= b.y12) && all(b.y12 .<= b.yN)
            @test all(b.yI .<= b.y21) && all(b.y21 .<= b.yN)
        end
    end
end


# ==============================================================================
@testset "biUFLP2026 — dedicated three-phase solver (fast tier)" begin
    for (dname, fname, expected) in referenceFast
        runInstanceTests(dname, fname, expected)
    end
end

if get(ENV, "BIUFLP_FULL_TESTS", "false") == "true"
    @testset "biUFLP2026 — dedicated three-phase solver (full tier)" begin
        for (dname, fname, expected) in referenceFull
            runInstanceTests(dname, fname, expected)
        end
    end
else
    println("Full test tier skipped (set BIUFLP_FULL_TESTS=true to include it).")
end
