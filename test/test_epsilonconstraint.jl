# ==============================================================================
# fname : test/test_epsilonconstraint.jl
# ==============================================================================
#
# OPTIONAL cross-validation test: checks that the generic epsilon-constraint
# method (JuMP + Gurobi, see src/voptMOA.jl) finds the same number of
# nondominated points as the dedicated three-phase solver, on a couple of
# very small instances only — the generic method is orders of magnitude
# slower (e.g. ~327s on F50-56.txt, vs ~0.02s for the dedicated solver, see
# Table 4 of the paper), so only instances where it completes in well under
# a second are included here.
#
# Requires a working, licensed Gurobi installation. NOT run by
# test/runtests.jl — run it separately and only if you want to exercise this
# path:
#
#   julia --project=. test/test_epsilonconstraint.jl
# ==============================================================================

using Test

const srcDir  = joinpath(@__DIR__, "..", "src")
const dataDir = joinpath(@__DIR__, "..", "data")

const verboseProd = false

include(joinpath(srcDir, "datastru.jl"))
include(joinpath(srcDir, "parser.jl"))
include(joinpath(srcDir, "voptMOA.jl"))  # using JuMP, Gurobi

# (dataset, filename, expected #YN) — all sub-second on the generic method
const referenceEpsilonCst = [
    ("dataFernandez", "F56-57.txt", 3),
    ("dataFernandez", "F55-57.txt", 4),
    ("dataFernandez", "F55-56.txt", 5),
]

@testset "biUFLP2026 — generic epsilon-constraint cross-check (requires Gurobi)" begin
    for (dname, fname, expectedYN) in referenceEpsilonCst
        data = load2UFLP(joinpath(dataDir, dname), fname)
        YN, elapsedTime = epsilonConstraint(data, "GUROBI")
        @testset "$fname" begin
            @test elapsedTime >= 0.0   # epsilonConstraint returns -1.0 on timeout/infeasibility
            @test length(YN) == expectedYN
        end
    end
end
