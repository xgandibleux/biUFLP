# ==============================================================================
# fname : main_didactic.jl
# ==============================================================================
#
# Minimal, standalone entry point restricted to the didactic instance
# (didactic1.txt / didactic2.txt), intended for presentation and teaching
# purposes. It solves the instance with the three-phase algorithm and draws
# the resulting paving with the Luxor backend (graphicluxor.jl).
#
# The Luxor backend is only meant to work on this small instance: it uses a
# fixed canvas size and hardcoded box indices for some illustrative overlays
# (see graphicluxor.jl), which do not generalize to the benchmark datasets.
#
# For running the benchmark datasets (F, H), the numerical experiments, or
# the PyPlot graphics backend, use main.jl instead.
# ==============================================================================

println("Run bi-objective 0/1 UFLP solver (didactic instance, presentation mode)")
using Printf

# ==============================================================================
# working directory: fix to the location of this file so that relative paths
# (e.g. "../data/...") resolve identically whether run from a REPL, VS Code,
# or the command line — regardless of the shell's current directory at launch.
cd(@__DIR__)

# ==============================================================================
# parameters to control the app

const verboseProd  = true   # displays on the terminal the useful info
const verboseDev   = false  # displays on the terminal the development info
const verboseDev2  = false  # displays on the terminal the development info for phase 2 (labeling)

# ==============================================================================
# instance selection — fixed to the didactic dataset
# uncomment one line; leave only one active

#const dnameRun = "../data/dataDidactic" ; const fnameRun = "didactic1.txt"
const dnameRun = "../data/dataDidactic" ; const fnameRun = "didactic2.txt"

# ==============================================================================
# files of the app

include("datastru.jl")
include("parser.jl")
include("biUFLP.jl")
include("reduce.jl")
include("mopRoutines.jl")
include("pavingBnB.jl")
include("reducePaving.jl")
include("labelingBourrinOrdrevct.jl")
include("graphicluxor.jl")


# ==============================================================================
function main()

    verboseProd ? println("\n>>> [0] LOADING DATASET ") : nothing
    data = load2UFLP(dnameRun, fnameRun)

    paving,timePaving,timeReducing,nbBoxPaving,nbBoxReducing,timeLabeling,ND_YN = biUFLPsolver(data)

    println("run completed! \n")
    println("      fnames         tPav         tRed         tTot    #BoxPav    #BoxRed         tLab        #YN         tTOT")
    @printf(" %11s   %10.6f   %10.6f   %10.6f     %6d     %6d   %10.6f     %6d   %10.6f\n",
            fnameRun[1:end-4], timePaving, timeReducing, timePaving+timeReducing,
            nbBoxPaving, nbBoxReducing, timeLabeling, length(ND_YN),
            timePaving+timeReducing+timeLabeling)

    # draw the paving with the Luxor backend (see graphicluxor.jl for what is drawn)
    YN = []
    drawGraphics(paving, YN)

    return nothing
end

# ==============================================================================
main()

finish()
preview()
