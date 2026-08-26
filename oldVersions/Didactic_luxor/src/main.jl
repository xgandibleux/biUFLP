# ==============================================================================
# fname : main.jl
# March 2026
# ==============================================================================

println("Run bi-objective 0/1 UFLP solver")
using Printf

# ==============================================================================
# working directory: fix to the location of this file so that relative paths
# (e.g. "../data/...") resolve identically whether run from a REPL, VS Code,
# or the command line — regardless of the shell's current directory at launch.
cd(@__DIR__)

# ==============================================================================
# parameters to control the app

const verboseProd  = false  # displays on the terminal the useful info in production
const verboseDev   = false  # displays on the terminal the development info
const verboseDev2  = false  # displays on the terminal the development info for phase 2 (labeling)
const experiment   = false  # run a full numerical experiment

const vOptSolver   = false   # run vOptGeneric with ϵ-constraint and glpk or gurobi
const graphics     = false   # display on the screen graphics (PyPlot backend)

# ==============================================================================
# instance selection

# --- used when experiment = false (single instance run) ---------------------
# uncomment one (dname, fname) pair; leave only one pair active
#const dnameRun = "../data/dataDidactic"  ; const fnameRun = "didactic1.txt"
#const dnameRun = "../data/dataDidactic"  ; const fnameRun = "didactic2.txt"
const dnameRun  = "../data/dataFernandez" ; const fnameRun = "F50-56.txt"
#const dnameRun = "../data/dataHarris"    ; const fnameRun = "H10-2000.txt"
#const dnameRun = "../data/dataHarris"    ; const fnameRun = "H10-10000.txt"

# --- used when experiment = true (run all instances of one collection) ------
# uncomment one line; leave only one active
const dnameExperiment = "../data/dataHarris"
#const dnameExperiment = "../data/dataFernandez"

# ==============================================================================
# graphics output
# used only when graphics=true (single-instance run only; graphics is not
# meant to be combined with experiment=true, which would otherwise produce
# one figure per instance in the collection)
const figDir = "../output"  # directory where figures are saved; created if missing

# ==============================================================================
# files of the app

include("datastru.jl")
include("parser.jl")
#include("vopt.jl")
include("voptMOA.jl")
include("biUFLP.jl")
include("reduce.jl")
include("mopRoutines.jl")
include("pavingBnB.jl")
include("reducePaving.jl")
include("labelingBourrinOrdrevct.jl")


if graphics
    include("graphicpyplot.jl")
end


# ==============================================================================
function main()

    # prepare the run for a single resolution or a numerical experiment --------
    if experiment
        # prepare the experiment on a list of instances ------------------------

        dname = dnameExperiment

        # get all the names of the datasets for a given collection available in folder 'dnameE'
        fnames = getfnames(dname)

        # matrix storing the results
        results = Array{Any, 2}(undef, length(fnames), 7)

        # run a didactic instance to compile all the code
        data = load2UFLP("../data/dataDidactic","didactic1.txt")
        if !vOptSolver
            biUFLPsolver(data)
        else
            #vOpt_MOA(data,"GUROBI")
            epsilonConstraint(data,"GUROBI") 
        end


    else
        # setup for a single instance -----------------------------------------

        dname   = dnameRun
        fnames  = [fnameRun]

        # matrix storing the results
        results = Array{Any, 2}(undef, 1, 7)
    end

    # call the solver ---------------------------------------------------------
    for i in eachindex(fnames)
        verboseProd ? println("\n>>> [0] LOADING DATASET ") : nothing

        data  = load2UFLP(dname,fnames[i])
        if vOptSolver
            #getTime = time()
            #YN, timevOPt = vOpt_MOA(data,"GUROBI") 
            YN, timevOPt = epsilonConstraint(data,"GUROBI") 
            #timevOPt = round(time()- getTime, digits=4)
            results[i,1] = fnames[i]
            results[i,6] = timevOPt
            results[i,7] = length(YN)
            if graphics
                thecolor = "green"; themarkerstyle = "o"; themarkersize = 15
                displayYN(YN,thecolor,themarkerstyle,themarkersize)
            end
            #@show YN
        else
            paving,timePaving,timeReducing,nbBoxPaving,nbBoxReducing,timeLabeling, ND_YN=biUFLPsolver(data) # TODO: also return YN and XE once integrated into vOptSolver
            results[i,1] = fnames[i]
            results[i,2] = timePaving
            results[i,3] = timeReducing
            results[i,4] = nbBoxPaving
            results[i,5] = nbBoxReducing
            results[i,6] = timeLabeling
            results[i,7] = length(ND_YN)
            #@show ND_YN

            if graphics
                setupGraphic(fnames[1], data, paving)
                displayPaving(paving)
                #displayAllRunningCosts(paving)
                thecolor = "red"; themarkerstyle = "o"; themarkersize = 15
                displayYN(ND_YN,thecolor,themarkerstyle,themarkersize)

                # save the figure to file (PyPlot does not display
                # interactively in every environment, e.g. some VS Code
                # setups, but savefig works regardless)
                mkpath(figDir)
                savefig(joinpath(figDir, fnames[1][1:end-4] * ".png"))
            end

        end
    end

    println("run completed! \n")

    if !vOptSolver
        println("      fnames         tPav         tRed         tTot    #BoxPav    #BoxRed         tLab        #YN         tTOT")
        for i in 1:size(results,1)
            @printf(" %11s   %10.6f   %10.6f   %10.6f     %6d     %6d   %10.6f     %6d   %10.6f\n", results[i,1][1:end-4], results[i,2], results[i,3], results[i,2]+results[i,3], results[i,4], results[i,5], results[i,6], results[i,7],results[i,2]+results[i,3]+results[i,6])
        end

    else
        println("      fnames         tOpt        #YN")
        for i in 1:size(results,1)
            @printf(" %11s   %10.6f     %6d\n", results[i,1][1:end-4], results[i,6], results[i,7])
        end
    end

    return nothing
end

# ==============================================================================
main()
