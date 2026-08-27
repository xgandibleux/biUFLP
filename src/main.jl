# ==============================================================================
# bi-objective 0/1 UFLP solver  
#   fname         : main.jl
#   author        : Xavier Gandibleux
#   last revision : August 2026
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

const experiment   = false  # run only one instance or a full numerical experiment
const epsilonCst   = false  # run JuMP with ϵ-constraint and Gurobi
const graphics     = false  # plot graphics (PyPlot backend)

const verboseProd  = false  # display on the terminal the useful info in production
const verboseDev   = false  # display on the terminal the development info
const verboseDev2  = false  # display on the terminal the development info for phase 2 (labeling)

# ==============================================================================
# instance selection

# --- used when experiment = false (single instance run) ---------------------
# uncomment one (dname, fname) pair; leave only one pair active
#const dnameRun = "../data/dataDidactic"  ; const fnameRun = "didactic1.txt"
#const dnameRun = "../data/dataDidactic"  ; const fnameRun = "didactic2.txt"
const dnameRun  = "../data/dataFernandez" ; const fnameRun = "F50-56.txt"
#const dnameRun = "../data/dataHarris"    ; const fnameRun = "H10-2000.txt"
#const dnameRun  = "../data/dataBeasley" ; const fnameRun = "Capa-capb-red-90-50.txt"

# --- used when experiment = true (run all instances of one collection) ------
# uncomment one line; leave only one active
const dnameExperiment = "../data/dataFernandez"
#const dnameExperiment = "../data/dataHarris"
#const dnameExperiment = "../data/dataBeasley/"

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
include("voptMOA.jl")
include("biUFLP.jl")
include("reduce.jl")
include("mopRoutines.jl")
include("pavingBnB.jl")
include("reducePaving.jl")
include("labeling.jl")

if graphics
    include("graphicpyplot.jl")
end


# =============================================================================
function main()


    # 1. Run a didactic instance to compile all the solver --------------------

    println("  Compiling the solver...")
    if !epsilonCst
        # solve with the specific UFLP solver
        data = load2UFLP("../data/dataFernandez","F50-51.txt")
        biUFLPsolver(data)
    else
        # solve with the generic ϵ-constraint solver
        data = load2UFLP("../data/dataDidactic","didactic1.txt")
        epsilonConstraint(data,"GUROBI") 
    end


    # 2. Grab the folder with instances and instance(s) names ----------------- 

    if experiment
        # prepare the experiment on a list of instances 
        dname = dnameExperiment
        # get all the names of the datasets for a given collection available in folder 'dname'
        fnames = getfnames(dname)
        # matrix storing the results
        results = Array{Any, 2}(undef, length(fnames), 7)
    else
        dname   = dnameRun
        fnames  = [fnameRun]
        # matrix storing the results
        results = Array{Any, 2}(undef, 1, 7)
    end


    # 3. Call the solver ------------------------------------------------------

    for i in eachindex(fnames)

        verboseProd ? println("\n>>> [0] LOADING DATASET ") : nothing
        data  = load2UFLP(dname,fnames[i])

        if epsilonCst
            # Call the generic ϵ-constraint solver
            YN, timevOPt = epsilonConstraint(data,"GUROBI") 

            results[i,1] = fnames[i]
            results[i,6] = timevOPt
            results[i,7] = length(YN)
            if graphics
                thecolor = "green"; themarkerstyle = "o"; themarkersize = 15
                displayYN(YN,thecolor,themarkerstyle,themarkersize)

                mkpath(figDir)
                savefig(joinpath(figDir, fnames[1][1:end-4] * "_generic.png"))
            end

        else
            # Call the specific UFLP solver
            paving,timePaving,timeReducing,nbBoxPaving,nbBoxReducing,timeLabeling, ND_YN=biUFLPsolver(data)

            results[i,1] = fnames[i]
            results[i,2] = timePaving
            results[i,3] = timeReducing
            results[i,4] = nbBoxPaving
            results[i,5] = nbBoxReducing
            results[i,6] = timeLabeling
            results[i,7] = length(ND_YN)

            if graphics
                setupGraphic(fnames[1], data, paving)
                displayPaving(paving)
                #displayAllRunningCosts(paving)
                thecolor = "red"; themarkerstyle = "o"; themarkersize = 15
                displayYN(ND_YN,thecolor,themarkerstyle,themarkersize)

                mkpath(figDir)
                savefig(joinpath(figDir, fnames[1][1:end-4] * "_specific.png"))
            end

        end
    end
    println("  Run completed! \n")


    # 4. Display results ------------------------------------------------------

    if !epsilonCst
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
