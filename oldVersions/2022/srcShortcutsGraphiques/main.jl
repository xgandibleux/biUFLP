# ==============================================================================
# fname : main.jl
# August 2022
# ==============================================================================

println("Run biUFLP solver")
using Printf

# ==============================================================================
# parameters to control the app

const verboseProd  = true  # displays on the terminal the useful info in production
const verboseDev   = false  # displays on the terminal the development info
const verboseDev2  = false  # displays on the terminal the development info for phase 2 (labeling)
const experiment   = false  # run a full numerical experiment
const vOptSolver   = false  # run vOptGeneric with ϵ-constraint and glpk or gurobi
const graphics     = false   # display on the screen graphics 
const backendGR    = :PyPlot
#const backendGR   = :Luxor

# ==============================================================================
# files of the app

include("dataStru.jl")
include("parser.jl")
include("vopt.jl")
include("biUFLP.jl")
include("mopRoutines.jl")
include("pavingBnB.jl")
include("reducePaving.jl")
include("labeling.jl")

if graphics
    if backendGR == :PyPlot
        include("graphicpyplot.jl")
    elseif backendGR == :Luxor
        include("graphicluxor.jl")
    end
end


# ==============================================================================
function main()

    # prepare the run for a single resolution or a numerical experiment --------
    if experiment
        # prepare the experiment on a list of instances ------------------------

        #dname = "../data/dataFernandez"
        dname = "../data/dataHarris"

        # get all the names of the datasets for a given collection available in folder 'dnameE'
        fnames = getfnames(dname)

        # matrix storing the results
        results = Array{Any, 2}(undef, length(fnames), 7)

        # run a didactic instance to compile all the code 
        data = load2UFLP("../data/dataDidactic","didactic1.txt")
        if !vOptSolver
            biUFLPsolver(data)
        else
            vOpt(data,"GUROBI") # GLPK or GUROBI
        end


    else
        # setup for a single instance -----------------------------------------

        #dname = "../data/dataDidactic"
        #fnames = ["didactic1.txt"]
        #fnames = ["didactic2.txt"]

        #dname = "../data/dataFernandez"
        #fnames  = ["F50-51.txt"]
        #fnames  = ["F55-56.txt"]
        #fnames  = ["F54-57.txt"]

        dname = "../data/dataHarris"
        fnames  = ["H10-2000.txt"]
        

        # matrix storing the results
        results = Array{Any, 2}(undef, 1, 7)
    end

    # call the solver ---------------------------------------------------------
    for i in eachindex(fnames)
        verboseProd ? println("\n>>> [0] LOADING DATASET ") : nothing
    
        data  = load2UFLP(dname,fnames[i])
        if vOptSolver
            getTime = time()
            YN=vOpt(data,"GUROBI") # GLPK or GUROBI
            timevOPt = round(time()- getTime, digits=4)
            results[i,1] = fnames[i]        
            results[i,6] = timevOPt
            results[i,7] = length(YN)                             
            if graphics
                if backendGR == :PyPlot
                    thecolor = "green"; themarkerstyle = "o"; themarkersize = 15
                    displayYN(YN,thecolor,themarkerstyle,themarkersize)
                end
            end
        else
            paving,timePaving,timeReducing,nbBoxPaving,nbBoxReducing,timeLabeling, ND_YN=biUFLPsolver(data) # ajouter le retour avec YN et XE lors integration dans vOptSolver
            results[i,1] = fnames[i]
            results[i,2] = timePaving
            results[i,3] = timeReducing                        
            results[i,4] = nbBoxPaving
            results[i,5] = nbBoxReducing   
            results[i,6] = timeLabeling
            results[i,7] = length(ND_YN) 

            if graphics
                if backendGR == :PyPlot
                    setupGraphic(fnames[1], data, paving)
                    displayPaving(paving)
                    #displayAllRunningCosts(paving)
                    thecolor = "green"; themarkerstyle = "o"; themarkersize = 15
                    displayYN(ND_YN,thecolor,themarkerstyle,themarkersize)

                elseif backendGR == :Luxor
                    # draw the results (currently works only for didactic instance)
                    YN=[]
                    drawGraphics(paving,YN)
                end
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

if graphics
    if backendGR == :Luxor
        finish()
        preview()
    end
end