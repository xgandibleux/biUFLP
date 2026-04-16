# ==============================================================================
# fname : biUFLP.jl
# August 2022
# ==============================================================================


# ==============================================================================
"""
    biUFLPsolver(data::Instance)

    Solver of bi-objective 0/1 uncapacited facility location problem.
    It principle is to build a paving by boxes on the objective space Y 
    which cover the nondominated points YN, and to analyse each box to 
    generate all y ∈ YN. 
"""
function biUFLPsolver(data::Instance)

    # -------------------------------------------------------------------------
    # Information about the instance to solve

    verboseProd ? println("\n    [1] INSTANCE \n") : nothing

    verboseProd ? println("        filename      : $(data.fname)") : nothing
    verboseProd ? println("        nI (users)    : $(data.nI)") : nothing
    verboseProd ? println("        nJ (services) : $(data.nJ)") : nothing


    # -------------------------------------------------------------------------
    # compute a paving of the instance

    verboseProd ? println("\n    [2] PAVING STAGE \n") : nothing

    getTime = time()
    paving = computePavingBranchAndBound(data)
    timePaving = round(time()- getTime, digits=4)
    nbBoxPaving = length(paving)
    verboseProd ? println("\n        Time(paving)  : ", timePaving, " sec \n") : nothing

    if verboseProd
        print("        J1 : ")
        for ib in reverse(eachindex(paving))
            print(paving[ib].J1)
        end
        println(" ")
    end

    # -------------------------------------------------------------------------
    # reduce a paving 

    verboseProd ? println("\n    [3] REFINING STAGE \n") : nothing

    getTime = time() 
    reducePaving!(data, paving)
    timeReducing = round(time()- getTime, digits=4)
    nbBoxReducing = length(paving)
    verboseProd ? println("\n        Time(refining)  : ", timeReducing, " sec \n") : nothing

    if verboseProd
        print("        J1 : ")
        for ib in reverse(eachindex(paving))
            print(paving[ib].J1)
        end
        println(" \n")
    end    

    # -------------------------------------------------------------------------
    # labeling the paving 

    verboseProd ? println("\n    [3] LABELING STAGE \n") : nothing
    
    getTime = time() 
    ND_YN, all_YN = labelingPaving!(data, paving)
    timeLabeling = round(time()- getTime, digits=4)
    verboseProd ? println("\n        Time(labeling)  : ", timeLabeling, " sec \n") : nothing    
    #@show paving
    #@show all_YN
    #@show ND_YN

    return paving, timePaving, timeReducing, nbBoxPaving, nbBoxReducing, timeLabeling, ND_YN
end
