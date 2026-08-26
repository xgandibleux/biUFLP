# ==============================================================================
# fname : main.jl
# August 2022
# ==============================================================================

include("dataStru.jl")
include("parser.jl")
include("model.jl")
include("paving.jl")
include("box.jl")
include("filtering.jl")
include("supportedPtsBox.jl")
include("labeling.jl")
include("ensBornant.jl")


# ==============================================================================
# parametres de controle des affichages

const draw_EBpointsIdeaux     = false
const draw_nondominatedPoints = true
const draw_allRunningCosts    = true
const draw_paving             = true
const draw_CRplusCA12andCA21  = true
const draw_yNS                = false
const draw_decompositionBox   = false
const draw_nondominatedPointsALGO = true


# ==============================================================================
# Load an instance of 2UFLP
dname = "../data/dataDidactic"
#dname = "../data/dataFernandez"
fname = "didactic1.txt"
#fname = "F55-56.txt"
data  = load2UFLP(dname,fname)


# ==============================================================================
# Compute Y_N with vOptSolver

println("\n\n>>> COMPUTATION OF Y_N WITH VOPTSOLVER \n\n")

# Define the 2UFLP model
mod2UFLP = createProblem2UFLP(GLPK.Optimizer, data)

# Optimize the model with vOptSolver using an ϵ-constraint method
vSolve(mod2UFLP, method=:epsilon, step = 1.0)

# Get the results
Y_N = getY_N( mod2UFLP )


# ==============================================================================
# compute an exhaustive paving of the instance

println("\n\n>>> EXHAUSTIVE GENERATION STAGE \n\n")
paving = generateExhaustivePaving(data)


# ==============================================================================
# filtering the paving

println("\n\n>>> FILTERING STAGE \n\n")
filteringPaving!(data, paving)


# ==============================================================================
# reducing and pruning the paving 

function reducingPaving!(data::Instance, paving::Vector{Box})

    # vector to mark the pruned boxes; the dominated boxes are physically deleted after the loop
    pruned::Vector{Bool} = fill(false,length(paving))
    nbPruned ::Dict{Symbol, Int64} = Dict(:REDUCTION => 0)

    # sort the boxes on base of yI[1] in the perspective of selecting relevant boxes (tri non encore utilise dans la strategie)
    #sort!(paving, by = v -> v.yI[1])

    # if it is not already done, generate some yNS for the box (verifier si listyNS vide peut vraiment apparaitre apres filtrage)
    for ib in eachindex(paving)
        if length(paving[ib].listyNS) == 0
            println("  Preliminary stage: generate some yNS for box $ib")
            computeyNS!(data, paving[ib])
            sort!(paving[ib].listyNS, by = v -> v[1]) # yNS are sorted by increasing values of f1
        end
    end

    # consider two by two all the boxes currently in the paving ---------------
    
    for ib1 in eachindex(paving), ib2 in eachindex(paving)

        if (ib1≠ib2) && (!pruned[ib1]) && (!pruned[ib2]) 
            # ib1 and ib2 are 2 indexes on 2 different boxes, neither ib1 nor ib2 are already pruned

            println("----------------------------------------------------")
            println(length(paving), " $ib1 $ib2 : ")

            # try to reduce and prune the box B2 
            tryToReducePruneBox!(data, paving, ib1, ib2, pruned, nbPruned)

        end
    end


    # Summary -----------------------------------------------------------------
    # 1) Remaining boxes (active boxes) in the paving
    println("\nSummary: ")
    print("  Active box(es): ")
    for i in eachindex(pruned)
        if !pruned[i]
            print("$i  ")
        end
    end
    println(" ")

    # 2) Effectiveness of the reduction/pruning procedure
    nbBoxes=length(pruned)
    percent = round((nbPruned[:REDUCTION])/nbBoxes*100; digits=2)
    println("  ",nbBoxes-nbPruned[:REDUCTION]," box(es) active || pruned: ",nbPruned[:REDUCTION]," box(es) | $percent %)\n\n")


    # Delete from the paving all the boxes useless ----------------------------
    deleteBoxes!(paving, pruned)

end


# REDUCTION NO+SE DE LA BOITE B2 SUR DOMINATION DE SUPPORTES DE B2 PAR SUPPORTES DE B1
# ET FILTRAGE SI UN SUPPORTE OU OPTIMA LEX DE B1 DOMINE UN IDEAL DE LA BOITE B2 

function tryToReducePruneBox!(data::Instance, paving::Vector{Box}, ib1::Int64, ib2::Int64, pruned::Vector{Bool}, nbPruned::Dict{Symbol, Int64})

    # -------------------------------------------------------------------------
    # Preliminary stage -------------------------------------------------------
    # -------------------------------------------------------------------------

    println(paving[ib1]);
    println(paving[ib2]);
    println(">>>> length(paving[ib1].listyNS) ", length(paving[ib1].listyNS), "  length(paving[ib2].listyNS) ",length(paving[ib2].listyNS))

    # -------------------------------------------------------------------------
    # Tentative of reduction of the box 2 by the North-West -------------------
    # -------------------------------------------------------------------------

    println("\nReduction by the North-West")

    # check if y12[1] for box 2 is on the right of y12[1] for box 1 
    if paving[ib2].y12[1] >= paving[ib1].y12[1] 

        # trying the reduction by NW is then allowed 

        ip1 = 1                  # iterator on the points yNS for box 1 set to the first
        stopReductionNW = false  # set to true if the reduction procedure by NW must be stopped 
        boxReduced = false       # set to true if the box has been reduced

        # while the condition for reducing are verified and the lists of yNS for boxes 1 & 2 are not empty
        while (!stopReductionNW) && (length(paving[ib1].listyNS)!=0) && (length(paving[ib2].listyNS)!=0)

            # test between a point yNS from box 1 and the first yNS from box 2
            if isDominates(paving[ib1].listyNS[ip1],paving[ib2].listyNS[begin])

                # domination stated => application of the reduction of the box 2

                println("Box2 reduced: $(paving[ib2].y12) ⇢ $(paving[ib2].listyNS[begin])")
                boxReduced = true
                paving[ib2].y12   = paving[ib2].listyNS[begin]    # update the z12 value of box 2
                paving[ib2].yI[1] = paving[ib2].listyNS[begin][1] # update the zI value of box 2
                paving[ib2].yN[2] = paving[ib2].listyNS[begin][2] # update the zN value of box 2
                popfirst!(paving[ib2].listyNS)                    # discard from the list the yNS used for updating z12 and continue with the next yNS of box 2

                if length(paving[ib2].listyNS)==0
                    # no remaining yNS point to test in the list
                    stopReductionNW = true
                end

            else
                # no domination stated => continue with the next yNS point of the box 1

                ip1+=1
                if ip1 > length(paving[ib1].listyNS)
                    # all the points yNS point of the box 1 have been tried => stop the NW reduction
                    stopReductionNW = true 

                elseif paving[ib1].listyNS[ip1][1] > paving[ib2].listyNS[begin][1]
                    # cannot reduce anymore the NW part of box 2 (current yNS[1] from box 1 is on the right of yNS[1] from box 2) => stop the NW reduction
                    stopReductionNW = true
                end
            end

        end          

    end # trying the reduction by NW

    # -------------------------------------------------------------------------
    # Tentative of pruning the box 2 (reduced or not) by box 1 ----------------
    # -------------------------------------------------------------------------

    # test if yI of b2 can be dominated now by one y12 or y21 of b1, and then to prune entirely b2
    if isDominates(paving[ib1].y12, paving[ib2].yI) || isDominates(paving[ib1].y21, paving[ib2].yI)
        # found!
        println("  pruning condition (y12/y21) found; box ib=$ib2 fully pruned !!!")
        pruned[ib2] = true
        nbPruned[:REDUCTION]+=1

    else  
        # test if yI of b2 can be dominated now by one yNS ∈ listyNS of b1, and then to prune entirely b2

        for ip in eachindex(paving[ib1].listyNS)
            if isDominates(paving[ib1].listyNS[ip], paving[ib2].yI)
                # found one yNS of b1 who D/W/= yI of b2
                println("  pruning condition (yNS) found; box ib=$ib2 fully pruned !!!")
                pruned[ib2] = true
                nbPruned[:REDUCTION]+=1
                break # dirty way to leave the for-loop when b2 is pruned
            end
        end

    end

    # -------------------------------------------------------------------------
    # Tentative of reduction of the box 2 by the South-East -------------------
    # -------------------------------------------------------------------------

    # check if the box 2 has not been pruned after the tentative of reduction by the NW
    if pruned[ib2] == false

        println("\nReduction by the South-East")

        # check if y21[2] for box 2 is over the bottom of y21[2] for box 1 
        if paving[ib2].y21[2] >= paving[ib1].y21[2]

            # trying the reduction by SE is then allowed 

            ip1 = length(paving[ib1].listyNS) # iterator on the points yNS for box 1 set to the last
            stopReductionSE = false           # set to true if the reduction procedure by NW must be stopped
            boxReduced = false                # set to true if the box has been reduced   
        
            # while the condition for reducing are verified and the lists of yNS for boxes 1 & 2 are not empty
            while (!stopReductionSE) && (length(paving[ib1].listyNS)!=0) && (length(paving[ib2].listyNS)!=0)    

                # test between a point yNS from box 1 and the last yNS from box 2    
                if isDominates(paving[ib1].listyNS[ip1],paving[ib2].listyNS[end])

                    # domination stated => reduction of the box 2

                    println("Box2 reduced: $(paving[ib2].y21) ⇢ $(paving[ib2].listyNS[end])")
                    boxReduced = true
                    paving[ib2].y21   = paving[ib2].listyNS[end]    # update the z21 value of box 2
                    paving[ib2].yI[2] = paving[ib2].listyNS[end][2] # update the zI value of box 2
                    paving[ib2].yN[1] = paving[ib2].listyNS[end][1] # update the zN value of box 2
                    pop!(paving[ib2].listyNS)                       # discard from the list the yNS used for updating z21 and continue with the previous yNS of box 2

                    if length(paving[ib2].listyNS)==0
                        # no remaining yNS point to test in the list
                        stopReductionSE = true
                    end

                else
                    # no domination stated => continue with the previous yNS point of the box 1

                    ip1-=1
                    if ip1 == 0
                        # all the points yNS point of the box 1 have been tried => stop the NW reduction
                        stopReductionSE = true 

                    elseif paving[ib1].listyNS[ip1][2] > paving[ib2].listyNS[end][2]
                        # cannot reduce anymore the SE part of box 2 (current yNS[1] from box 1 is over of yNS[1] from box 2) => stop the SE reduction
                        stopReductionSE = true
                    end
                end

            end

        end # trying the reduction by SE

        # -------------------------------------------------------------------------
        # Tentative of pruning the box 2 (reduced or not) by box 1 ---------------- (passage de code identique -> factoriser eventuellement)
        # -------------------------------------------------------------------------

        # test if yI of b2 can be dominated now by one y12 or y21 of b1, and then to prune entirely b2
        if isDominates(paving[ib1].y12, paving[ib2].yI) || isDominates(paving[ib1].y21, paving[ib2].yI)
            # found!
            println("  pruning condition found; box ib=$ib2 fully pruned (y12 or y21) !!! ")
            pruned[ib2] = true
            nbPruned[:REDUCTION]+=1

        else  
            # test if zI of b2 can be dominated now by one yNS ∈ listyNS of b1, and then to prune entirely b2
        
            for ip in eachindex(paving[ib1].listyNS)
                if isDominates(paving[ib1].listyNS[ip], paving[ib2].yI)
                    # found one yNS of b1 who D/W/= yI of b2
                    println("  pruning condition found; box ib=$ib2 fully pruned (yNS) !!!")
                    pruned[ib2] = true
                    nbPruned[:REDUCTION]+=1
                    break # dirty way to leave the for-loop when b2 is pruned
                end
            end
        end

    end    

end

println("\n\n>>> REDUCING+PRUNING STAGE \n\n")

reduction = true
if reduction    
    reducingPaving!(data, paving)
end # reduction


# ==============================================================================
# compute the bound set based on yI

ND_yI, BoundSet_yI = elaborateBoundSet_yI(paving) 


# ==============================================================================
# ==============================================================================
# labeling 
ND_YN,all_YN = labelingPaving!(data, paving)


# =============================================================================
# =============================================================================

# version (incomplete) avec une matrice

#=
# intitialize the matrix of labels
labeling = Matrix{Union{Nothing,Vector{Label}}}(nothing, length(J1),data.nI)
#labeling[1,1] = [Label()]

# 1st user => (0,0) + (c1,c2)
i=1
for (pos,j) in enumerate(J1nd[i])
    # compute the labels for user i - service j∈J1nd
    c1 = data.c1[i,j]
    c2 = data.c2[i,j]
    println("$i $j [$c1 $c2 ; $j]")
    labeling[j,i] = [Label([c1,c2],[j])]
end

for i in 2:data.nI
    # compute the labels for the user i
    for (pos,j) in enumerate(J1nd[i])
        # compute the labels for user i - service j∈J1nd
        println("$i $j")
    end
end

# display the matrix of labels
for j in 1:length(J1nd)
    for i in 1:data.nI
        if labeling[j,i] === nothing
            print("($j,$i) x   ")
        else 
            print("($j,$i) ",length(labeling[j,i]) , "   ")
        end
    end
    println(" ")
end
=#

# ==============================================================================
# display graphically the results


backend = :pyplot
#backend = :luxor

if backend == :pyplot
    include("graphicpyplot.jl")
elseif backend == :luxor
    include("graphicluxor.jl")
end