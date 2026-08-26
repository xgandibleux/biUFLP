# =============================================================================
# filtering.jl
# August 2022
# =============================================================================

# =============================================================================
"""
    isOverlap(b1::Box, b2::Box)

    Test if an overlap exists between 2 boxes
"""
function isOverlap(b1::Box, b2::Box)

    if (b1.y12[1] >= b2.y21[1] || b2.y12[1] >= b1.y21[1])
        # If one rectangle is on left side of other
        return false

    elseif (b1.y21[2] >= b2.y12[2] || b2.y21[2] >= b1.y12[2])
        # If one rectangle is above other
        return false

    else
        return true
    end

end


# ==============================================================================
"""
    isDominates(yA::Vector{Int64}, yB::Vector{Int64})

    Test if yA stronglyDominates/weaklyDominates/isEqualTo yB
"""
function isDominates(yA::Vector{Int64}, yB::Vector{Int64})

    #print("  yA=$yA  yB=$yB ⟶  ")

    if (yA[1] < yB[1]) && (yA[2] < yB[2])
        # strict dominance between yA and yB
        #println("yA=$(yA) D yB=$(yB)")
        return true

    elseif ((yA[1] ≤ yB[1]) && (yA[2] < yB[2])) || ((yA[1] < yB[1]) && (yA[2] ≤ yB[2]))
        # weak dominance between yA and yB
        #println("yA=$(yA) W yB=$(yB)")
        return true        

    elseif (yA[1] == yB[1]) && (yA[2] == yB[2])
        # egality between yA and yB
        #println("yA=$(yA) = yB=$(yB)")  
        return true    

    else
        #println("no D|W|= of yA over yB")
        return false
    end
end


# ==============================================================================
"""
    tryToPruneBox!(data::Instance, paving::Vector{Box}, ib1::Int64, ib2::Int64, pruned::Vector{Bool}, nbPruned::Dict{Symbol, Int64})

    Given two boxes B1 and B2 indexed by ib1 and ib2, try to prune the box B2 by applying hierarchically two strategies: 
      (1) trying with y12,y21 ∈ b1 and 
      (2) trying with a yNS ∈ b1
"""
function tryToPruneBox!(data::Instance, paving::Vector{Box}, ib1::Int64, ib2::Int64, pruned::Vector{Bool}, nbPruned::Dict{Symbol, Int64})

    # Try to prune entirely b2 with y12,y21 ∈ b1 ------------------------------

    # Test if yI ∈ b2 is in the cone of one feasible point y12,y21 ∈ b1 
    if (isDominates(paving[ib1].y12, paving[ib2].yI)) || (isDominates(paving[ib1].y21, paving[ib2].yI))
        # found y12 or y21 of b1 who D/W/= yI of b2
        println("  pruning condition by y12 or y21 found => $ib2 fully pruned !!!")
        pruned[ib2] = true
        nbPruned[:OPTIMALEX]+=1

    # try to prune entirely b2 with yNS ∈ b1 ----------------------------------
    
    elseif isOverlap(paving[ib1],paving[ib2])

        println("Overlap")

        # if it is not already done, generate some yNS for b1
        if length(paving[ib1].listyNS) == 0
            println("  generate some yNS for box $ib1")
            computeyNS!(data, paving[ib1])  # depthlMax is a global constant
            sort!(paving[ib1].listyNS, by = v -> v[1]) # yNS are sorted by increasing values of f1
        end

        # Try to prune entirely b2 with one yNS ∈ listyNS
        for ip in eachindex(paving[ib1].listyNS)
            #  Test if yI ∈ b2 is in the cone of one feasible point yNS ∈ listyNS
            if isDominates(paving[ib1].listyNS[ip], paving[ib2].yI)
                # found one yNS of b1 who D/W/= yI of b2
                println("  pruning condition found; $ib2 fully pruned")
                pruned[ib2] = true
                nbPruned[:SUPPORTED]+=1
                break # dirty way to leave the for-loop when b2 is pruned
            end
        end

    # nothing to claim between b1 and b2 at this stage ------------------------
    else
        println("---")
    end
end


"""
    deleteBoxes!(paving::Vector{Box}, pruned::Vector{Bool})

    Delete from the paving all the boxes useless
"""
function deleteBoxes!(paving::Vector{Box}, pruned::Vector{Bool})

    for i in reverse(1:length(pruned))
        if pruned[i]
            #println("box pruned : $i")
            deleteat!(paving,i)
        end
    end
end


# ==============================================================================
"""
    filteringPaving!(data::Instance, paving::Vector{Box})

    Refine the paving in filtering useless boxes thanks to 2 tests applied hierarchically
"""
function filteringPaving!(data::Instance, paving::Vector{Box})

    # initialize the vector to mark the pruned boxes; the dominated boxes are physically deleted after the loop
    pruned=fill(false,length(paving))
    nbPruned = Dict(:OPTIMALEX => 0 , :SUPPORTED => 0)

    # consider two by two all the boxes currently in the paving ---------------
    for ib1 in eachindex(paving), ib2 in eachindex(paving)

        if (ib1≠ib2) && (!pruned[ib1]) && (!pruned[ib2]) 
            # ib1 and ib2 are 2 indexes on 2 different boxes, neither ib1 nor ib2 are already pruned

            print("$ib1 $ib2 : ")

            # try to prune the box B2 
            # 1) if a lex optimal point from box 1 dominates the ideal point of Box2
            # 2) when the two boxes are overlapping, if a supported point of box 1 dominates the ideal point of Box2
            tryToPruneBox!(data, paving, ib1, ib2, pruned, nbPruned)

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

    # 2) Effectiveness of the filtering procedure
    nbBoxes=length(pruned)
    percent = round((nbPruned[:OPTIMALEX]+nbPruned[:SUPPORTED])/nbBoxes*100; digits=2)
    println("  ",nbBoxes-(nbPruned[:OPTIMALEX]+nbPruned[:SUPPORTED])," box(es) active || pruned: ",nbPruned[:OPTIMALEX]," box(es)/lex | ",nbPruned[:SUPPORTED]," box(es)/yNS | $percent %)")


    #for ib1 in eachindex(paving)
    #    print(">>>> length(paving[ib1].listyNS) $ib1 : ", length(paving[ib1].listyNS))
    #    if pruned[ib1] println(" pruned") else println(" ") end
    #end

    # Delete from the paving all the boxes useless ----------------------------
    deleteBoxes!(paving, pruned)

end