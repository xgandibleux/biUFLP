# ==============================================================================
# fname : reducePaving.jl
# August 2022
# ==============================================================================


# =============================================================================
# actuellement fonction non appelee 
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
    tryToReducePruneBox!(paving::Vector{Box}, ib1::Int64, ib2::Int64, pruned::Vector{Bool}, nbPruned::Dict{Symbol, Int64})

    Given two boxes indexed by ib1 and ib2, try to reduce box 2 and after check if b2 can be pruned
    - The reduction is first tempted on the North-West part of b2, and next on the South-East part of b2.
    - The reduction exploits some nondominated supported points (YNS) of b1 and check if y ∈ YNS of b2 are dominated
    - Box 2 is pruned is exists y12,y21 ∈ YNS of b1 or a y ∈ YNS of b1 who dominates yI of b2
"""
function tryToReducePruneBox!(paving::Vector{Box}, ib1::Int64, ib2::Int64, pruned::Vector{Bool}, nbPruned::Dict{Symbol, Int64})

    verboseDev ? println(paving[ib1]) : nothing
    verboseDev ? println(paving[ib2]) : nothing
    verboseDev ? println(">>>> length(paving[ib1].listyNS) ", length(paving[ib1].listyNS), "  length(paving[ib2].listyNS) ",length(paving[ib2].listyNS)) : nothing

    # -------------------------------------------------------------------------
    # Tentative of reduction of the box 2 by the North-West -------------------
    # -------------------------------------------------------------------------

    verboseDev ? println("\nReduction by the North-West") : nothing

    # check if y12[1] for box 2 is on the right of y12[1] for box 1 
    if paving[ib2].y12[1] >= paving[ib1].y12[1] 

        # trying the reduction by NW is then allowed 

        ip1 = 1                  # iterator on the points yNS for box 1 set to the first
        stopReductionNW = false  # set to true if the reduction procedure by NW must be stopped 
        boxReduced = false       # set to true if the box has been reduced

        # while the condition for reducing are verified and the lists of yNS for boxes 1 & 2 are not empty
        while (!stopReductionNW) && (length(paving[ib1].listyNS)!=0) && (length(paving[ib2].listyNS)!=0)

            # test between a point yNS from box 1 and the local ideal between z12 and the first yNS from box 2
            yIlocal = Vector{Int64}(undef, 2);   
            yIlocal[1] = paving[ib2].y12[1];   yIlocal[2] = paving[ib2].listyNS[begin][2] 

            if isDominates(paving[ib1].listyNS[ip1],yIlocal)             
#            if isDominates(paving[ib1].listyNS[ip1],paving[ib2].listyNS[begin]) # correction apportee

                # domination stated => application of the reduction of the box 2

                verboseDev ? println("Box2 reduced: $(paving[ib2].y12) ⇢ $(paving[ib2].listyNS[begin])") : nothing
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

                elseif paving[ib1].listyNS[ip1][1] > paving[ib2].y12[1]
                    # cannot reduce anymore the NW part of box 2 (current yNS[1] from box 1 is on the right of y12[1] from box 2) => stop the NW reduction
                    stopReductionNW = true                    

                elseif paving[ib1].listyNS[ip1][1] > paving[ib2].listyNS[begin][1] # encore activable avec le test precedent?
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
        verboseDev ? println("  pruning condition (y12/y21) found; box ib=$ib2 fully pruned !!!") : nothing
        pruned[ib2] = true
        nbPruned[:REDUCTION]+=1

    else  
        # test if yI of b2 can be dominated now by one yNS ∈ listyNS of b1, and then to prune entirely b2

        for ip in eachindex(paving[ib1].listyNS)
            if isDominates(paving[ib1].listyNS[ip], paving[ib2].yI)
                # found one yNS of b1 who D/W/= yI of b2
                verboseDev ? println("  pruning condition (yNS) found; box ib=$ib2 fully pruned !!!") : nothing
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

        verboseDev ? println("\nReduction by the South-East") : nothing

        # check if y21[2] for box 2 is over the bottom of y21[2] for box 1 
        if paving[ib2].y21[2] >= paving[ib1].y21[2]

            # trying the reduction by SE is then allowed 

            ip1 = length(paving[ib1].listyNS) # iterator on the points yNS for box 1 set to the last
            stopReductionSE = false           # set to true if the reduction procedure by NW must be stopped
            boxReduced = false                # set to true if the box has been reduced   
        
            # while the condition for reducing are verified and the lists of yNS for boxes 1 & 2 are not empty
            while (!stopReductionSE) && (length(paving[ib1].listyNS)!=0) && (length(paving[ib2].listyNS)!=0)    
                
                # test between a point yNS from box 1 and the local ideal between the first yNS from box 2 and z21
                yIlocal = Vector{Int64}(undef, 2);   
                yIlocal[1] = paving[ib2].listyNS[end][1]; yIlocal[2] = paving[ib2].y21[2];   
 
                if isDominates(paving[ib1].listyNS[ip1],yIlocal)                 
                #if isDominates(paving[ib1].listyNS[ip1],paving[ib2].listyNS[end])  # correction apportee

                    # domination stated => reduction of the box 2

                    verboseDev ? println("Box2 reduced: $(paving[ib2].y21) ⇢ $(paving[ib2].listyNS[end])") : nothing
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

                    elseif paving[ib1].listyNS[ip1][2] > paving[ib2].y21[2]
                        # cannot reduce anymore the SE part of box 2 (current yNS[2] from box 1 is over of y21[2] from box 2) => stop the SE reduction
                        stopReductionSE = true

                    elseif paving[ib1].listyNS[ip1][2] > paving[ib2].listyNS[end][2]
                        # cannot reduce anymore the SE part of box 2 (current yNS[2] from box 1 is over of yNS[2] from box 2) => stop the SE reduction
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
            verboseDev ? println("  pruning condition found; box ib=$ib2 fully pruned (y12 or y21) !!! ") : nothing
            pruned[ib2] = true
            nbPruned[:REDUCTION]+=1

        else  
            # test if zI of b2 can be dominated now by one yNS ∈ listyNS of b1, and then to prune entirely b2
        
            for ip in eachindex(paving[ib1].listyNS)
                if isDominates(paving[ib1].listyNS[ip], paving[ib2].yI)
                    # found one yNS of b1 who D/W/= yI of b2
                    verboseDev ? println("  pruning condition found; box ib=$ib2 fully pruned (yNS) !!!") : nothing
                    pruned[ib2] = true
                    nbPruned[:REDUCTION]+=1
                    break # dirty way to leave the for-loop when b2 is pruned
                end
            end
        end

    end    

end


# ==============================================================================
"""
    tryToReducePruneBoxUnique!(paving::Vector{Box}, ib1::Int64, ib2::Int64, pruned::Vector{Bool}, nbPruned::Dict{Symbol, Int64})

    Given two boxes indexed by ib1 and ib2, try to reduce box 2 and after check if b2 can be pruned
    - The reduction is first tempted on the North-West part of b2, and next on the South-East part of b2.
    - The reduction exploits y12 of b1 and check if y ∈ YNS of b2 are dominated
    - Box 2 is pruned if y12,y21 of b1 dominates yI of b2
"""
function tryToReducePruneBoxUnique!(paving::Vector{Box}, ib1::Int64, ib2::Int64, pruned::Vector{Bool}, nbPruned::Dict{Symbol, Int64})

    verboseDev ? println(paving[ib1]) : nothing
    verboseDev ? println(paving[ib2]) : nothing
    verboseDev ? println(">>>> UNIQUE length(paving[ib1].listyNS) ", length(paving[ib1].listyNS), "  length(paving[ib2].listyNS) ",length(paving[ib2].listyNS)) : nothing

    # -------------------------------------------------------------------------
    # Tentative of reduction of the box 2 by the North-West -------------------
    # -------------------------------------------------------------------------

    verboseDev ? println("\nReduction by the North-West") : nothing

    # check if y12[1] for box 2 is on the right of y12[1] for box 1 
    if paving[ib2].y12[1] >= paving[ib1].y12[1] 

        # trying the reduction by NW is then allowed 

        stopReductionNW = false  # set to true if the reduction procedure by NW must be stopped 
        boxReduced = false       # set to true if the box has been reduced

        # while the condition for reducing are verified and the lists of yNS for box 2 is not empty
        while (!stopReductionNW) && (length(paving[ib2].listyNS)!=0)

            # test between point y12 from box 1 and the local ideal between z12 and the first yNS from box 2
            yIlocal = Vector{Int64}(undef, 2);   
            yIlocal[1] = paving[ib2].y12[1];   yIlocal[2] = paving[ib2].listyNS[begin][2] 

            if isDominates(paving[ib1].y12,yIlocal)             

                # domination stated => application of the reduction of the box 2

                verboseDev ? println("Box2 reduced: $(paving[ib2].y12) ⇢ $(paving[ib2].listyNS[begin])") : nothing
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
                # no domination stated => stop                
                stopReductionNW = true
            end

        end          

    end # trying the reduction by NW

    # -------------------------------------------------------------------------
    # Tentative of pruning the box 2 (reduced or not) by box 1 ----------------
    # -------------------------------------------------------------------------

    # test if yI of b2 can be dominated now by one y12 or y21 of b1, and then to prune entirely b2
    if isDominates(paving[ib1].y12, paving[ib2].yI)
        # found!
        verboseDev ? println("  pruning condition (y12) found; box ib=$ib2 fully pruned !!!") : nothing
        pruned[ib2] = true
        nbPruned[:REDUCTION]+=1
    end

    # -------------------------------------------------------------------------
    # Tentative of reduction of the box 2 by the South-East -------------------
    # -------------------------------------------------------------------------

    # check if the box 2 has not been pruned after the tentative of reduction by the NW
    if pruned[ib2] == false

        verboseDev ? println("\nReduction by the South-East") : nothing

        # check if y21[2] for box 2 is over the bottom of y21[2] for box 1 
        if paving[ib2].y21[2] >= paving[ib1].y21[2]

            # trying the reduction by SE is then allowed 

            stopReductionSE = false           # set to true if the reduction procedure by NW must be stopped
            boxReduced = false                # set to true if the box has been reduced   
        
            # while the condition for reducing are verified and the lists of yNS for boxes 1 & 2 are not empty
            while (!stopReductionSE) && (length(paving[ib2].listyNS)!=0)    
                
                # test between a point yNS from box 1 and the local ideal between the first yNS from box 2 and z21
                yIlocal = Vector{Int64}(undef, 2);   
                yIlocal[1] = paving[ib2].listyNS[end][1]; yIlocal[2] = paving[ib2].y21[2];   
 
                if isDominates(paving[ib1].y12,yIlocal)                 

                    # domination stated => reduction of the box 2

                    verboseDev ? println("Box2 reduced: $(paving[ib2].y21) ⇢ $(paving[ib2].listyNS[end])") : nothing
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
                    # no domination stated => stop
                    stopReductionSE = true

                end

            end

        end # trying the reduction by SE

        # -------------------------------------------------------------------------
        # Tentative of pruning the box 2 (reduced or not) by box 1 ---------------- (passage de code identique -> factoriser eventuellement)
        # -------------------------------------------------------------------------

        # test if yI of b2 can be dominated now by one y12 of b1, and then to prune entirely b2
        if isDominates(paving[ib1].y12, paving[ib2].yI) 
            # found!
            verboseDev ? println("  pruning condition found; box ib=$ib2 fully pruned (y12 or y21) !!! ") : nothing
            pruned[ib2] = true
            nbPruned[:REDUCTION]+=1
        end

    end    

end


# ==============================================================================
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
    reducePaving!(data::Instance, paving::Vector{Box})

    Given a paving, try to reduce all the boxes and prune useless boxes
"""
function reducePaving!(data::Instance, paving::Vector{Box})

    # vector to mark the pruned boxes; the dominated boxes are physically deleted after the loop
    pruned::Vector{Bool} = fill(false,length(paving))
    nbPruned ::Dict{Symbol, Int64} = Dict(:REDUCTION => 0)

    # sort the boxes on base of yI[1] in the perspective of selecting relevant boxes (tri non encore utilise dans la strategie)
    #sort!(paving, by = v -> v.yI[1])

    # if it is not already done, generate some yNS for the box (verifier si listyNS vide peut vraiment apparaitre apres filtrage)
    for ib in eachindex(paving)
        if length(paving[ib].listyNS) == 0
            verboseDev ? println("  Preliminary stage: generate some yNS for box $ib") : nothing
            computeyNS!(data, paving[ib])
            sort!(paving[ib].listyNS, by = v -> v[1]) # yNS are sorted by increasing values of f1
        end
    end

    # consider two by two all the boxes currently in the paving ---------------
    
    for ib1 in eachindex(paving), ib2 in eachindex(paving)

        if (ib1≠ib2) && (!pruned[ib1]) && (!pruned[ib2]) 
            # ib1 and ib2 are 2 indexes on 2 different boxes, neither ib1 nor ib2 are already pruned

            verboseDev ? println("----------------------------------------------------") : nothing
            verboseDev ? println(length(paving), " $ib1 $ib2 : ") : nothing

            # try to reduce and prune the box B2 
            # les 2 fonctions suivantes pourraient faire qu'une seule pour factoriser le code au prix de quelques tests supplementaires
            if length(paving[ib1].J1)==1
                # unique service (=> no supported points available in the "box")
                tryToReducePruneBoxUnique!(paving, ib1, ib2, pruned, nbPruned)
            else
                tryToReducePruneBox!(paving, ib1, ib2, pruned, nbPruned)
            end

        end
    end


   # summary of the activity -------------------------------------------------  
   nbBoxes=length(pruned)
   verboseProd ? println("        #initial boxes          : ", nbBoxes) : nothing 
   verboseProd ? println("        #pruned boxes           : ", nbPruned[:REDUCTION]) : nothing  
   verboseProd ? println("        #active boxes           : ", nbBoxes-nbPruned[:REDUCTION]) : nothing       
   verboseProd ? println("        ratio of boxes pruned   : ", round((nbPruned[:REDUCTION])/nbBoxes*100; digits=2)  , "%") : nothing  


    # Delete from the paving all the boxes useless ----------------------------
    deleteBoxes!(paving, pruned)

end