# Version with a matrix as datastructure

const withTraceAllocation = false


# ==============================================================================
function layer(i::Int64)
    return mod(i+1,2)+1
end


# ==============================================================================
function computeListNondominatedAllocationJ1oneUser(data::Instance, iUser::Int64, J1::Vector{Int64})

    # For a given user ('iUser') and a list J1 of services to consider, assemble in a matrix 
    # the 2D allocation costs user-service + its index of service. 
    # convention  ->  line :: user i  | column ::  1=>c1(i,j)   2=>c2(i,j)   3=>j 
    vctCosts = Array{Int64,2}(undef,length(J1),3) 
    for j in 1:length(J1)
        vctCosts[j,1] = data.c1[iUser,J1[j]]
        vctCosts[j,2] = data.c2[iUser,J1[j]]
        vctCosts[j,3] = J1[j] 
    end

    # sort on dim 1 (lines) the matrix 'vctCosts'  => lexicographically sort
    vctCosts = sortslices(vctCosts,dims=1)

    # list of indexes of nondominated couples
    J1nd = (Int)[vctCosts[1,3]] # after the sort, the first couple is not dominated => save the value of service j 

    # save the index for couple of costs nondominated 
    minpoints_y = vctCosts[1,2]
    for j in 2:length(J1)
        if vctCosts[j,2] < minpoints_y
            # nondominated
            minpoints_y = vctCosts[j,2]
            push!(J1nd,vctCosts[j,3])      
        end
    end

    return J1nd
end


# ==============================================================================
function computeListNondominatedAllocationJ1allUsers(data::Instance, J1::Vector{Int64})

    vctNDalloc = Vector{Vector{Int64}}(undef,data.nI)

    # for all the users, compute the list of services corresponding to non-dominated allocation
    for iUser in 1:data.nI

        if verboseDev2
            println("i = $iUser")
            for j in 1:length(J1)
                println("$(J1[j])  ( $(data.c1[iUser,J1[j]])  $(data.c2[iUser,J1[j]]) ) ") 
            end
            println(" ")
        end

        vctNDalloc[iUser] = computeListNondominatedAllocationJ1oneUser(data, iUser, J1)

        # summary -----------------------------------------------------------------
        if verboseDev2
            J1nd = vctNDalloc[iUser]
            if length(vctNDalloc[iUser])==1  
                println("$(J1nd[begin])  ( $(data.c1[iUser,J1nd[begin]])  $(data.c2[iUser,J1nd[begin]]) ) dominating all others") 
            else
                for (pos,j) in enumerate(vctNDalloc[iUser])
                    # pos : position in the list -> not used
                    println("$j  ( $(data.c1[iUser,j])  $(data.c2[iUser,j]) ) nondominated") 
                end
            end
            println("--------")
        end
    end
    return vctNDalloc
end


# ==============================================================================
# Labeling procedure for a given box
function labelingOneBox!(dataFull::Instance, b::Box, all_YN)

    verboseDev2 ? println("=================================================================") : nothing
    verboseDev2 ? println("J1 $(b.J1)") : nothing


    # -------------------------------------------------------------------------    
    # Compute the nondominated allocations for all users ---------------------- 

  #=  nbsingleAlloc = 0
    #if verboseDev2 
        for iUser in 1:data.nI
            if length(J1nd[iUser]) == 1
            println("i:$iUser  j:$(J1nd[iUser])")
            nbsingleAlloc+=1
            end
        end
    #end
    @show nbsingleAlloc
=#

    data = deepcopy(dataFull)

    #J1nd = computeListNondominatedAllocationJ1allUsers(data, b.J1); k1=0; k2=0
    J1nd, k1, k2 = computeLocalDominatedAllocation(data::Instance, b::Box)
    #@show J1nd
    #@show k1, k2, data.nI


    # The datastructure for representing the graphs of allocations ------------
    labeling = Array{Union{Nothing,Vector{Vector{Int64}}}}(nothing,data.nJ,2)

    # Labels candidates generated
    candidates = Array{Union{Nothing,Vector{Vector{Int64}}}}(nothing,data.nJ)


    # -------------------------------------------------------------------------
    # Compute the first layer of the graph ------------------------------------
    # As for the first user (i.e first layer) the dummy label [(0,0)[]] has to be
    # the first layer is obtained in computing (0,0) + (c1,c2) and []//service
    
    i=1
    verboseDev2 ? println("LAYER $i ==================================================") : nothing
    for (pos,j) in enumerate(J1nd[i])
        labeling[j,layer(i)] = [ [ data.c1[i,j] , data.c2[i,j] ] ]
        verboseDev2 ? println("$i $j [$(data.c1[i,j]) $(data.c2[i,j]) ; $j]") : nothing
    end

    # -------------------------------------------------------------------------    
    # Compute the others layers of the graph ----------------------------------
    # Each label of the previous layer is propagated to the current layer.
    # The obtained label is checked to detect if
    # - it is not dominated by an other label of the layer; if yes => pruned
    # - it is not dominating an other label of the layer; if yes => discard the dominated label
    # - it is not dominated by a bound know on the box ; if yes => pruned
    # - it is not dominated by a bound know on a other box ; if yes => pruned

    for i in 2:data.nI
        verboseDev2 ? println("LAYER $i ==================================================") : nothing

        # prepare the layer for the user i in filling the structure with nothing
        labeling[:,layer(i)] = fill(nothing,data.nJ)

        for (pos,jPrec) in enumerate( J1nd[ i-1 ] )

            candidates[:] = fill(nothing,data.nJ)
              
            # forall nondominated services opened on the preceding users ------
            verboseDev2 ? println("  i(previous)=$(i-1) j=$jPrec ") : nothing


            if labeling[jPrec,layer(i-1)] != nothing

                for (pos,label) in enumerate( labeling[jPrec,layer(i-1)] )

                    verboseDev2 ? println("      label to propagate : ", label) : nothing
                    # forall existing labels --------------------------------------

                    for (pos,j) in enumerate(J1nd[i])
                        # compute the labels  user i and service j∈J1nd -----------

                        Cnew = [ label[1] + data.c1[i,j] , label[2] + data.c2[i,j] ]
                        verboseDev2 ? println("      >>>  allocation cost of i=$i at j=$j : ($(data.c1[i,j]) $(data.c2[i,j]))$j  --> giving [$(Cnew[1]) $(Cnew[2])]") : nothing
 
                        # test the propagated cost label.CA+C -------------------------

                        # ---------------------------------------------------------
                        # 1. test with a local bound (known on the current box) 
                        if  (    (b.CR[2] + Cnew[2]+k2 > b.y12[2]) # out of the box at NW
                              || (b.CR[1] + Cnew[1]+k1 > b.y21[1]) # out of the box at SE
                            )
                            # new label is out of the bounds of its reduced box => discarded
                            verboseDev2 ? println("test 1 true: OUT OF BOX") : nothing

                        else
                            if candidates[j] != nothing
                                # labels already  available -> add the label to the list of label(s)
                                push!(candidates[j],Cnew)
                            else
                                # label not yet present -> initiate the list with the label
                                candidates[j] = [Cnew]
                            end
                        end
                    end
                end
            end
            #@show candidates 

            # ---------------------------------------------------------
            # 3. tests with labels already available on the current layer            
            for (pos,jLayerCur) in enumerate(J1nd[i])
                if candidates[jLayerCur] != nothing
                    nbAdded = zeros(Int64, data.nJ)
                    for (pos,cand) in enumerate(candidates[jLayerCur])

                        addNewLabel = true

                        # test with all others labels already on the layer
                        for (pos,jND) in enumerate(J1nd[i])

                            if labeling[jND,layer(i)] != nothing

                                #nbdeleted = 0
                                long = length(labeling[jND,layer(i)])
                                dern = long - nbAdded[jND]

                                for iLabel in reverse(1:dern)
                                    if isDominates(labeling[jND,layer(i)][iLabel],cand)
                                        #println(" 1) Oui en $iLabel. $(labeling[jND,layer(i)][iLabel]) $(cand)")
                                        addNewLabel = false
                                        break # again a dirty way to leave the loop
                                    else
                                        if isDominates(cand,labeling[jND,layer(i)][iLabel])
                                            #println(" 2) Oui en $iLabel. $(labeling[jND,layer(i)][iLabel]) $(cand)")
                                            deleteat!(labeling[jND,layer(i)],iLabel)
                                        #    labeling[jND,layer(i)][iLabel] = labeling[jND,layer(i)][long-nbdeleted]
                                        #    nbdeleted+=1
                                        end
                                    end
                                end
                                #resize!(labeling[jND,layer(i)],long-nbdeleted)
                              
                            end
                        end # test with all others labels already on the layer

                        if addNewLabel               
                            if labeling[jLayerCur,layer(i)] != nothing
                                # labels already  present -> add the label to the list of label(s)
                                push!(labeling[jLayerCur,layer(i)],cand)
                            else
                                # label not yet present -> initiate the list with the label
                                labeling[jLayerCur,layer(i)] = [cand]
                            end
                            nbAdded[jLayerCur]+=1
                            verboseDev2 ? println("           ADD : ",cand) : nothing
                        else 
                            verboseDev2 ? println("           DISCARD (new label is dominated) : ",cand) : nothing
                        end 
                    end
                    #@show nbAdded
                end
            end

        end # nondominated services opened on the previous layer


    end # users         


    # -----------------------------------------------------------------------------
    # Extract all the labels available on the last layer --------------------------
    for (pos,jND) in enumerate(J1nd[data.nI])
        # jND : services nondominated opened on the last layer

        if labeling[jND,layer(data.nI)] != nothing

            for ilabel in 1:length(labeling[jND,layer(data.nI)])
                verboseDev2 ? println(labeling[jND,layer(data.nI)][ilabel]) : nothing
                verboseDev2 ? println(labeling[jND,layer(data.nI)][ilabel]+b.CR) : nothing
                push!(all_YN, (b.CR[1]+labeling[jND,layer(data.nI)][ilabel][1]+k1, b.CR[2]+labeling[jND,layer(data.nI)][ilabel][2]+k2) )

                #println(J1," ", labeling2[layer(data.nI)][jND][ilabel].CA+b.CR, labeling2[layer(data.nI)][jND][ilabel])
                # FAIRE : stocker dans la box (champ a ajouter) pour archiver CA et allocation correspondante
            end
        end

    end

    #@show all_YN
    #@assert false "stop"
    return nothing
end


# ==============================================================================
function labelingPaving!(data::Instance, paving::Vector{Box})

    #@show data.c1  # allocation users-services cost 1 (matrix nI x nJ)
    #@show data.c2  # allocation users-services cost 2 (matrix nI x nJ)

    all_YN=(Tuple{Int64, Int64})[]
    verboseDev2 ? println("Number of boxes to develop: ", length(paving)) : nothing
    numbox=0

    # Develop all the boxes in the paving -------------------------------------
    for ib in eachindex(paving)
        numbox+=1
        #verboseDev2 ? 
        println("  Develop box number: $numbox", paving[ib].J1) #: nothing
        if length(paving[ib].J1) > 1
            labelingOneBox!(data, paving[ib], all_YN)
        else
            # box reduced to a point => added to all_YN and nothing to do
            push!(all_YN, (paving[ib].y12[1], paving[ib].y12[2]) )
        end
    end # 

    verboseDev2 ? println("all_YN", all_YN) : nothing

    # Extract the global set of nondominated points from the local sets of nondominated points (from all boxes) 
    ND_YN = getNonDominatedPoints(all_YN)
    #verboseDev2 ? 
    #println("ND_YN", ND_YN) #: nothing

    return ND_YN, all_YN
end
