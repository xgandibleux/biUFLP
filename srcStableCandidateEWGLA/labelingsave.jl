const withTraceAllocation = false

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

        #old
        #J1nd = computeListNondominatedAllocationJ1oneUser(data, iUser, J1)
        #vctNDalloc[iUser] = copy(J1nd)
        #new
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
function labelingOneBox!(data::Instance, b::Box, all_YN)

    verboseDev2 ? println("=================================================================") : nothing
    verboseDev2 ? println("J1 $(b.J1)") : nothing


    # -------------------------------------------------------------------------    
    # Compute the nondominated allocations for all users ---------------------- 
    J1nd = computeListNondominatedAllocationJ1allUsers(data, b.J1)

    if verboseDev2 
        for iUser in 1:data.nI
            println("i:$iUser  j:$(J1nd[iUser])")
        end
    end


    # The datastructure for representing the graphs of allocations ------------
    # services-users is a vector of dictionaries where
    # - keys are nondominated allocations services-users for a given level
    # - values are a list of labels
    labeling2=Vector{Dict{Int, Vector{Label}}}(undef,data.nI)


    # -------------------------------------------------------------------------
    # Compute the first layer of the graph ------------------------------------
    # As for the first user (i.e first layer) the dummy label [(0,0)[]] has to be
    # the first layer is obtained in computing (0,0) + (c1,c2) and []//service
    i=1
    verboseDev2 ? println("LAYER $i ==================================================") : nothing
    
    # Create one empty dictionary for user i
    labeling2[layer(i)] = Dict{Int,Vector{Label}}()

    for (pos,j) in enumerate(J1nd[i])
        # compute the labels for user i - service j∈J1nd
        if withTraceAllocation
            @inbounds labeling2[layer(i)][j] = [Label([ data.c1[i,j] , data.c2[i,j] ], Vector{UInt16}(undef,1))]
            @inbounds labeling2[layer(i)][j][1].listJ1[1] = UInt16(j)
        else
            @inbounds labeling2[layer(i)][j] = [Label([ data.c1[i,j] , data.c2[i,j] ], [])]
        end

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
        #print(i," ")

        # Create one empty dictionary for user i 
        labeling2[layer(i)] = Dict{Int,Vector{Label}}()

        for jPrec in keys(labeling2[layer(i-1)])
            # forall nondominated services opened on the preceding users ------
            verboseDev2 ? println("  i(previous)=$(i-1) j=$jPrec ") : nothing

            for label in values(labeling2[layer(i-1)][jPrec])
                # forall existing labels --------------------------------------

                verboseDev2 ? println("      label to propagate : ", label.CA, " ", label.listJ1) : nothing

                for (pos,j) in enumerate(J1nd[i])
                    # compute the labels  user i and service j∈J1nd -----------

                    @inbounds Cnew = [ label.CA[1] + data.c1[i,j] , label.CA[2] + data.c2[i,j] ]
                    verboseDev2 ? println("      >>>  allocation cost of i=$i at j=$j : ($(data.c1[i,j]) $(data.c2[i,j]))$j  --> giving [$(Cnew[1]) $(Cnew[2])]") : nothing

                    # test the propagated cost label.CA+C -------------------------

                    addNewLabel = true

                    # ---------------------------------------------------------
                    # 1. test with a local bound (known on the current box) 
                    if  (    (b.CR[2] + Cnew[2] > b.y12[2]) # out of the box at NW
                          || (b.CR[1] + Cnew[1] > b.y21[1]) # out of the box at SE
                        )
                        # new label is out of the bounds of its reduced box => discarded
                        verboseDev2 ? println("test 1 true: OUT OF BOX") : nothing
                        addNewLabel = false
                       end

                    # ---------------------------------------------------------                       
                    # 2. test with a global bound (known on the others boxes) -
                        # todo

                    # ---------------------------------------------------------
                    # 3. tests with labels already available on the current layer
                    if addNewLabel

                        for jND in keys(labeling2[layer(i)])
                            # j : services nondominated opened on the current layer

                           !addNewLabel && break # dirty way to leave the 'for jND loop' when the new label is dominated

                            # test : verify if it exists one label of 'labeling2[layer(i)][jND]' who dominates 'Cnew'
                            # if yes, the new label is pruned
                            verboseDev2 ? println("               listCA :", labeling2[layer(i)][jND]) : nothing
                            verboseDev2 ? println("               TEST 1 $jND (isDominates(x.CA,Cnew)?) : ",any(x -> isDominates(x.CA,Cnew), labeling2[layer(i)][jND]) ) : nothing

                            #println(length(labeling2[layer(i)][jND]))
                            # REV1 : begin
                            # old
                            addNewLabel = !any(x -> isDominates(x.CA,Cnew), labeling2[layer(i)][jND])
                            # new
                            #for iLabel in eachindex(labeling2[layer(i)][jND])
                            #    if isDominates(labeling2[layer(i)][jND][iLabel].CA,Cnew)
                            #        #println(" Oui en $i. $(labeling2[layer(i)][jND][iLabel].CA) $(Cnew)")
                            #        addNewLabel = false
                            #        break # again a dirty way to leave the loop
                            #    end
                            #end
                            # REV1 : end


                            # REV1 : begin
                            # old                            
                            if addNewLabel != false
                                # test : remove all labels of 'labeling2[layer(i)][jND]' dominated by 'Cnew'                      
                                verboseDev2 ? println("               TEST 2 $jND (isDominates(Cnew,x.CA)?) : ",any(x -> isDominates(Cnew,x.CA,), labeling2[layer(i)][jND])) : nothing
                                filter!(x -> !isDominates(Cnew,x.CA), labeling2[layer(i)][jND])
                            end
                            #new
                            #if addNewLabel
                            #    for iLabel in reverse(eachindex(labeling2[layer(i)][jND]))
                            #        if isDominates(Cnew,labeling2[layer(i)][jND][iLabel].CA)
                            #          deleteat!(labeling2[layer(i)][jND],iLabel)
                            #        end
                            #    end                         
                            #end
                            # REV1 : end  
                            

                            # REV1 : begin  
                            # old
                            # any et filter v deleteat   
                            #for iLabel in reverse(eachindex(labeling2[layer(i)][jND]))
                            #    if isDominates(labeling2[layer(i)][jND][iLabel].CA,Cnew)
                            #        #println(" Oui en $i. $(labeling2[layer(i)][jND][iLabel].CA) $(Cnew)")
                            #        addNewLabel = false
                            #        break # again a dirty way to leave the loop
                            #    else
                            #        if isDominates(Cnew,labeling2[layer(i)][jND][iLabel].CA)
                            #            deleteat!(labeling2[layer(i)][jND],iLabel)
                            #        end                                    
                            #    end
                            #end                       
                            # REV1 : end  


                            # REV1 : begin  
                            # old
                            # any et filter  v copy
                            #nbpruned = 0
                            #lg = length(labeling2[layer(i)][jND])
                            #iLabel = 1  
                            #while iLabel <= lg-nbpruned
                            #    if isDominates(labeling2[layer(i)][jND][iLabel].CA,Cnew)
                            #        #println(" Oui en $i. $(labeling2[layer(i)][jND][iLabel].CA) $(Cnew)")
                            #        addNewLabel = false
                            #        break # again a dirty way to leave the loop
                            #    else
                            #        if isDominates(Cnew,labeling2[layer(i)][jND][iLabel].CA)
                            #            # copy the reference only
                            #            labeling2[layer(i)][jND][iLabel] = labeling2[layer(i)][jND][end-nbpruned]
                            #            nbpruned+=1
                            #        else
                            #            iLabel+=1
                            #        end                                    
                            #    end
                            #end 
                            #if nbpruned > 0
                            #    resize!(labeling2[layer(i)][jND],lg-nbpruned)
                            #end                      
                            # REV1 : end  

                            verboseDev2 ? println(" ") : nothing              
                        end

                    end # test 3

                    # ---------------------------------------------------------                    
                    # if lnew is nondominated it is added in the list of labels 
                    if addNewLabel               
                        if withTraceAllocation 
                            #lnew = Label(Cnew, copy(label.listJ1))
                            lnew = Label(Cnew, Vector{UInt16}(undef,i))
                            lnew.listJ1[1:i-1] = label.listJ1[1:i-1]
                            @inbounds lnew.listJ1[i] = UInt16(j)     # allocation user i to service j
                        else 
                            lnew = Label(Cnew, [])  
                        end                             
                        verboseDev2 ? println("           ADD : ",lnew) : nothing

                        if j in keys(labeling2[layer(i)])
                            # key  present -> add the label to the list of label(s)
                            push!(labeling2[layer(i)][j],lnew)
                        else
                            # key not yet present -> initiate the list with the label
                            labeling2[layer(i)][j] = [lnew]
                        end
                    else 
                        verboseDev2 ? println("           DISCARD (new label is dominated) : ",Cnew) : nothing
                        #@assert false "stop"
                    end #addNewLabel

                end # nondominated services opened on the current layer

            end # labels of nondominated services j on previous layer

            verboseDev2 ? println(" ") : nothing
        end # nondominated services opened on the previous layer

    end # users


    # -----------------------------------------------------------------------------
    # Extract all the labels available on the last layer --------------------------
    for jND in keys(labeling2[layer(data.nI)])
        # jND : services nondominated opened on the last layer

        for ilabel in 1:length(labeling2[layer(data.nI)][jND])
            verboseDev2 ? println(labeling2[layer(data.nI)][jND][ilabel]) : nothing
            verboseDev2 ? println(labeling2[layer(data.nI)][jND][ilabel].CA+b.CR) : nothing
            push!(all_YN, (b.CR[1]+labeling2[layer(data.nI)][jND][ilabel].CA[1], b.CR[2]+labeling2[layer(data.nI)][jND][ilabel].CA[2]) )

            #println(J1," ", labeling2[layer(data.nI)][jND][ilabel].CA+b.CR, labeling2[layer(data.nI)][jND][ilabel])
            # FAIRE : stocker dans la box (champ a ajouter) pour archiver CA et allocation correspondante
        end

    end

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
        #println("  Develop box number: $numbox", paving[ib].J1) #: nothing
        labelingOneBox!(data, paving[ib], all_YN)
    end # 

    verboseDev2 ? println("all_YN", all_YN) : nothing

    # Extract the global set of nondominated points from the local sets of nondominated points (from all boxes) 
    ND_YN = getNonDominatedPoints(all_YN)
    verboseDev2 ? println("ND_YN", ND_YN) : nothing
    #println("ND_YN", ND_YN)

    return ND_YN, all_YN
end
