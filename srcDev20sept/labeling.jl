# ==============================================================================
# fname : labeling.jl
# September 2022
# ==============================================================================

# Version with only 2 layers and one list of ND labels per layer

const withTraceAllocation = false


# ==============================================================================
"""
    layer(i::Int64)

    Switch between the two layers depending of the value of 'i'
"""
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
"""
    filter2!(NDlabels::Vector{Label}, yCandidate0::Vector{Int64}, j::Int64)

    Given a candidate label, maintain the list of nondominated labels on a layer
"""

#cpt1 = 0
#cpt2 = 0
#cpt3 = 0
#cpt4 = 0
#cpt5 = 0
#cpt6 = 0
#nbdeleted = 0

function filter2!(NDlabels::Vector{Label}, yCandidate0::Vector{Int64}, j::Int64)

    nbrNDlabels :: Int64 = length(NDlabels)

    if nbrNDlabels == 0
        # the list is empty -> insert -----------------------------------------
        # F50-51: 576 times
        # H2000 : 514 times
        push!(NDlabels,Label( copy(yCandidate0),[UInt16(j)] ))
        #global cpt1+=1

    else
        # search the position for insering the candidate label ----------------
        # F50-51:    18 925 300 times
        # H2000 : 1 591 731 696 times
        #ipos = 1
        #while (ipos <= nbrNDlabels) && ( @inbounds NDlabels[ipos].CA[1] < @inbounds yCandidate0[1] )
        #    ipos+=1
            #global cpt2+=1
        #end
        ipos = searchsortedfirst(NDlabels, yCandidate0[1], lt=(r,x)->r.CA[1]<x)
        
        # check for a possible insertion in first position --------------------
        # F50-51: 88 times
        # H2000 :  0 times
        if (ipos == 1)  && (! isDominates(NDlabels[ipos].CA,yCandidate0) )
            pushfirst!(NDlabels,Label(copy(yCandidate0),[UInt16(j)]))
            ipos+=1
            #global cpt3+=1

        # check for a possible insertion in middle ----------------------------
        # F50-51:  49 355 times
        # H2000 : 268 404 times
        elseif (ipos <= nbrNDlabels) && ! isDominates(NDlabels[ipos].CA,yCandidate0) && (  @inbounds NDlabels[ipos-1].CA[2] >   @inbounds yCandidate0[2] )
            splice!(NDlabels,ipos:ipos-1,[Label(copy(yCandidate0),[UInt16(j)])])
            ipos+=1	
            #global cpt4+=1	

        # check for a possible insertion in last position ---------------------
        # F50-51:  81 044 times
        # H2000 : 523 002 times
        elseif (ipos > nbrNDlabels)  &&  ! isDominates(NDlabels[ipos-1].CA,yCandidate0) 
            push!(NDlabels,Label(copy(yCandidate0),[UInt16(j)]))
            ipos+=1		
            #global cpt5+=1
        end	
        
        # remove from the list the dominated labels if they exist -------------
        # F50-51:  44 086 times
        # H2000 : 257 461 times

        #searchsortedfirst(NDlabels, yCandidate0[2], by=u->u.CA[2], lt=(u,v)->u<v)

        #filter!(x -> x.CA[2] >  yCandidate0[2], NDlabels[ipos:end])

        #=
        while (ipos <= nbrNDlabels) && (  @inbounds NDlabels[ipos].CA[2] >   @inbounds yCandidate0[2] )
            deleteat!(NDlabels,ipos)
            nbrNDlabels-=1
            global cpt6+=1
        end
        =#
        
        iposStart = ipos
        while (ipos <= nbrNDlabels) && (  @inbounds NDlabels[ipos].CA[2] >   @inbounds yCandidate0[2] )
            ipos+=1
        end
        if iposStart != ipos  
            ipos-=1
            deleteat!(NDlabels,iposStart:ipos)
        end        
        
    end
    return nothing
end


function addLabel!(NDlabels::Vector{Label}, yCandidate0::Vector{Int64}, j::Int64)

    if length(NDlabels) == 0
        # the list is empty -> insert -----------------------------------------
        push!(NDlabels,Label( copy(yCandidate0),[UInt16(j)] ))

    else
        # search the position for insering the candidate label ----------------
        ipos = searchsortedfirst(NDlabels, yCandidate0[1], lt=(r,x)->r.CA[1]<x)
        splice!(NDlabels,ipos:ipos-1,[Label(copy(yCandidate0),[UInt16(j)])])   
    end
    return nothing
end

#nbipos=0
#nbiposnext=0

function eraseDominatedlabels!(NDlabels::Vector{Label})

    nbrNDlabels :: Int64 = length(NDlabels)
    ipos        :: Int64 = 1
    iposNext    :: Int64 = 2

    while (iposNext <= nbrNDlabels)
        if @inbounds NDlabels[ipos].CA[1] == @inbounds NDlabels[iposNext].CA[1]
            if @inbounds NDlabels[ipos].CA[2] > @inbounds NDlabels[iposNext].CA[2]
                #NDlabels[setdiff(1:end, ipos)]
                deleteat!(NDlabels,ipos)
                #splice!(NDlabels,ipos)
                nbrNDlabels-=1
                #println("Delete ipos")
                #global nbipos+=1
            else
                #NDlabels[setdiff(1:end, iposNext)]
                deleteat!(NDlabels,iposNext)   
                #splice!(NDlabels,iposNext) 
                nbrNDlabels-=1
                #println("Delete iposnext")
                #global nbiposnext+=1
            end             
        else
            if @inbounds NDlabels[ipos].CA[2] > @inbounds NDlabels[iposNext].CA[2]
                ipos+=1
                iposNext+=1
            else
                #NDlabels[setdiff(1:end, iposNext)]
                deleteat!(NDlabels,iposNext)
                #splice!(NDlabels,iposNext)
                nbrNDlabels-=1
                #println("Delete iposnext")
                #global nbiposnext+=1
            end
        end
    end
    
    return nothing

end




# ==============================================================================
# Labeling procedure for a given box
function labelingOneBox!(dataFull::Instance, b::Box, all_YN)

    verboseDev2 ? println("=================================================================") : nothing
    verboseDev2 ? println("J1 $(b.J1)") : nothing

    # -------------------------------------------------------------------------    
    # Compute the nondominated allocations for all users ---------------------- 

    data = deepcopy(dataFull)

    #J1nd = computeListNondominatedAllocationJ1allUsers(data, b.J1); k1=0; k2=0
    J1nd, k1, k2 = computeLocalDominatedAllocation!(data::Instance, b::Box)

    # The datastructure for representing the non-dominated labels in the graph of allocation ------------
    labeling = Array{Vector{Label}}(undef,2)

    # -------------------------------------------------------------------------
    # Compute the first layer of the graph ------------------------------------
    
    i=1
    verboseDev2 ? println("LAYER $i ==================================================") : nothing

    # initialize the list of nondominated labels for the user i
    labeling[layer(i)] = (Vector{Label})[]

    for (pos,j) in enumerate(J1nd[i])
        push!( labeling[layer(i)] , Label( [ b.CR[1] + k1 + data.c1[i,j] , b.CR[2] + k2 + data.c2[i,j] ], [ j ] ) )
        verboseDev2 ? println("$i $j [$(data.c1[i,j]) $(data.c2[i,j]) ; $j]") : nothing
    end


    # -------------------------------------------------------------------------    
    # Compute the others layers of the graph ----------------------------------
    # Each label of the previous layer is propagated to the current layer.
    # The obtained label is checked to detect if
    # - it is not dominated by a bound know on the box ; if yes => pruned
    # - it is not dominated by an other label of the layer; if yes => pruned
    # - it is not dominating an other label of the layer; if yes => discard the dominated label

    #candidate = Label( [ 0 , 0 ], [ 0 ] )
    yCandidate :: Vector{Int64} = [ 0 , 0 ]

    for i in 2:data.nI
        verboseDev2 ? println("LAYER $i ==================================================") : nothing

        # initialize the list of nondominated labels for the user i 
        labeling[layer(i)] = (Vector{Label})[]

        # for all allocations j∈J1nd to consider for the user i
        for (pos,j) in enumerate(J1nd[i])
            #@show i,j

            # for all nondominated labels available on the previous layer 
            for (pos,l) in enumerate( labeling[layer(i-1)] )

                # propagate a label l
                yCandidate[1] = l.CA[1] + data.c1[i,j]     
                yCandidate[2] = l.CA[2] + data.c2[i,j] 
                #candidate.listJ1[1] =  j 

                # ---------------------------------------------------------
                # 1. test with a local bound (known on the current box) 
                if  (  (yCandidate[2] > b.y12[2]) # out of the box at NW
                    || (yCandidate[1] > b.y21[1]) # out of the box at SE
                    )
                    # new label is out of the bounds of its reduced box => discarded
                    verboseDev2 ? println("test 1 true: OUT OF BOX") : nothing
                else
                    # update the list of nondominated labels on layer i with the candidate label
                    addLabel!(labeling[layer(i)], yCandidate, j)
                    #@show yCandidate
                    #@show labeling[layer(i)]
                end              
            end     
        end
        #println(" ")
        eraseDominatedlabels!(labeling[layer(i)])
        #@assert false "stop" 
    end


    # -----------------------------------------------------------------------------
    # Extract all the labels available on the last layer --------------------------

    for (pos,label) in enumerate( labeling[layer(data.nI)] )
        push!(all_YN, ( label.CA[1] , label.CA[2] ) )
    end

    return nothing
end


# ==============================================================================
function labelingPaving!(data::Instance, paving::Vector{Box})

    #@show data.c1  # allocation users-services cost 1 (matrix nI x nJ)
    #@show data.c2  # allocation users-services cost 2 (matrix nI x nJ)

    all_YN = (Tuple{Int64, Int64})[] # temporary solution to collect the ND points
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
    #@show cpt1
    #@show cpt2
    #@show cpt3
    #@show cpt4
    #@show cpt5
    #@show cpt6     
    #@show nbipos
    #@show nbiposnext              

    return ND_YN, all_YN
end
