# ==============================================================================
# fname : labelingBourrin.jl
# September 2022
# ==============================================================================

using DataStructures
import Base.isless

mutable struct zInd
    z1::Int64
    z2::Int64
    ind::Int64
end

isless(a::zInd,b::zInd) = 10000*a.z1+ a.z2 < 10000*b.z1 + b.z2

#=
"""
    updatexg!{T}(h::MutableBinaryHeap{T}, i::Int, v::T)
Replace the element at index `i` in heap `h` with `v`.
"""
using Base.Order: Forward, Ordering, lt
function updatexg!(h::MutableBinaryHeap{T}, i::Int, z1,z2,ind) where T
    nodes = h.nodes
    nodemap = h.node_map
    ordering = h.ordering

    nd_id = nodemap[i]
    v0 = nodes[nd_id].value
    #x = convert(T, v)
    #@show nodes[nd_id].value.z1
    nodes[nd_id].value.z1 = z1
    nodes[nd_id].value.z2 = z2    
    nodes[nd_id].value.ind = ind
    # = MutableBinaryHeapNode(x, i)
    #if Base.lt(ordering, x, v0)
    if 10000*z1+ z2 < 10000*v0.z1 + v0.z2
        _heap_bubble_up!(ordering, nodes, nodemap, nd_id)
    else
        _heap_bubble_down!(ordering, nodes, nodemap, nd_id)
    end
end
struct MutableBinaryHeapNode{T}
    value::T
    handle::Int
end
function _heap_bubble_up!(ord::Ordering,
    nodes::Vector{MutableBinaryHeapNode{T}}, nodemap::Vector{Int}, nd_id::Int) where T

    @inbounds nd = nodes[nd_id]
    v::T = nd.value

    swapped = true  # whether swap happens at last step
    i = nd_id

    while swapped && i > 1  # nd is not root
        p = i >> 1
        @inbounds nd_p = nodes[p]

        if Base.lt(ord, v, nd_p.value)
            # move parent downward
            @inbounds nodes[i] = nd_p
            @inbounds nodemap[nd_p.handle] = i
            i = p
        else
            swapped = false
        end
    end

    if i != nd_id
        nodes[i] = nd
        nodemap[nd.handle] = i
    end
end

function _heap_bubble_down!(ord::Ordering,
    nodes::Vector{MutableBinaryHeapNode{T}}, nodemap::Vector{Int}, nd_id::Int) where T

    @inbounds nd = nodes[nd_id]
    v::T = nd.value

    n = length(nodes)
    last_parent = n >> 1

    swapped = true
    i = nd_id

    while swapped && i <= last_parent
        il = i << 1

        if il < n   # contains both left and right children
            ir = il + 1

            # determine the better child
            @inbounds nd_l = nodes[il]
            @inbounds nd_r = nodes[ir]

            if Base.lt(ord, nd_r.value, nd_l.value)
                # consider right child
                if Base.lt(ord, nd_r.value, v)
                    @inbounds nodes[i] = nd_r
                    @inbounds nodemap[nd_r.handle] = i
                    i = ir
                else
                    swapped = false
                end
            else
                # consider left child
                if Base.lt(ord, nd_l.value, v)
                    @inbounds nodes[i] = nd_l
                    @inbounds nodemap[nd_l.handle] = i
                    i = il
                else
                    swapped = false
                end
            end

        else  # contains only left child
            nd_l = nodes[il]
            if Base.lt(ord, nd_l.value, v)
                @inbounds nodes[i] = nd_l
                @inbounds nodemap[nd_l.handle] = i
                i = il
            else
                swapped = false
            end
        end
    end

    if i != nd_id
        @inbounds nodes[i] = nd
        @inbounds nodemap[nd.handle] = i
    end
end

=#

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


# ==============================================================================
# Labeling procedure for a given box
function labelingOneBox!(dataFull::Instance, b::Box, labeling::Array{Vector{Label}}, iVctLabels::Vector{Int64}, all_YN)

    verboseDev2 ? println("=================================================================") : nothing
    verboseDev2 ? println("J1 $(b.J1)") : nothing

    # -------------------------------------------------------------------------
    # Compute the nondominated allocations for all users ----------------------

    data = deepcopy(dataFull)

    #J1nd = computeListNondominatedAllocationJ1allUsers(data, b.J1); k1=0; k2=0
    J1nd, k1, k2 = computeLocalDominatedAllocation!(data::Instance, b::Box)


    # -------------------------------------------------------------------------
    # Compute the first layer of the graph ------------------------------------

    i=1
    # initialize the list of nondominated labels for the user i
    iVctLabels[layer(i)] = 0

    for j in J1nd[i]

        iVctLabels[layer(i)]+=1
        labeling[layer(i)][iVctLabels[layer(i)]].CA[1] = b.CR[1] + k1 + data.c1[i,j]
        labeling[layer(i)][iVctLabels[layer(i)]].CA[2] = b.CR[2] + k2 + data.c2[i,j]
        labeling[layer(i)][iVctLabels[layer(i)]].listJ1[1] = UInt16(0)   

    end


    # -------------------------------------------------------------------------
    # Compute the others layers of the graph ----------------------------------
    # Each label of the previous layer is propagated to the current layer.
    # The obtained label is checked to detect if
    # - it is not dominated by a bound know on the box ; if yes => pruned
    # - it is not dominated by an other label of the layer; if yes => pruned
    # - it is not dominating an other label of the layer; if yes => discard the dominated label

    indices::Vector{Int64} = Vector{Int64}(undef,data.nJ) # ligne bidon pour faire une déclaration des indices service->client
    Tas::MutableBinaryMinHeap{zInd} = MutableBinaryMinHeap{zInd}() # Déclaration de la variable Tas
    somme::Int64  = 0 # Déclaration de la somme des indices (utilisée pour la condition d'arrêt de l'énumération)
    ztete::zInd   = zInd(0,0,0)
    handle::Int64 = 0
    zmaj::zInd = zInd(0,0,0)

     for i in 2:data.nI
        #verboseDev2 ? println("LAYER $i ==================================================") : nothing
        indices = ones(Int64,length(J1nd[i])) # Tous les indices sont initialement à 1

        Tas = MutableBinaryMinHeap{zInd}([ zInd( labeling[layer(i-1)][1].CA[1] + data.c1[i,j] , 
                                                 labeling[layer(i-1)][1].CA[2] + data.c2[i,j] , 
                                                 pos
                                               ) for (pos,j) in enumerate(J1nd[i])
                                         ]
                                        ) # Initialisation du Tas

        # initialize the list of nondominated labels for the user i
        iVctLabels[layer(i)] = 0

        # Lancement de la boucle principale
        somme = length(J1nd[i])
        fin = length(J1nd[i]) * (iVctLabels[layer(i-1)] + 1)

        while (somme < fin)
            # On prend toujours l'élément en tête pour instancier le nouveu label (qui sera le plus petit lexicographiquement des labels restant à générer)
            ztete, handle = top_with_handle(Tas)

            # ---------------------------------------------------------
            # 1. test with a local bound (known on the current box)
            if ( (ztete.z1 > b.y21[1]) || (ztete.z2 > b.y12[2]) )
                # new label is out of the bounds of its reduced box => discarded
                verboseDev2 ? println("test 1 true: OUT OF BOX") : nothing
            else
                # Test de dominance du nouveau label avant insertion éventuelle
                if ( (iVctLabels[layer(i)] == 0) || (ztete.z2 < labeling[layer(i)][iVctLabels[layer(i)]].CA[2]) )
 
                    iVctLabels[layer(i)]+=1
                    labeling[layer(i)][iVctLabels[layer(i)]].CA[1] = ztete.z1
                    labeling[layer(i)][iVctLabels[layer(i)]].CA[2] = ztete.z2
                    labeling[layer(i)][iVctLabels[layer(i)]].listJ1[1] = UInt16(0)         

                end
            end

            indices[ztete.ind] += 1
            somme += 1
            if (indices[ztete.ind] <= iVctLabels[layer(i-1)])

            #    updatexg!(Tas, handle,  labeling[layer(i-1)][indices[ztete.ind]].CA[1] + data.c1[i,J1nd[i][ztete.ind]], 
            #                            labeling[layer(i-1)][indices[ztete.ind]].CA[2] + data.c2[i,J1nd[i][ztete.ind]] , 
            #                            ztete.ind)
                
                update!(Tas, handle, zInd( labeling[layer(i-1)][indices[ztete.ind]].CA[1] + data.c1[i,J1nd[i][ztete.ind]], 
                                           labeling[layer(i-1)][indices[ztete.ind]].CA[2] + data.c2[i,J1nd[i][ztete.ind]] , 
                                           ztete.ind
                                         )
                       )
         #              @show Tas
            else
                pop!(Tas)
            end
        end  

    end


    # -----------------------------------------------------------------------------
    # Extract all the labels available on the last layer --------------------------

    for ilabel in 1:iVctLabels[layer(data.nI)]
        push!(all_YN, ( labeling[layer(data.nI)][ilabel].CA[1] , labeling[layer(data.nI)][ilabel].CA[2] ) )
    end

    return nothing
end


# ==============================================================================
function labelingPaving!(data::Instance, paving::Vector{Box})

    all_YN = (Tuple{Int64, Int64})[] # temporary solution to collect the ND points
    verboseDev2 ? println("Number of boxes to develop: ", length(paving)) : nothing
    numbox=0

    # The datastructure for representing the non-dominated labels in the graph of allocation ------------
    labeling = Array{Vector{Label}}(undef,2)
    labeling[1] = Vector{Label}(undef,250000)
    labeling[2] = Vector{Label}(undef,250000)   
    println("allocation a priori")
    @time for i in 1:250000
        labeling[1][i]=Label([typemax(Int64),typemax(Int64)],[typemax(UInt16)]) 
        labeling[2][i]=Label([typemax(Int64),typemax(Int64)],[typemax(UInt16)])      
    end  
    iVctLabels = [0,0]

    # Develop all the boxes in the paving -------------------------------------
    for ib in eachindex(paving)
        numbox+=1
        #verboseDev2 ?
        println("  Develop box number: $numbox", paving[ib].J1) #: nothing
        if length(paving[ib].J1) > 1
            labelingOneBox!(data, paving[ib], labeling, iVctLabels, all_YN)
        else
            # box reduced to a point => added to all_YN and nothing to do
            push!(all_YN, (paving[ib].y12[1], paving[ib].y12[2]) )
        end
    end 

    verboseDev2 ? println("all_YN", all_YN) : nothing

    # Extract the global set of nondominated points from the local sets of nondominated points (from all boxes)
    ND_YN = getNonDominatedPoints(all_YN)
    #verboseDev2 ? println("ND_YN", ND_YN) : nothing

    return ND_YN, all_YN
end
