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

#global cpt1=0
global nbmaxlabel = 0
global inbmaxlabel = 1

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

    for j in J1nd[i]
        push!( labeling[layer(i)] , Label( [ b.CR[1] + k1 + data.c1[i,j] , b.CR[2] + k2 + data.c2[i,j] ], [ j ] ) )
        verboseDev2 ? println("$i $j [$(data.c1[i,j]) $(data.c2[i,j]) ; $j]") : nothing
    end

    global nbmaxlabel = max(nbmaxlabel,length(labeling[layer(i)]))


    # -------------------------------------------------------------------------
    # Compute the others layers of the graph ----------------------------------
    # Each label of the previous layer is propagated to the current layer.
    # The obtained label is checked to detect if
    # - it is not dominated by a bound know on the box ; if yes => pruned
    # - it is not dominated by an other label of the layer; if yes => pruned
    # - it is not dominating an other label of the layer; if yes => discard the dominated label

    #candidate = Label( [ 0 , 0 ], [ 0 ] )
    #yCandidate :: Vector{Int64} = [ 0 , 0 ]
    indices::Vector{Int64} = Vector{Int64}(undef,0) # ligne bidon pour faire une déclaration des indices service->client
    Tas::MutableBinaryMinHeap{zInd} = MutableBinaryMinHeap{zInd}() # Déclaration de la variable Tas
    somme::Int64 = 0 # Déclaration de la somme des indices (utilisée pour la condition d'arrêt de l'énumération)

    @time for i in 2:data.nI
        verboseDev2 ? println("LAYER $i ==================================================") : nothing
        indices = ones(Int64,length(J1nd[i])) # Tous les indices sont initialement à 1
        #lprec = labeling[layer(i-1)] # Alias pour éviter une indirection
        #=if (i == 2)
            println("2 : $lprec")
            println(J1nd[i])
            for j in J1nd[i]
                println("($(data.c1[i,j]),$(data.c2[i,j]))")
            end
        end=#
        l = labeling[layer(i-1)][1] # On part du premier label de la colonne précédente (puisqu'il est le plus petit lexicographiquement)
        Tas = MutableBinaryMinHeap{zInd}([zInd(l.CA[1] + data.c1[i,j],l.CA[2] + data.c2[i,j],pos) for (pos,j) in enumerate(J1nd[i])]) # Initialisation du Tas
        #=if (i == 2)
            println("Tas initial : $Tas")
        end=#

        # initialize the list of nondominated labels for the user i
        labeling[layer(i)] = (Vector{Label})[]
        #lsuiv = labeling[layer(i)] # Alias pour éviter une indirection

        # Lancement de la boucle principale
        somme = length(J1nd[i])
        fin = length(J1nd[i]) * (length(labeling[layer(i-1)]) + 1)

        while (somme < fin)
            # On prend toujours l'élément en tête pour instancier le nouveu label (qui sera le plus petit lexicographiquement des labels restant à générer)
            ztete, handle = top_with_handle(Tas)
            #if (i == 2) println("ztete = $ztete") end
            # ---------------------------------------------------------
            # 1. test with a local bound (known on the current box)
            if ( (ztete.z1 > b.y21[1]) || (ztete.z2 > b.y12[2]) )
                # new label is out of the bounds of its reduced box => discarded
                verboseDev2 ? println("test 1 true: OUT OF BOX") : nothing
                #if (i == 2) println("suppression immédiate") end
            else
                # Test de dominance du nouveau label avant insertion éventuelle
                if ( (length(labeling[layer(i)]) == 0) || (ztete.z2 < labeling[layer(i)][end].CA[2]) )
                    push!(labeling[layer(i)],Label([ztete.z1,ztete.z2],[UInt16(0)])) # Je ne comprends pas vraiment ce dernier champ dans Label, il est complété uniquement pour la compilation
                    #if (i == 2) println("MAJ lsuiv : $lsuiv") end
                end
            end

            indices[ztete.ind] += 1
            somme += 1
            #global cpt1+=1
            if (indices[ztete.ind] <= length(labeling[layer(i-1)]))
                #l = labeling[layer(i-1)][indices[ztete.ind]]
                update!(Tas, handle, zInd(labeling[layer(i-1)][indices[ztete.ind]].CA[1] + data.c1[i,J1nd[i][ztete.ind]], labeling[layer(i-1)][indices[ztete.ind]].CA[2] + data.c2[i,J1nd[i][ztete.ind]],ztete.ind))
            else
                pop!(Tas)
            end
            #if (i == 2) println("MAJ Tas : $Tas") end
        end
        if nbmaxlabel < length(labeling[layer(i)])
          global nbmaxlabel = length(labeling[layer(i)])
          global inbmaxlabel = i
        end
    end


    # -----------------------------------------------------------------------------
    # Extract all the labels available on the last layer --------------------------

    #@show labeling[layer(data.nI)]
    #@assert false "stop" 
    for label in labeling[layer(data.nI)]
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
    @show nbmaxlabel
    @show inbmaxlabel

    return ND_YN, all_YN
end
