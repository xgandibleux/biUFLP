# ==============================================================================
# fname : labeling.jl
# September 2022
# ==============================================================================

const withTraceAllocation = false


# ==============================================================================
"""
    TasEchange!(tas::Vector{zInd}, i::Int64, j::Int64)

    Function exchanging two elements in a binary heap
"""
function TasEchange!(tas::Vector{zInd}, i::Int64, j::Int64)

    @inbounds  tas[i].value, tas[j].value =  tas[j].value, tas[i].value
    @inbounds  tas[i].z1,    tas[j].z1    =  tas[j].z1,    tas[i].z1
    @inbounds  tas[i].z2,    tas[j].z2    =  tas[j].z2,    tas[i].z2
    @inbounds  tas[i].ind,   tas[j].ind   =  tas[j].ind,   tas[i].ind

    return nothing
end


# ==============================================================================
"""
    TasDescend!(tas::Vector{zInd}, i::Int64, tailletas::Int64)

    Function used to put an element at its correct position in the binary heap
"""
function TasDescend!(tas::Vector{zInd}, i::Int64, tailletas::Int64)

    min::Int64 = 0
    if @inbounds ( (2 * i <= tailletas) && (tas[2 * i].value < tas[i].value) )
        min = 2 * i
    else
        min = i
    end

    if @inbounds ( (2 * i + 1 <= tailletas) && (tas[2 * i + 1].value < tas[min].value) )
        min = 2 * i + 1
    end

    if (min != i)
        TasEchange!(tas, i, min)
        TasDescend!(tas, min, tailletas)
    end
end


# ==============================================================================
"""
    layer(i::Int64)

    Switch between the two layers depending of the value of 'i'
"""
function layer(i::Int64)
    return mod(i+1,2)+1
end


# ==============================================================================
"""
    computeListNondominatedAllocationJ1oneUser(data::Instance, iUser::Int64, J1::Vector{Int64})

    Compute the list of nondominated allocation from J1 to one user
"""
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
"""
    computeListNondominatedAllocationJ1allUsers(data::Instance, J1::Vector{Int64})

    Compute the list of nondominated allocation from J1 to all users
"""
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
"""
    labelingOneBox!(dataFull::Instance, b::Box, labeling::Array{Vector{Label}}, all_YN)

    Labeling procedure for a given box
"""
function labelingOneBox!(dataFull::Instance, b::Box, labeling::Array{Vector{Label}}, maxLabels::Int64, all_YN)

    verboseDev2 ? println("=================================================================") : nothing
    verboseDev2 ? println("J1 $(b.J1)") : nothing

    # -------------------------------------------------------------------------
    # Compute the nondominated allocations for all users ----------------------

    nIfull = dataFull.nI
    data = deepcopy(dataFull)
    #J1nd = computeListNondominatedAllocationJ1allUsers(data, b.J1); k1=0; k2=0
    J1nd, k1, k2, dominantallocations = computeLocalDominatedAllocation!(data::Instance, b::Box)
    #@show allocations


    # -------------------------------------------------------------------------
    # Compute the first layer of the graph ------------------------------------

    # Initialize the list of nondominated labels for the user 1 

    ilabels::Int64 = 0
    for j in J1nd[1]

        ilabels +=1
        labeling[layer(1)][ilabels].CA[1]     =  b.CR[1] + k1 + data.c1[1,j]
        labeling[layer(1)][ilabels].CA[2]     =  b.CR[2] + k2 + data.c2[1,j]
        labeling[layer(1)][ilabels].listJ1[1] =  UInt16(j)   

    end    
    @show labeling[layer(1)]

    # -------------------------------------------------------------------------
    # Compute the others layers of the graph 
    # Each label of the previous layer is propagated to the current layer.

    k::Int64   = 0
    pos::Int64 = 0

    indices::Vector{Int64} = Vector{Int64}(undef,0) # ligne bidon pour faire une déclaration des indices service->client

    # Defining and initializing the heap
    Tas::Vector{zInd} = Vector{zInd}(undef,data.nJ) # Déclaration et allocation de la variable Tas
    for i in 1:data.nJ
        Tas[i] = zInd(0,0,0,0) # Allocation mémoire pour chaque élément du tas
    end
    tailleTas::Int64 = 0

    somme::Int64      = 0 # Déclaration de la somme des indices (utilisée pour la condition d'arrêt de l'énumération)
    fin::Int64        = 0
    longueur::Int64   = 0
    indc::Int64       = 0
    indloc::Int64     = 0
    tailleprec::Int64 = length(J1nd[1])

    @time for i in 2:data.nI
        #verboseDev2 ? println("LAYER $i ==================================================") : nothing

        @inbounds indices = ones(Int64,length(J1nd[i])) # All indexes are initially set to 1
        @inbounds prevlabels = labeling[layer(i-1)] # labels on the previous layer
        
        # We start with 'labeling[layer(i-1)][1]', the first label on the previous layer as it is the smallest lexicographically
        # Building the heap with the sum of elements indexed by 0 and the first label on 'prevlabels'
        tailleTas = 0
        @inbounds for j in 1:length(J1nd[i])
            tailleTas += 1
            @inbounds Tas[tailleTas].z1    = prevlabels[1].CA[1] + data.c1[i,J1nd[i][j]];
            @inbounds Tas[tailleTas].z2    = prevlabels[1].CA[2] + data.c2[i,J1nd[i][j]];
            @inbounds Tas[tailleTas].value = 10000*Tas[tailleTas].z1 + Tas[tailleTas].z2;
            @inbounds Tas[tailleTas].ind   = j;
            k = tailleTas
            while @inbounds ( (k > 1) && (Tas[div(k,2)].value > Tas[k].value) )
                TasEchange!(Tas,k,div(k,2))
                k = div(k,2)
            end
        end


        currlabels = @inbounds labeling[layer(i)] # labels on the current layer
        longueur = 0

        # Lancement de la boucle principale -----------------------------------
        @inbounds somme = length(J1nd[i])
        @inbounds fin = length(J1nd[i]) * (tailleprec + 1)
        while (somme < fin)

            # On prend toujours l'élément en tête pour instancier le nouveu label (qui sera le plus petit lexicographiquement des labels restant à générer)
            @inbounds ztete = Tas[1]

            # ---------------------------------------------------------
            # 1. test with a local bound (known on the current box)
            if ( (ztete.z1 > b.y21[1]) || (ztete.z2 > b.y12[2]) )
                # new label is out of the bounds of its reduced box => discarded
            #    verboseDev2 ? println("test 1 true: OUT OF BOX") : nothing
            else

                # Test de dominance du nouveau label avant insertion éventuelle
                @inbounds if ( (longueur == 0) || (ztete.z2 < currlabels[longueur].CA[2]) )
                    #println("Insertion de la tête du tas")
                    longueur += 1
                    
                    if longueur > maxLabels

                        # reallocation
                        resize!(labeling[1], 2*maxLabels)
                        resize!(labeling[2], 2*maxLabels) 
                        println("reallocation a priori")
                        for i in maxLabels+1:2*maxLabels
                            labeling[1][i]=Label([typemax(Int64),typemax(Int64)],[typemax(UInt16)]) 
                            labeling[2][i]=Label([typemax(Int64),typemax(Int64)],[typemax(UInt16)])      
                        end 
                        maxLabels*=2

                    end
                    
                    @inbounds currlabels[longueur].CA[1]     = ztete.z1
                    @inbounds currlabels[longueur].CA[2]     = ztete.z2
                    @inbounds currlabels[longueur].listJ1[1] = UInt16(ztete.ind)
                end

            end

            indc = ztete.ind
            @inbounds indices[indc] += 1
            @inbounds indloc = indices[indc]
            somme += 1

            if @inbounds (indloc <= tailleprec)

                # Remplacement du premier élément du tas (en place) -----------
                @inbounds label = labeling[layer(i-1)][indloc]
                @inbounds Tas[1].z1    = label.CA[1] + data.c1[i,J1nd[i][indc]]
                @inbounds Tas[1].z2    = label.CA[2] + data.c2[i,J1nd[i][indc]]
                @inbounds Tas[1].value = 10000*Tas[1].z1 + Tas[1].z2
                @inbounds Tas[1].ind   = indc
                TasDescend!(Tas, 1, tailleTas)

            else

                # Suppression de la tête du tas sans remplacement -------------
                @inbounds Tas[1].z1    = Tas[tailleTas].z1
                @inbounds Tas[1].z2    = Tas[tailleTas].z2
                @inbounds Tas[1].value = Tas[tailleTas].value
                @inbounds Tas[1].ind   = Tas[tailleTas].ind
                tailleTas -= 1
                TasDescend!(Tas, 1, tailleTas)

            end
        end

        tailleprec = longueur
    end


    # -----------------------------------------------------------------------------
    # Extract all the labels available on the last layer --------------------------

    allocation :: Vector{Int64} = zeros(Int64,nIfull)

    # set the dominant allocations into the solution
    for a in dominantallocations
        allocation[a[1]] = a[2]
    end
    @show allocation    
    for ilabel in 1:longueur
        label = @inbounds labeling[layer(data.nI)][ilabel]
        @inbounds push!(all_YN, ( label.CA[1] , label.CA[2] ) )
    end

    return nothing
end

# ==============================================================================
function labelingPaving!(data::Instance, paving::Vector{Box})

    all_YN = (Tuple{Int64, Int64})[] # temporary mechanism to collect the ND points
    all_Sol = (Solution)[] # temporary mechanism to collect the solutions    
    verboseDev2 ? println("Number of boxes to develop: ", length(paving)) : nothing
    numbox=0

    # The datastructure for representing the non-dominated labels in the graph of allocation ------------
    maxLabels = 50000
    labeling = Array{Vector{Label}}(undef,2)
    labeling[1] = Vector{Label}(undef,maxLabels)
    labeling[2] = Vector{Label}(undef,maxLabels)   
    println("allocation a priori")
    @time for i in 1:maxLabels
        labeling[1][i]=Label([typemax(Int64),typemax(Int64)],[typemax(UInt16)]) 
        labeling[2][i]=Label([typemax(Int64),typemax(Int64)],[typemax(UInt16)])      
    end  

    # Develop all the boxes in the paving -------------------------------------
    for ib in eachindex(paving)
        numbox+=1
        #verboseDev2 ?
        println("  Develop box number: $numbox", paving[ib].J1) #: nothing
        if length(paving[ib].J1) > 1
            labelingOneBox!(data, paving[ib], labeling, maxLabels, all_YN)
        else
            # box reduced to a point => added to all_YN and nothing to do
            push!(all_YN, ( paving[ib].y12[1], paving[ib].y12[2] ) )
            if withTraceAllocation
                x = fill(paving[ib].J1[1], data.nJ)
                s = paving[ib].J1
                y = [ paving[ib].y12[1], paving[ib].y12[2] ]
                #@show  x, s, y
                push!(all_Sol, Solution(x,s,y))
            end
        end
    end 

    verboseDev2 ? println("all_YN", all_YN) : nothing

    # Extract the global set of nondominated points from the local sets of nondominated points (from all boxes)
    ND_YN = getNonDominatedPoints(all_YN)
    #verboseDev2 ? println("ND_YN", ND_YN) : nothing

    return ND_YN, all_YN, all_Sol
end