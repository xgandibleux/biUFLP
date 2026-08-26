# =============================================================================
# pavingBB.jl
# August 2022
# =============================================================================


# =============================================================================
"""
    computeRunningCost(data::Instance, J1:: Vector{Int64})

    Compute the running cost vector y^R (Definition 1 of the paper) resulting
    from opening the set of facilities J1.
    - ↓ data   : the instance (provides the running costs r1, r2)
    - ↓ J1     : indices of the open facilities
    - ↑ return : a 2-element vector CR = (sum of r1[j], sum of r2[j]) for j ∈ J1
"""
function computeRunningCost(data::Instance, J1::Vector{Int64})

    CR=zeros(Int64,2)
    for j in J1
        CR[1]+=data.r1[j]
        CR[2]+=data.r2[j]
    end

    return CR
end


# =============================================================================
"""
    evaluateTest1(bU::Box, listBE::Vector{Box})

    apply the test 1 on bU knowing listBE with 
    - bU: an unexpended box
    - listBE: a list of expended boxes 
"""
function evaluateTest1(bU::Box, listBE::Vector{Box})

    verboseDev ? print(bU.CR," <-> ") : nothing
    for ib in eachindex(listBE)

        verboseDev ? println(listBE[ib].y12, " ", listBE[ib].y21) : nothing

        # Test if CR ∈ bU is in the cone of one feasible point y12,y21 of b1 ∈ listBE 
        if (isDominates(listBE[ib].y12, bU.CR)) || (isDominates(listBE[ib].y21, bU.CR))
            # found y12 or y21 of a bE who D/W/= CR of bU
            verboseDev ? println("  pruning condition by y12 or y21 found => bU (",bU.J1,") pruned !!! \n") : nothing
            return true
        end

    end
    verboseDev ? println("x") : nothing

    return false
end


# =============================================================================
"""
    expand!(data::Instance, b::Box)

    Expand a box in computing its characteristic points y12, y21, yI, yN
"""
function expand!(data::Instance, b::Box)

    # vector of costs corresponding to the lexOptimal allocation
    b.CA12=zeros(Int64,2)
    b.CA21=zeros(Int64,2)    

    for i in 1:data.nI

        # decide of the allocation for the user i -----------------------------
        vmin12 = typemax(Int64); jmin12 = -1
        vmin21 = typemax(Int64); jmin21 = -1        

        # search the service with minimal cost for user i for f1 and for f2
        for j in b.J1
            if data.c1[i,j]<vmin12
                vmin12 = data.c1[i,j];  jmin12 = j
            end
            if data.c2[i,j]<vmin21
                vmin21 = data.c2[i,j];  jmin21 = j
            end  
        end    
        b.CA12[1]+=data.c1[i,jmin12];  b.CA12[2]+=data.c2[i,jmin12] 
        b.CA21[1]+=data.c1[i,jmin21];  b.CA21[2]+=data.c2[i,jmin21]
    end
     
    b.y12 = b.CR + b.CA12
    b.y21 = b.CR + b.CA21
    b.yI  = Vector{Int64}(undef, 2);   b.yI[1] = b.y12[1];   b.yI[2] = b.y21[2]     
    b.yN  = Vector{Int64}(undef, 2);   b.yN[1] = b.y21[1];   b.yN[2] = b.y12[2]
end


# =============================================================================
"""
    evaluateTest2(bE::Box, listBE::Vector{Box})

    apply the test 2: yI of bE dominated by  y12,y21 of b ∈ listBE ?
"""
function evaluateTest2(bE::Box, listBE::Vector{Box})

    verboseDev ? print(bE.yI," <-> ") : nothing
    for ib in eachindex(listBE)

        verboseDev ? println(listBE[ib].y12, " ", listBE[ib].y21) : nothing

        # Test if yI ∈ bE is in the cone of one feasible point y12,y21 of b1 ∈ listBE 
        if (isDominates(listBE[ib].y12, bE.yI)) || (isDominates(listBE[ib].y21, bE.yI))
            # found y12 or y21 of b' ∈ listBE who D/W/= yI of bE
            verboseDev ? println("  pruning condition by y12 or y21 found => bE (",bE.J1,") pruned !!! \n") : nothing
            return true
        end

    end
    verboseDev ? println("x") : nothing

    return false    
end


# =============================================================================
"""
    evaluateTest3!(bE::Box, listBE::Vector{Box}, nbPruned::Dict{Symbol, Int64})

    apply the test 3: yI of b ∈ listBE dominated by  y12,y21 of bE ?
"""
function evaluateTest3!(bE::Box, listBE::Vector{Box}, nbPruned::Dict{Symbol, Int64})

    verboseDev ? print(bE.y12, " ", bE.y21," <-> ") : nothing


    for ib in reverse(2:length(listBE))
        verboseDev ? println(listBE[ib].yI) : nothing

        # Test if yI of b' ∈ listBE is in the cone of one feasible point y12,y21 of bE 
        if (isDominates(bE.y12, listBE[ib].yI)) || (isDominates(bE.y21, listBE[ib].yI))
            # found y12 or y21 of a bE who D/W/= yI of b' ∈ listBE
            verboseDev ? println("  pruning condition by y12 or y21 found => bE (",bE.J1," ",ib,") pruned !!! \n") : nothing
            nbPruned[:TEST3]+=1
            deleteat!(listBE,ib)
        end

    end

    return nothing
end


# =============================================================================
"""
    computePavingBranchAndBound(data::Instance)

    compute a paving with the branch-and-bound algorithm following this principle:
    - it is based on a breadth-first search
    - the tree is implemented by two lists of nodes, 'listL' of unexpended nodes and 'listB' of expended nodes
    - generation of (at most) all subsets of indexes (number of subset = sum_{p=1}^{n} binomial(n,p) )
    - application of test1, test2, test3 to prune useless (unexpended or expended) nodes
"""
function computePavingBranchAndBound(data::Instance)

    # Preliminary definitions and datastructures:
    # - unexpended node: a box defined only by J1,CR
    # - expended node: a box defined by J1,CR,CA12,CA21,yI,yN,y12,y21

    listL = (Box)[]   # list of unexpended nodes
    listB = (Box)[]   # list of expended nodes    

    # vector for counting the number of times that a test has been triggered
    nbPruned = Dict(:TEST1 => 0 , :TEST2 => 0, :TEST3 => 0)


    # Root node of the branch-and-bound ---------------------------------------
    # - each level of the tree is composed of nodes, saved consecutively in a list named 'listL'
    # - each element of 'listL' is composed by a list of indexes corresponding to J1 (opened services)
    # - 'head' indicates the first node for the current level of the tree
    # - 'nexthead' indicates the first node for the next level of the tree
    # - the root of the tree is initialized with J1=∅ (i.e. [])

    # INIT: create the root node which is an unexpended node
    J1    = (Int64)[]
    CR    = [0,0]

    # ENQUEUE: add in the tree the root node 
    push!(listL, Box(J1, CR)) 
    
    # setup the index on the first node for the current level of the tree
    head = 1

    # All others nodes of the branch-and-bound --------------------------------

    # k is an index on the level in the tree
    for k in 1:data.nJ+1 # off-by-one bug fixed (+1)

        nexthead = length(listL)+1   # index on the (future) node serving as head in the next level 

        # DEQUEUE: ------------------------------------------------------------        
        # j is an index on the node from listL to be processed at the current level of the tree 
        for j in head:length(listL)
           # for lll in 1:length(listL)
           #     println(" >>> ", listL[lll].J1)#$head  $j")
           # end
           # println("")

            # TEST 1 ----------------------------------------------------------
            # test if CR is weakly dominated

            test1 = false

            verboseDev ? print(j," test 1: ") : nothing
            if head == 1
                # test1 is not triggered (root node has no predecessor)
                verboseDev ? println("-") : nothing
            else
                test1 = evaluateTest1(listL[j],listB)
            end

            if test1

                # the child nodes are NOT generated                
                nbPruned[:TEST1]+=1

            else
                # the child nodes are generated

                # get the last index in J1 of the father node
                if length(listL[j].J1)==0
                    # the father is the root (no index)
                    lastIndex = 0
                else
                    # the father is a node with a set of indexes; keep the last index
                    lastIndex = listL[j].J1[end]
                end  

                # generate the children nodes and place them in the next level of the tree  
                for i in lastIndex+1:data.nJ

                    # create an unexpended node
                    # `copy` (shallow) is enough and faster than `deepcopy` here:
                    # J1 is a flat Vector{Int64} with no nested mutable structure
                    # to protect against aliasing.
                    J1 = copy(listL[j].J1)
                    push!(J1,i)
                    CR = computeRunningCost(data,J1) # TODO: possible improvement here, discuss with the Fernandez-Puerto bound from the literature

                    # ENQUEUE: add in the tree the current node 
                    push!(listL, Box(J1, CR))

                end

                # -------------------------------------------------------------
                # EXPAND: 
                if head == 1
                    # expand is not triggered (root node cannot be expanded)
                    verboseDev ? println("-") : nothing
                else
                    # bU becomes bE
                    # Every element of `listL` is created with the 2-argument
                    # Box(J1,CR) constructor (see above and the root init),
                    # which always leaves CA12, CA21, yI, yN, y12, y21, listyNS
                    # as empty vectors. Reconstructing via the same constructor
                    # (with copies of J1 and CR) reproduces exactly what
                    # `deepcopy(listL[j])` would produce, without the overhead
                    # of deep-copying 7 vectors that are always empty here.
                    b = Box(copy(listL[j].J1), copy(listL[j].CR))
                    expand!(data::Instance, b::Box)

                    verboseDev ? print(j," test 2: ") : nothing
                    # TEST 2 --------------------------------------------------
                    # test if b is weakly dominated
                    test2 = evaluateTest2(b,listB)
                    
                    if test2

                        # the expanded nodes is killed                
                        nbPruned[:TEST2]+=1
                    
                    else
                        # ENQUEUE: add in the tree ('listB') the current expended node 'b'
                        pushfirst!(listB, b)

                        verboseDev ? print(j," test 3: ") : nothing
                        # TEST 3 ----------------------------------------------
                        # test if b' ∈ listB is weakly dominated by b
                        evaluateTest3!(b,listB,nbPruned) 
        
                    end

                end
                
            end


        end

        # update the index for the next level
        head = nexthead
    end

    # summary of the activity -------------------------------------------------  
    verboseProd ? println("        #test 1 activated : ", nbPruned[:TEST1]) : nothing 
    verboseProd ? println("        #test 2 activated : ", nbPruned[:TEST2]) : nothing  
    verboseProd ? println("        #test 3 activated : ", nbPruned[:TEST3]) : nothing     
   
    verboseProd ? println("\n        #max of boxes           : ", 2^data.nJ-1) : nothing 
    verboseProd ? println("        #active boxes           : ", length(listB)) : nothing 
    verboseProd ? println("        #unexpended boxes       : ", length(listL)) : nothing        
    verboseProd ? println("        ratio of boxes pruned   : ", round((1.0-length(listB)/length(listL))*100; digits=2)  , "%") : nothing  
 
    return listB
end