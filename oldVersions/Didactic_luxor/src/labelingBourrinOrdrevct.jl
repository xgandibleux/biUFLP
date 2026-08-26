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

# Lexicographic order on (z1,z2): z1 dominates the comparison, z2 breaks ties.
# Compared directly as a tuple rather than packed into a single scalar, so the
# order is correct regardless of the magnitude of z1/z2 (no hidden assumption
# on the cost range of the instance).
isless(a::zInd,b::zInd) = (a.z1, a.z2) < (b.z1, b.z2)

# Maximum number of nondominated labels that can be stored per layer during the
# labeling algorithm (phase 3, generation). This is a heuristic capacity, not a
# formally derived bound — see labelingOneBox! for the explicit check raised
# when this capacity is exceeded, and labelingPaving! where the buffers are
# preallocated to this size.
const maxLabelsPerLayer = 250_000


# Version with only 2 layers and one list of ND labels per layer


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
    J1nd = (Int64)[vctCosts[1,3]] # after the sort, the first couple is not dominated => save the value of service j

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
    ensureLabelInitialized!(labeling::Array{Vector{Label}}, maxInitializedLabel::Vector{Int64}, layerIdx::Int64, idx::Int64)

    Lazily construct the Label object at `labeling[layerIdx][idx]` on its first
    use. `labeling[layerIdx]` is created with `Vector{Label}(undef, ...)`
    (uninitialized slots), so this must be called before mutating the fields of
    `labeling[layerIdx][idx]` for the first time. Once constructed, an object
    is reused (its fields mutated in place) across all subsequent boxes that
    reach that index — the placeholder values written here are never read
    before being overwritten, so any valid Int64/UInt16 works.

    `maxInitializedLabel[layerIdx]` tracks the highest index ever initialized
    in that layer, across all boxes processed so far in the current
    labelingPaving! call. Within one box's processing of one layer, indices are
    always written consecutively from 1 upward (iVctLabels only ever
    increases by 1 before each write), so a plain integer comparison against
    this high-water mark is equivalent to (and cheaper than) an `isassigned`
    check on every single label write — the difference matters once most of
    the buffer ends up used, i.e. for instances with a very large number of
    nondominated points (e.g. dataset H).
"""
function ensureLabelInitialized!(labeling::Array{Vector{Label}}, maxInitializedLabel::Vector{Int64}, layerIdx::Int64, idx::Int64)
    if idx > maxInitializedLabel[layerIdx]
        # placeholder matches the original eager-init sentinel: never read as a
        # meaningful value (always overwritten before use), except in the
        # defensive read below, where a huge value is naturally filtered out
        # by the box-bound test — exactly as the original code tolerated it.
        labeling[layerIdx][idx] = Label([typemax(Int64),typemax(Int64)],[typemax(UInt16)])
        maxInitializedLabel[layerIdx] = idx
    end
    return nothing
end


# ==============================================================================
"""
    labelingOneBox!(dataFull::Instance, b::Box, labeling::Array{Vector{Label}}, iVctLabels::Vector{Int64}, maxInitializedLabel::Vector{Int64}, all_YN)

    Solve the bi-objective allocation subproblem associated to box `b` with a
    two-layer labeling algorithm (dynamic programming over customers), and push
    the resulting locally nondominated points into `all_YN`.
    - ↓ dataFull            : the full instance (a working copy is reduced in place)
    - ↓ b                   : the box to develop (its J1 defines the open facilities)
    - ↕ labeling            : preallocated buffer of Label vectors, reused across boxes (2 layers)
    - ↕ iVctLabels          : number of labels currently stored in each of the 2 layers
    - ↕ maxInitializedLabel : highest index ever constructed in each of the 2 layers (see ensureLabelInitialized!)
    - ↕ all_YN              : accumulator of nondominated points (Tuple{Int64,Int64}) found so far
"""
function labelingOneBox!(dataFull::Instance, b::Box, labeling::Array{Vector{Label}}, iVctLabels::Vector{Int64}, maxInitializedLabel::Vector{Int64}, all_YN)

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
        iVctLabels[layer(i)] > maxLabelsPerLayer && error(
            "labelingOneBox!: number of nondominated labels for box J1=$(b.J1) " *
            "exceeds the preallocated capacity (maxLabelsPerLayer=$maxLabelsPerLayer) " *
            "on user $i. Increase maxLabelsPerLayer in labelingBourrinOrdrevct.jl and retry."
        )
        ensureLabelInitialized!(labeling, maxInitializedLabel, layer(i), iVctLabels[layer(i)])
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

    indices::Vector{Int64} = Vector{Int64}(undef,data.nJ) # per-source pointer into J1nd[i]; sized on nJ as a safe upper bound
    heap::MutableBinaryMinHeap{zInd} = MutableBinaryMinHeap{zInd}() # binary min-heap merging the |J1nd[i]| candidate label streams
    generated::Int64 = 0 # number of candidate labels generated so far for user i (loop stop condition)
    zTop::zInd    = zInd(0,0,0)
    handle::Int64 = 0

     for i in 2:data.nI
        #verboseDev2 ? println("LAYER $i ==================================================") : nothing
        indices = ones(Int64,length(J1nd[i])) # all pointers start at position 1

        # defensive: guarantee labeling[layer(i-1)][1] exists before reading it
        # below, even in the (believed unreachable, but not proven) edge case
        # where the previous layer ended up with zero labels.
        ensureLabelInitialized!(labeling, maxInitializedLabel, layer(i-1), 1)

        heap = MutableBinaryMinHeap{zInd}([ zInd( labeling[layer(i-1)][1].CA[1] + data.c1[i,j] , 
                                                 labeling[layer(i-1)][1].CA[2] + data.c2[i,j] , 
                                                 pos
                                               ) for (pos,j) in enumerate(J1nd[i])
                                         ]
                                        ) # heap initialized with the first candidate of each source

        # initialize the list of nondominated labels for the user i
        iVctLabels[layer(i)] = 0

        # main k-way merge loop
        generated = length(J1nd[i])
        totalToGenerate = length(J1nd[i]) * (iVctLabels[layer(i-1)] + 1)

        while (generated < totalToGenerate)
            # always take the head of the heap to instantiate the new label (the lexicographically smallest among the remaining candidates)
            zTop, handle = top_with_handle(heap)

            # ---------------------------------------------------------
            # 1. test with a local bound (known on the current box)
            if ( (zTop.z1 > b.y21[1]) || (zTop.z2 > b.y12[2]) )
                # new label is out of the bounds of its reduced box => discarded
                verboseDev2 ? println("test 1 true: OUT OF BOX") : nothing
            else
                # dominance test of the new label before its possible insertion
                if ( (iVctLabels[layer(i)] == 0) || (zTop.z2 < labeling[layer(i)][iVctLabels[layer(i)]].CA[2]) )
 
                    iVctLabels[layer(i)]+=1
                    iVctLabels[layer(i)] > maxLabelsPerLayer && error(
                        "labelingOneBox!: number of nondominated labels for box J1=$(b.J1) " *
                        "exceeds the preallocated capacity (maxLabelsPerLayer=$maxLabelsPerLayer) " *
                        "on user $i. Increase maxLabelsPerLayer in labelingBourrinOrdrevct.jl and retry."
                    )
                    ensureLabelInitialized!(labeling, maxInitializedLabel, layer(i), iVctLabels[layer(i)])
                    labeling[layer(i)][iVctLabels[layer(i)]].CA[1] = zTop.z1
                    labeling[layer(i)][iVctLabels[layer(i)]].CA[2] = zTop.z2
                    labeling[layer(i)][iVctLabels[layer(i)]].listJ1[1] = UInt16(0)         

                end
            end

            indices[zTop.ind] += 1
            generated += 1
            if (indices[zTop.ind] <= iVctLabels[layer(i-1)])

                update!(heap, handle, zInd( labeling[layer(i-1)][indices[zTop.ind]].CA[1] + data.c1[i,J1nd[i][zTop.ind]], 
                                           labeling[layer(i-1)][indices[zTop.ind]].CA[2] + data.c2[i,J1nd[i][zTop.ind]] , 
                                           zTop.ind
                                         )
                       )
            else
                pop!(heap)
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
"""
    labelingPaving!(data::Instance, paving::Vector{Box})

    Apply the labeling algorithm (phase 3, generation) to every box of the
    (refined) paving, and merge the locally nondominated points obtained in
    each box into the global nondominated set YN of the 2-UFLP instance.
    - ↓ data   : the instance to solve
    - ↓ paving : the paving obtained after phases 1 (paving) and 2 (refinement)
    - ↑ return : (ND_YN, all_YN) — ND_YN is the filtered set of nondominated
                 points, all_YN is the raw (unfiltered, possibly dominated)
                 collection of points gathered across all boxes
"""
function labelingPaving!(data::Instance, paving::Vector{Box})

    all_YN = (Tuple{Int64, Int64})[] # temporary solution to collect the ND points
    verboseDev2 ? println("Number of boxes to develop: ", length(paving)) : nothing
    numbox=0

    # The datastructure for representing the non-dominated labels in the graph of allocation ------------
    # maxLabelsPerLayer (defined above) is a heuristic upper bound on the number
    # of nondominated labels expected per layer for the instances handled by this
    # solver (see dataset F and H in the paper); it is not derived from a formal
    # bound on the instance size. labelingOneBox! raises a clear error if an
    # instance ever produces more nondominated labels in one layer than this
    # capacity, instead of silently indexing out of bounds.
    #
    # PERFORMANCE NOTE: the two buffers are only *sized* to maxLabelsPerLayer
    # here; individual Label objects are constructed lazily on first use (see
    # ensureLabelInitialized!) instead of eagerly filling all 2*maxLabelsPerLayer
    # slots upfront. Eager initialization used to dominate the total runtime for
    # instances with few nondominated points (most of dataset F), since it pays
    # for ~500,000 small allocations regardless of how many labels the instance
    # actually needs. maxInitializedLabel tracks, per layer, the highest index
    # already constructed (a plain integer high-water mark, cheaper than an
    # `isassigned` check on every single label write once most of the buffer
    # ends up used — relevant for the very large instances of dataset H).
    labeling = Array{Vector{Label}}(undef,2)
    labeling[1] = Vector{Label}(undef,maxLabelsPerLayer)
    labeling[2] = Vector{Label}(undef,maxLabelsPerLayer)
    iVctLabels = [0,0]
    maxInitializedLabel = [0,0]

    # Develop all the boxes in the paving -------------------------------------
    for ib in eachindex(paving)
        numbox+=1
        verboseDev2 ? println("  Develop box number: $numbox", paving[ib].J1) : nothing
        if length(paving[ib].J1) > 1
            labelingOneBox!(data, paving[ib], labeling, iVctLabels, maxInitializedLabel, all_YN)
        else
            # box reduced to a point => added to all_YN and nothing to do
            push!(all_YN, (paving[ib].y12[1], paving[ib].y12[2]) )
        end
    end 

    verboseDev2 ? println("all_YN", all_YN) : nothing

    # Extract the global set of nondominated points from the local sets of nondominated points (from all boxes)
    ND_YN = getNonDominatedPoints(all_YN)

    return ND_YN, all_YN
end
