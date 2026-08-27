# =============================================================================
# datastru.jl
# August 2022
# =============================================================================


# =============================================================================
# Global constants

# -----------------------------------------------------------------------------
# Maximum recursion depth allowed when computing nonsupported points (yNS) for
# a box (see computeyNS! in mopRoutines.jl). Each level of recursion splits a
# box in two, so depth 8 caps the number of yNS computed per box at 2^8-1=255.
# This is a heuristic trade-off between the quality of the refinement (phase 2)
# and its computational cost; it is not derived from a formal bound.
const depthlMax = 8


# =============================================================================
"""
    Instance

    Data of a bi-objective 0/1 uncapacitated facility location problem (2-UFLP):
    `nI` users (customers) to assign to a subset of `nJ` candidate services
    (facilities), with two independent objectives (1 and 2) each combining an
    assignment cost matrix and a facility running (opening) cost vector.
"""
mutable struct Instance
    fname :: String          # name of the file
    nI    :: Int64           # number of users
    nJ    :: Int64           # number of services
    c1    :: Array{Int64,2}  # assignment costs users-services for objective 1
    c2    :: Array{Int64,2}  # assignment costs users-services for objective 2
    r1    :: Array{Int64,1}  # running costs for objective 1
    r2    :: Array{Int64,1}  # running costs for objective 2
end


# =============================================================================
"""
    Box

    A box `B_J1` of the paving, associated to a set `J1` of open facilities
    (see Definition 2 of the paper). It is a 2-dimensional subset of the
    objective space bounded by the lexicographic optimal points `y12` and
    `y21` of the allocation subproblem defined by `J1`, together with the
    corresponding ideal point `yI`, nadir point `yN`, running cost `CR`, and
    the assignment costs `CA12`/`CA21` of the lexicographic solutions. A box
    may also carry `listyNS`, a list of nonsupported nondominated points
    computed for it during the refinement phase (phase 2).
"""
mutable struct Box
    J1   :: Vector{Int64}        # services opened
    #
    CR   :: Vector{Int64}        # running cost for J1
    CA12 :: Vector{Int64}        # assigment cost corresponding to optLex(1;2)
    CA21 :: Vector{Int64}        # assigment cost corresponding to optLex(2;1)
    #
    yI   :: Vector{Int64}        # ideal point
    yN   :: Vector{Int64}        # nadir point
    y12  :: Vector{Int64}        # optLex(1;2) point
    y21  :: Vector{Int64}        # optLex(2;1) point
    #
    listyNS :: Vector{Vector{Int64}}  # some nodominated supported points
    #
    
    function Box(J1,CR,CA12,CA21,yI,yN,y12,y21) # inner constructor
        box         = new()
        box.J1      = J1
        box.CR      = CR 
        box.CA12    = CA12
        box.CA21    = CA21     
        box.yI      = yI 
        box.yN      = yN
        box.y12     = y12
        box.y21     = y21
        box.listyNS = []
        return box
    end
    
    function Box(J1,CR) # inner constructor
        box         = new()
        box.J1      = J1
        box.CR      = CR 
        box.CA12    = []
        box.CA21    = []     
        box.yI      = [] 
        box.yN      = []
        box.y12     = []
        box.y21     = []
        box.listyNS = []
        return box
    end
end   


# =============================================================================
"""
    Label

    A (partial) label of the labeling algorithm (phase 3, generation): the
    bi-objective allocation cost `CA` accumulated so far, together with the
    (partial) user-to-facility allocation `listJ1` that produced it.
"""
 mutable struct Label
    CA     :: Vector{Int64}   # (partial) allocation cost 
    listJ1 :: Vector{UInt16}   # (partial) allocation of user-service !! Max 65535 users
    #
    function Label(CA,listJ1) # inner constructor
        label        = new()
        label.CA     = CA
        label.listJ1 = listJ1
        return label
    end
end
