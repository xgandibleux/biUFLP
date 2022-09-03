# =============================================================================
# datastru.jl
# August 2022
# =============================================================================


# =============================================================================
# Global constants

# -----------------------------------------------------------------------------
# size of the drawing area (square)
const CanvasSize = 700

# -----------------------------------------------------------------------------
# depth limitation on the recursivity when yNS for a box are computed
const depthlMax = 8


# =============================================================================
# structure of an instance
mutable struct Instance
    fname :: String        # name of the file
    nI    :: Int64         # number of users
    nJ    :: Int64         # number of services
    c1    :: Array{Int,2}  # assignment costs users-services for objective 1
    c2    :: Array{Int,2}  # assignment costs users-services for objective 2
    r1    :: Array{Int,1}  # running costs for objective 1
    r2    :: Array{Int,1}  # running costs for objective 2
end


# =============================================================================
# structure of a box
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
# structure of a label
mutable struct Label
    CA     :: Vector{Int64}   # (partiel) allocation cost 
    #listJ1 :: Vector{Int64}   # (partial) allocation of user-service
    #
    function Label() # inner constructor
        label        = new()
        label.CA     = [0,0]
        #label.listJ1 = []
        return label
    end
    #
    function Label(CA)#,listJ1) # inner constructor
        label        = new()
        label.CA     = CA
        #label.listJ1 = listJ1
        return label
    end
end