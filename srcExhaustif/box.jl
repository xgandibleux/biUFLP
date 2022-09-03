# =============================================================================
# box.jl
# August 2022
# =============================================================================


# ==============================================================================
"""
    computeWeight(b::Box)

    For a given box, compute the weight in direction of the ideal point 
"""
function computeWeight(b::Box)

    λ1 = b.y12[2]-b.y21[2]
    λ2 = b.y21[1]-b.y12[1]

    return λ1/(λ1+λ2)

end


# ==============================================================================
"""
    splitBox(b::Box, yλ::Vector{Int64})

    Given a point yλ in the box b, split the box in two subboxes b_u and b_l
"""
function splitBox(b::Box, yλ::Vector{Int64}, CAλ::Vector{Int64})

    # b_u, upper box created
    yI    = Vector{Int64}(undef, 2);  
    yI[1] = b.y12[1];  
    yI[2] = yλ[2]     
    
    yN    = Vector{Int64}(undef, 2);  
    yN[1] = yλ[1];  
    yN[2] = b.y12[2]
    b_u   = Box(b.J1, b.CR, b.CA12, CAλ, yI, yN, b.y12, yλ)

    # b_l, lower box created
    yI    = Vector{Int64}(undef, 2);  
    yI[1] = yλ[1];  
    yI[2] = b.y21[2]     
    
    yN    = Vector{Int64}(undef, 2);  
    yN[1] = b.y21[1];  
    yN[2] = yλ[2]
    b_l   = Box(b.J1, b.CR, CAλ, b.CA21, yI, yN, yλ, b.y21)

    return b_u, b_l
end