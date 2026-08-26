# ==============================================================================
# fname : mopRoutines.jl
# August 2022
# ==============================================================================


# ==============================================================================
"""
    isDominates(yA::Vector{Int64}, yB::Vector{Int64})

    Test if yA stronglyDominates/weaklyDominates/isEqualTo yB
"""
function isDominates(yA::Vector{Int64}, yB::Vector{Int64}) :: Bool

    ifelse((@inbounds (yA[1] - yB[1]) ≤ 0) && (@inbounds (yA[2] - yB[2]) ≤ 0), true, false)

end



# -----------------------------------------------------------------------------
"""
    getNonDominatedPoints(lstPoints::Vector{Tuple{Int64, Int64}})

    Algorithm of Kung: it extracts S_N from a static set of points S of IR^2
    Remark: the weakly nondominated points are filtered in this implementation
"""
function getNonDominatedPoints(lstPoints::Vector{Tuple{Int64, Int64}})

    SN = (Tuple{Int64, Int64})[] 

    if length(lstPoints) > 1

        sort!(lstPoints, by = x -> x[1])
        #@show lstPoints
        push!(SN, lstPoints[1]) ; minpoints_y = lstPoints[1][2]

        for i = 2:length(lstPoints)
            if  (lstPoints[i][2] < minpoints_y)
                if (SN[end][1] == lstPoints[i][1]) && (SN[end][2] > lstPoints[i][2])
                    # weakly nondominated point identified => removed from SN
                    #println("\nFND identifie : ", SN[end], " ", lstPoints[i])
                    #@show SN
                    pop!(SN)
                    #@show SN
                    #println(" ")
                end 
                push!(SN, lstPoints[i]) ; minpoints_y = lstPoints[i][2]
                #@show SN
            end
        end
        #@show SN
    end

    return SN

end


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


# ==============================================================================
"""
    computeyNS!(data::Instance, b_ref::Box, b_cur::Box, depth::Int64)

    Compute recursively some nonsupported points for a given box.
    The limitation is a parameter on the depth of the recursivity. 
    The root box given in paramater is modified.

    The (private) recursive function
    - data   : the data of the instance
    - b_ref  : the parent box for which the computation of yNS is applied
    - b_cur  : the current box candidate to be split
    - depth  : current depth of recursivity
"""
function computeyNS!(data::Instance, b_ref::Box, b_cur::Box, depth::Int64)

    espace :: String = ""
    for i=1:depth
        espace=espace*"  "
    end

    verboseDev ? println(depth, espace, "yI=", b_cur.yI, "  yN=", b_cur.yN) : nothing

    if depth > depthlMax
        verboseDev ? println(depth, espace, "Recursivity limit reached") : nothing
        # recursive depth limit reached => not decomposable; nothing to do with b_cur 
        return nothing
    end
    # the limited recursive depth is not reached

    if  b_cur.y12 ==  b_cur.y21
        verboseDev ? println(depth, espace, "Box reduced to a unique point") : nothing
        # box reduced to a unique point => not decomposable; nothing to do with b_cur 
        return nothing
    end
    # box not reduced to a unique point 

    # compute a point with λ (weighted sum of the two objectives) -----------------
    λ   = computeWeight(b_cur)
    cλ  = λ * data.c1 + (1.0-λ) * data.c2
    CAλ = zeros(Int,2)
    yNS = zeros(Int64,2)

    # determine an assignment for all users
    for i=1:data.nI
        vmin = typemax(Int64); jmin = -1
        for j in b_cur.J1
            # for a given user i, identify the smaller assigment cost(i,j) forall j∈J1
            if cλ[i,j] < vmin
                vmin = cλ[i,j]; jmin = j
            end
        end
        CAλ[1]+=data.c1[i,jmin]
        CAλ[2]+=data.c2[i,jmin]  
    end
    yNS[1] = b_cur.CR[1]+CAλ[1]
    yNS[2] = b_cur.CR[2]+CAλ[2]
    verboseDev ? println(depth, espace, "yNS=", yNS) : nothing

    # check if the obtained point is not already known ------------------------
    if  (b_cur.y12 == yNS) || (b_cur.y21 == yNS) 
        verboseDev ? println(depth, espace, "no new point found") : nothing
        # nothing to do with b_cur => not decomposable; returned as is
        return nothing
    end
    # the point obtained is a new one => the box is decomposable

    # store the nondominated supported point
    push!(b_ref.listyNS, yNS)

    # split the box based on yNS ----------------------------------------------
    b_u,b_l = splitBox(b_cur, yNS, CAλ)

    verboseDev ? println(depth, espace, "upper") : nothing
    computeyNS!(data, b_ref, b_u, depth+1)
    verboseDev ? println(depth, espace, "lower") : nothing
    computeyNS!(data, b_ref, b_l, depth+1)
    return nothing

end


# ==============================================================================
"""
    computeyNS!(data::Instance, b::Box)

    Compute recursively some nonsupported points for a given box.
        The limitation is a parameter on the depth of the recursivity. 
        The root box given in paramater is modified.
    
        The (public) function
        - data   : the data of the instance
        - b      : the box for which the computation of yNS is applied
"""
function computeyNS!(data::Instance, b::Box)
    depth :: Int64 = 1
    computeyNS!(data, b, b, depth)
    return nothing
end
