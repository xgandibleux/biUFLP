# ==============================================================================

# variable globale !! pour besoins d'affichage graphique
decomposition_yNS=[]


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

    println(depth, espace, "yI=", b_cur.yI, "  yN=", b_cur.yN)

    if depth > depthlMax
        println(depth, espace, "Recursivity limit reached")
        # recursive depth limit reached => not decomposable; nothing to do with b_cur 
        return nothing
    end
    # the limited recursive depth is not reached

    if  b_cur.y12 ==  b_cur.y21
        println(depth, espace, "Box reduced to a unique point")
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
    println(depth, espace, "yNS=", yNS)

    # check if the obtained point is not already known ------------------------
    if  (b_cur.y12 == yNS) || (b_cur.y21 == yNS) 
        println(depth, espace, "no new point found")
        # nothing to do with b_cur => not decomposable; returned as is
        return nothing
    end
    # the point obtained is a new one => the box is decomposable

    # store the nondominated supported point
    push!(b_ref.listyNS, yNS)

    # split the box based on yNS ----------------------------------------------
    b_u,b_l = splitBox(b_cur, yNS, CAλ)

    #push!(decomposition_yNS, b_u)   # global variable!! 
    #push!(decomposition_yNS, b_l)   # global variable!! 

    println(depth, espace, "upper")
    computeyNS!(data, b_ref, b_u, depth+1)
    println(depth, espace, "lower")    
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