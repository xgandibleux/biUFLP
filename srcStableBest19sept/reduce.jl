# =============================================================================
# fname : reduce.jl
# September 2022
# =============================================================================


# -----------------------------------------------------------------------------      
"""
    computeLocalDominatedAllocation!(data::Instance, b::Box)

    1) Compute the nondominated allocations for all users
    2) Allocate the dominant allocation and compute (k1,k2) the corresponding constant
    3) Remove into the allocation costs the lines corresponding to the dominant allocation
"""
function computeLocalDominatedAllocation!(data::Instance, b::Box)

    # Compute the nondominated allocations for all users ----------------------
    J1nd = computeListNondominatedAllocationJ1allUsers(data, b.J1)

    iAllocated = findall(x->length(x)==1,J1nd)
    jOpened = Set{Int64}([])
    for i in eachindex(iAllocated)
        push!(jOpened,J1nd[i][1])
    end
    #println("nb local single allocations user-service: ", length(iAllocated))


    # Compute the allocation cost for all dominant allocations ----------------
    k1::Int64=0
    k2::Int64=0

    for (pos,i) in enumerate(iAllocated)
        j = J1nd[i][1]
        #println("Set allocation i=$i j=$j ",data.c1[i,j]," ", data.c2[i,j])
        k1 += data.c1[i,j] 
        k2 += data.c2[i,j] #+ data.r2[j]
    end


    # Reduce the cost matrices consequently -----------------------------------

    data.c1 = data.c1[setdiff(1:end, iAllocated),:]
    data.c2 = data.c2[setdiff(1:end, iAllocated),:]
    data.nI -= length(iAllocated)
    J1nd = J1nd[setdiff(1:end, iAllocated)]

    return  J1nd, k1, k2
end


# ----------------------------------------------------------------------------- 
function computeGlobalDominatedAllocation(data::Instance)

    # Compute the nondominated allocations for all users ----------------------
    J1nd = computeListNondominatedAllocationJ1allUsers(data, collect(1:data.nJ))

    
    iAllocated = findall(x->length(x)==1,J1nd)
    jOpened = Set{Int64}([])
    for i in eachindex(iAllocated)
        push!(jOpened,J1nd[i][1])
    end
    #@show iAllocated
    #@show jOpened

    #k1::Int64=0
    #k2::Int64=0
    #all_YN=(Tuple{Int64, Int64})[]

    #for (pos,i) in enumerate(iAllocated)
    #    j=J1nd[i][1]
        #println("Set allocation i=$i j=$j ",data.c1[i,j]," ", data.c2[i,j])
        #k1 = data.c1[i,j] 
        #k2 = data.c2[i,j] #+ data.r2[j]
    #end

    println("nb global single allocations user-service: ", length(iAllocated))

    #println( data.c1)
    #println( data.c2)

    #data.c1 = data.c1[setdiff(1:end, iAllocated),:]
    #data.c2 = data.c2[setdiff(1:end, iAllocated),:]
    #println( data.c1)
    #println( data.c2)


    return  nothing
end