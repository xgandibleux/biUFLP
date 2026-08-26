# ==============================================================================
# fname : vopt.jl
# August 2022 - Revision 2026
# ==============================================================================

print("Loading and compiling JuMP, Gurobi, MultiObjectiveAlgorithms...")
using JuMP, Gurobi
import MultiObjectiveAlgorithms as MOA
println(" done!")


# ==============================================================================
"""
    createProblem2UFLP_JuMP(optimizer::DataType, data::Instance)

    Create the JuMP model of 2UFLP 
"""
function createProblem2UFLP_JuMP(optimizer::DataType, data::Instance)

    model = Model(() -> MOA.Optimizer(optimizer))    
    
    @variable(model, x[1:data.nI,1:data.nJ], Bin)
    @variable(model, s[1:data.nJ], Bin)

    @expression( model, f1, sum(data.c1[i,j]*x[i,j] for i in 1:data.nI, j in 1:data.nJ) + sum(data.r1[j]*s[j] for j in 1:data.nJ) )
    @expression( model, f2, sum(data.c2[i,j]*x[i,j] for i in 1:data.nI, j in 1:data.nJ) + sum(data.r2[j]*s[j] for j in 1:data.nJ) )
    @objective( model, Min, [f1,f2])

    @constraint( model, [i=1:data.nI], sum(x[i,j] for j in 1:data.nJ) == 1 )
    @constraint( model, [i=1:data.nI,j=1:data.nJ], x[i,j] <= s[j] )    

    return model

end


# ==============================================================================
"""
    vOpt_MOA(data::Instance, mipOptimizerToUse::String)

    Compute Y_N with MOA
"""
function vOpt_MOA(data::Instance, mipoptimizerToUse::String)

    # -------------------------------------------------------------------------
    # Information about the instance to solve

    verboseProd ? println("\n    [A] INSTANCE \n") : nothing

    verboseProd ? println("        filename      : $(data.fname)") : nothing
    verboseProd ? println("        nI (users)    : $(data.nI)") : nothing
    verboseProd ? println("        nJ (services) : $(data.nJ)") : nothing


    # -------------------------------------------------------------------------
    # compute with vOptGeneric a paving of the instance

    verboseProd ? println("\n    [B] PAVING STAGE WITH JuMP+MOA (MIP=$mipOptimizerToUse) \n") : nothing

    optimizer = Gurobi.Optimizer
    time_limit::Float64 = 600.0

    # Define the 2UFLP model
    mod2UFLP = createProblem2UFLP_JuMP(optimizer, data)

    set_time_limit_sec(mod2UFLP, time_limit)
    set_silent(mod2UFLP)

    # Optimize the model with MOA using an ϵ-constraint method
    set_attribute(mod2UFLP, MOA.Algorithm(), MOA.EpsilonConstraint())
    #set_attribute(mod2UFLP, MOA.Algorithm(), MOA.TambyVanderpooten())

    getTime = time()
    optimize!(mod2UFLP)
    elapsedTime = round(time()- getTime, digits=4)

    # Get the results
    Y_N = (Vector{Int64})[]
    print("        filename      : $(data.fname) : ")
    if elapsedTime < time_limit  #is_solved_and_feasible(mod2UFLP)
        for i in 1:result_count(mod2UFLP)
            #print(i, ": x = ", round.(Int, value.(x; result = i)), " | ")
            #println("$i : z = ", round.(Int, objective_value(mod2UFLP; result = i)))
            push!(Y_N, round.(Int, objective_value(mod2UFLP; result = i)) )
        end
        println(result_count(mod2UFLP), " | ", elapsedTime)
    else
        Y_N = []
        println(" -1", " | ", elapsedTime)
    end
    

    # summary of the activity -------------------------------------------------  

    #@show Y_N
    verboseProd ? println("\n        #YN       : ", length(Y_N)) : nothing       
    verboseProd ? println("\n        Time(MOA) : ", elapsedTime, " sec \n") : nothing

    return Y_N, elapsedTime

end


# ========== dirty EpsilonConstraint

# ==============================================================================
"""
    build2UFLP(optimizer, data::Instance)
    Build a single-objective JuMP model (no MOA wrapper).
"""
function build2UFLP(optimizer, data::Instance)

    model = Model(optimizer)                       
    @variable(model, x[1:data.nI, 1:data.nJ], Bin)
    @variable(model, s[1:data.nJ], Bin)

    @expression(model, f1,
        sum(data.c1[i,j]*x[i,j] for i in 1:data.nI, j in 1:data.nJ) +
        sum(data.r1[j]*s[j]     for j in 1:data.nJ))
    @expression(model, f2,
        sum(data.c2[i,j]*x[i,j] for i in 1:data.nI, j in 1:data.nJ) +
        sum(data.r2[j]*s[j]     for j in 1:data.nJ))

    @constraint(model, [i=1:data.nI], sum(x[i,j] for j in 1:data.nJ) == 1)
    @constraint(model, [i=1:data.nI, j=1:data.nJ], x[i,j] <= s[j])

    return model
end


# ==============================================================================
"""
    epsilonConstraint(data::Instance, mipOptimizerToUse::String)

    Compute Y_N with a direct ε-constraint method (no MOA wrapper), filtering
    dominated points inline as the RHS ε decreases.
"""
function epsilonConstraint(data::Instance, mipOptimizerToUse::String)

    verboseProd ? println("\n    [A] INSTANCE \n") : nothing
    verboseProd ? println("        filename      : $(data.fname)") : nothing
    verboseProd ? println("        nI (users)    : $(data.nI)") : nothing
    verboseProd ? println("        nJ (services) : $(data.nJ)") : nothing

    optimizer   = Gurobi.Optimizer
    time_limit  = 600.0

    mod = build2UFLP(optimizer, data)
    set_time_limit_sec(mod, time_limit)
    MOI.set(mod, MOI.Silent(), true)

    Y_N         = (Vector{Int64})[]    
    getTime     = time()

# ------------------------------------------------------------------
# STEP 1 — min f2  →  gives z2_min  
# ------------------------------------------------------------------
    @objective(mod, Min, mod[:f2])
    optimize!(mod)
    if !is_solved_and_feasible(mod) || (time() - getTime >= time_limit)
        return [], -1.0
    end
    z2_min       = round(Int, value(mod[:f2]))
    extreme_right = [round(Int, value(mod[:f1])), z2_min]   # to add at the end

# ------------------------------------------------------------------
# STEP 2 — min f1  →  gives z2_max 
# ------------------------------------------------------------------
    @objective(mod, Min, mod[:f1])
    optimize!(mod)
    if !is_solved_and_feasible(mod) || (time() - getTime >= time_limit)
        return [], -1.0
    end
    z2_max = round(Int, value(mod[:f2]))
    push!(Y_N, [round(Int, value(mod[:f1])), z2_max])

# ------------------------------------------------------------------
# STEP 3 — main loop: min f1  s.t.  f2 <= ε, ε decreasing, filtering inline
# ------------------------------------------------------------------
    epsilon = z2_max - 1
    eps_con = @constraint(mod, mod[:f2] <= epsilon)

    while epsilon >= z2_min && (time() - getTime) < time_limit
        set_normalized_rhs(eps_con, epsilon)
        set_time_limit_sec(mod, time_limit - (time() - getTime))
        optimize!(mod)
        if !is_solved_and_feasible(mod) || (time() - getTime >= time_limit)
            return [], -1.0
        end
        z1_new = round(Int, value(mod[:f1]))
        z2_new = round(Int, value(mod[:f2]))
        if z1_new == Y_N[end][1] && z2_new < Y_N[end][2]
            Y_N[end] = [z1_new, z2_new]   # same value on f1, better value on f2 → replace
        else
            push!(Y_N, [z1_new, z2_new])  # different value on f1 → add
        end
        epsilon = z2_new - 1
    end

# ------------------------------------------------------------------    
# Add the last extreme point (min f2) if required
# ------------------------------------------------------------------
    if Y_N[end] != extreme_right
        if extreme_right[1] == Y_N[end][1]
            Y_N[end] = extreme_right      # same on f1 → replace
        else
            push!(Y_N, extreme_right)     # different value on f1 → add
        end
    end

    elapsedTime = round(time() - getTime, digits=4)

    print("        filename      : $(data.fname) : ")
    println(length(Y_N), " | ", elapsedTime)
    verboseProd ? println("\n        #YN       : ", length(Y_N)) : nothing
    verboseProd ? println("\n        Time(ε-c) : ", elapsedTime, " sec \n") : nothing

    if elapsedTime >= time_limit
        return [], -1.0
    else
        return Y_N, elapsedTime
    end
end