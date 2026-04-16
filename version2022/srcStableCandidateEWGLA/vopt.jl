# ==============================================================================
# fname : vopt.jl
# August 2022
# ==============================================================================

print("Loading and compiling vOptGeneric, JuMP, GLPK, Gurobi...")
using vOptGeneric, JuMP, GLPK, Gurobi
println(" done!")


# ==============================================================================
"""
    createProblem2UFLP(solver::DataType, data::Instance)

    Create the vOptGeneric model of 2UFLP 
"""
function createProblem2UFLP(solver::DataType, data::Instance)

    model = vModel( solver )
    
    @variable(model, x[1:data.nI,1:data.nJ], Bin)
    @variable(model, s[1:data.nJ], Bin)

    @addobjective( model, Min, sum(data.c1[i,j]*x[i,j] for i in 1:data.nI, j in 1:data.nJ) + sum(data.r1[j]*s[j] for j in 1:data.nJ) )
    @addobjective( model, Min, sum(data.c2[i,j]*x[i,j] for i in 1:data.nI, j in 1:data.nJ) + sum(data.r2[j]*s[j] for j in 1:data.nJ) )

    @constraint( model, [i=1:data.nI], sum(x[i,j] for j in 1:data.nJ) == 1 )
    @constraint( model, [i=1:data.nI,j=1:data.nJ], x[i,j] <= s[j] )    

    return model

end


# ==============================================================================
"""
    vOpt(data::Instance, mipSolverToUse::String)

    Compute Y_N with vOptSolver
"""
function vOpt(data::Instance, mipSolverToUse::String)

    # -------------------------------------------------------------------------
    # Information about the instance to solve

    verboseProd ? println("\n    [A] INSTANCE \n") : nothing

    verboseProd ? println("        filename      : $(data.fname)") : nothing
    verboseProd ? println("        nI (users)    : $(data.nI)") : nothing
    verboseProd ? println("        nJ (services) : $(data.nJ)") : nothing


    # -------------------------------------------------------------------------
    # compute with vOptGeneric a paving of the instance

    verboseProd ? println("\n    [B] PAVING STAGE WITH VOPTSOLVER (MIP=$mipSolverToUse) \n") : nothing

    solver = GLPK.Optimizer
    if mipSolverToUse == "GUROBI"
        solver = Gurobi.Optimizer
    end

    getTime = time()

    # Define the 2UFLP model
    mod2UFLP = createProblem2UFLP(solver, data)
    set_silent(mod2UFLP)

    # Optimize the model with vOptSolver using an ϵ-constraint method
    vSolve(mod2UFLP, method=:epsilon, step = 1.0)

    # Get the results
    Y_N = getY_N( mod2UFLP )

    # summary of the activity -------------------------------------------------  

    verboseProd ? println("\n        #YN               : ", length(Y_N)) : nothing       
    verboseProd ? println("\n        Time(vOptGeneric) : ", round(time()- getTime, digits=4), " sec \n") : nothing

    #printX_E( mod2UFLP )

    return Y_N

end