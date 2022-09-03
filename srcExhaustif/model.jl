print("Loading and compiling the bi-objective optimization model (vOptGeneric, JuMP, GLPK)...")
using vOptGeneric, JuMP, GLPK
println(" done!")


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

