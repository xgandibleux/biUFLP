# ==============================================================================

print("Loading and compiling the combinatoric package (Combinatorics)...")
using Combinatorics
println(" done!")


"""
    generateExhaustivePaving(data::Instance)

    Exhaustive generation of the 2^n-1 boxes giving a paving
"""
function generateExhaustivePaving(data::Instance)

    # the paving 
    paving = (Box)[]

    # genere les couts d'ouverture pour toutes les combinaisons de services -------

    ib ::Int64 = 0
    for J1 in combinations(collect(1:data.nJ))

        # vecteur des couts d'ouverture des services J1
        CR=zeros(Int,2)
        for j in J1
            CR[1]+=data.r1[j]
            CR[2]+=data.r2[j]
        end

        # vecteur des couts d'affectation lexOptimale
        CA12=zeros(Int,2)
        CA21=zeros(Int,2)    
        for i in 1:data.nI
            # traite l'usager i
            vmin12 = typemax(Int64); jmin12 = -1
            vmin21 = typemax(Int64); jmin21 = -1        
            # recherche le service de cout min pour l'usager i pour f1 et pour f2 
            for j in J1
                if data.c1[i,j]<vmin12
                    vmin12 = data.c1[i,j];  jmin12 = j
                end
                if data.c2[i,j]<vmin21
                    vmin21 = data.c2[i,j];  jmin21 = j
                end  
            end    
            CA12[1]+=data.c1[i,jmin12]; CA12[2]+=data.c2[i,jmin12] 
            CA21[1]+=data.c1[i,jmin21]; CA21[2]+=data.c2[i,jmin21]
        end
 
        y12 = CR+CA12
        y21 = CR+CA21
        yI  = Vector{Int64}(undef, 2); yI[1]=y12[1]; yI[2]=y21[2]     
        yN  = Vector{Int64}(undef, 2); yN[1]=y21[1]; yN[2]=y12[2]

        push!(paving, Box(J1, CR, CA12, CA21, yI, yN, y12, y21))   
        ib+=1

        print("$ib  J1=$J1  CR=($(CR[1]);$(CR[2]))  CA12=($(CA12[1]);$(CA12[2]))  CA21=($(CA21[1]);$(CA21[2]))")  
        println("  ||  y12=($(y12[1]);$(y12[2]))  y21=($(y21[1]);$(y21[2]))  yI=($(yI[1]);$(yI[2]))  yN=($(yN[1]);$(yN[2]))")
  
    end

    return paving

end