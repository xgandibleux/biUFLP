
# Code didactique qui
# 1) genere tous les sous-ensembles d'indices
# 2) affiche le graphe correspondant (procedure "manuelle")

# =============================================================================

function branchAndBound(n)

    # Generation of all subsets of indexes (number of subset = sum_{p=1}^{n} binomial(n,p) )


    # Root node of the branch-and-bound (implementation of the tree by a list FIFO of nodes)
    # - each level of the tree is composed of nodes, saved consecutively in a list named 'nodes'
    # - each element of 'nodes' is composed by a list of indexes corresponding to J1 (opened services)
    # - 'head' indicates the first node for the current level of the tree
    # - 'nexthead' indicates the first node for the next level of the tree
    # - the root of the tree is initialized with J1=∅ (i.e. [])

    nodes=(Vector{Int64})[]
    J1=(Int64)[]
    push!(nodes,J1)
    head = 1


    # k is an index on the level in the tree
    for k in 1:n

        nexthead = length(nodes)+1   # index on the (future) node serving as head in the next level 

        # j is an index on the nodes to visit at the current level of the tree 
        for j in head:length(nodes)

            if length(nodes[j])==0
                # the father is the root (no index)
                lastIndex = 0
            else
                # the father is a node with a set of indexes; keep the last index
                lastIndex = nodes[j][end]
            end  

            # generate the children nodes and place them in the next level of the tree  
            for i in lastIndex+1:n
                J1=copy(nodes[j])
                push!(J1,i)
                push!(nodes,J1)
            end
        end

        # update the index for the next level
        head = nexthead
    end

    return nodes
end


# =============================================================================

using GraphRecipes, Plots

function dessineArbreIndicesUFLPinitial(n)

    # calcule toutes les combinaisons de i elements dans n
    somme=0
    for i in 1:n
        somme = somme + binomial(n,i)
    end
    g=zeros(Int8,somme+1,somme+1)


    # Generation of all subsets of indexes (number of subset = sum_{p=1}^{n} binomial(n,p) )
    nodes = branchAndBound(n)

#=    
n = 5

nodes=(Vector{Int64})[]
J1=(Int64)[]
push!(nodes,J1)
head = 1

for k in 1:n

nexthead = length(nodes)+1
for j in head:length(nodes)
  if length(nodes[j])==0
    dernier = 0
  else
    dernier = nodes[j][end]
  end  

  for i in dernier+1:n
    J1=copy(nodes[j])
    push!(J1,i)
    push!(nodes,J1)
    #g[] = 1
  end
end
head = nexthead
println(" ")
nodes
end
=#
    # construit "en dur" la matrice representant l'arbre pour un probleme a 5 services
    g[1,2] = 1
    g[1,3] = 1
    g[1,4] = 1
    g[1,5] = 1
    g[1,6] = 1

    g[2,7] = 1
    g[2,8] = 1
    g[2,9] = 1
    g[2,10] = 1

    g[3,11] = 1
    g[3,12] = 1
    g[3,13] = 1

    g[4,14] = 1
    g[4,15] = 1

    g[5,16] = 1

    g[7,17] = 1
    g[7,18] = 1
    g[7,19] = 1

    g[8,20] = 1
    g[8,21] = 1

    g[9,22] = 1

    g[11,23] = 1
    g[11,24] = 1

    g[12,25] = 1

    g[14,26] = 1

    g[17,27] = 1
    g[17,28] = 1
        
    g[18,29] = 1
    g[20,30] = 1
       
    g[23,31] = 1

    g[27,32] = 1        


    # affiche l'arbre
    txt = copy(nodes)
    stxt=[]
    for i in eachindex(txt)
        push!(stxt, replace("{"*string(txt[i])[begin+1:end-1]*"}"," " => "" ) )
    end
    stxt[1]="{}"

    default(size=(1000, 1000))
    graphplot(g,root=:left, nodeshape=:ellipse, method=:buchheim, curves=:false, nodesize=0.08, nodecolor=:yellow, fontsize=9, names=stxt)

end

# Number of services
n=5
L=branchAndBound(n)
println(L)

dessineArbreIndicesUFLPinitial(n)