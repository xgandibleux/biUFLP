# -----------------------------------------------------------------------------
# Algorithme de Kung (extrait S_N d'un ensemble statique de points S de IR^2)

function getNonDominatedPoints(lstPoints::Vector{Tuple{Int64, Int64}})
    
    sort!(lstPoints, by = x -> x[1])
    SN=(Tuple{Int64, Int64})[] ; push!(SN, lstPoints[1]) ; minpoints_y = lstPoints[1][2]
    for i=2:length(lstPoints)
        if  (lstPoints[i][2] < minpoints_y)
            if (lstPoints[i][1] == lstPoints[i-1][1])
                pop!(SN)
            end 
            push!(SN, lstPoints[i]) ; minpoints_y = lstPoints[i][2]
        end
    end

    return SN

end


# -----------------------------------------------------------------------------
# Extraction des points non-domines de l'ensemble des points yI

function extractND_yI(paving::Vector{Box})

    all_yI=(Tuple{Int64, Int64})[]
    for b in eachindex(paving)
        push!(all_yI, (paving[b].yI[1],paving[b].yI[2]) )
    end
    ND_yI = getNonDominatedPoints(all_yI)

    return ND_yI

end


# -----------------------------------------------------------------------------
# Elaboration de l'ensemble bornant constitue des points non-domines de yI

function elaborateBoundSet_yI(paving::Vector{Box}) 

    ND_yI = extractND_yI(paving)
    BoundSet_yI = [(ND_yI[1][1], CanvasSize)]
    for i in eachindex(ND_yI)
        push!(BoundSet_yI, (ND_yI[i][1], ND_yI[i][2]))
    end
    push!(BoundSet_yI, (CanvasSize, ND_yI[end][2]))

    return ND_yI, BoundSet_yI

end