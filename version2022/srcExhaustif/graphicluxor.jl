print("Loading and compiling the graphic package (Luxor, MathTeXEngine)...")
using Luxor               # dessin vectoriel
using MathTeXEngine       # affichage avec commandes latex 
println(" done!")


# -----------------------------------------------------------------------------
# transforme la coordonnee du point pour correspondre a la reference (0,0) haut-gauche de luxor

function coordNO(p::Point)

    x = p.x
    y = (CanvasSize) - p.y

    return Point(x,y)

end


# -----------------------------------------------------------------------------
# draw a box 

function displayBoxJ1(b, name, color, legend)

    # contour et nom de la boite
	sethue(color)
    setline(1)
	box( coordNO(Point(b.yI[1], b.yI[2])), coordNO(Point(b.yN[1], b.yN[2])), action=:stroke)
    if legend == :legend
    	text(name, coordNO(Point(b.yI[1]+(b.yN[1]-b.yI[1])/2, b.yI[2]+(b.yN[2]-b.yI[2])/2)), halign=:center, valign=:baseline, angle=0)
    end
		
    # points caracteristiques (ideal nadir) de la boite
	setdash("solid")
	circle(coordNO(Point(b.yI[1], b.yI[2])), 3, action = :fill)
    if legend == :legend
    	text(L"y^I", coordNO(Point(b.yI[1]-10, b.yI[2]-20)), halign=:left, valign=:baseline, angle=0)
    end
	circle(coordNO(Point(b.yN[1], b.yN[2])), 3, action = :fill)
    if legend == :legend
        text(L"y^N", coordNO(Point(b.yN[1]+10, b.yN[2]+10)), halign=:center, valign=:baseline, angle=0)
    end
	
    # points correspondants aux solutions lexicographiquement optimales
	box(coordNO(Point(b.y12[1],b.y12[2])), 4, 4, action=:stroke)
    if legend == :legend
        text(L"y^{12}", coordNO(Point(b.y12[1]-10, b.y12[2]+10)), halign=:left, valign=:baseline, angle=0)
    end
	box(coordNO(Point(b.y21[1],b.y21[2])), 4, 4, action=:stroke)
    if legend == :legend
        text(L"y^{21}", coordNO(Point(b.y21[1]+10, b.y21[2]-20)), halign=:center, valign=:baseline, angle=0)
    end

    return nothing

end

function displayBoxJ1unique(b, name, color, legend)

    # contour et nom de la boite
	sethue(color)
    setline(1)
	box( coordNO(Point(b.yI[1], b.yI[2])), coordNO(Point(b.yN[1], b.yN[2])), action=:stroke)
    if legend == :legend
    	text(name, coordNO(Point(b.yI[1]+(b.yN[1]-b.yI[1])/2, 5+b.yI[2]+(b.yN[2]-b.yI[2])/2)), halign=:center, valign=:baseline, angle=0)
    end
		#=
    # points caracteristiques (ideal nadir) de la boite
	setdash("solid")
	circle(coordNO(Point(b.yI[1], b.yI[2])), 3, action = :fill)
    if legend == :legend
    	text(L"y^I", coordNO(Point(b.yI[1]-20, b.yI[2]-30)), halign=:left, valign=:baseline, angle=0)
    end
	circle(coordNO(Point(b.yN[1], b.yN[2])), 3, action = :fill)
    if legend == :legend
        text(L"y^N", coordNO(Point(b.yN[1]+20, b.yN[2]+20)), halign=:center, valign=:baseline, angle=0)
    end
	=#
    # points correspondants aux solutions lexicographiquement optimales
	box(coordNO(Point(b.y12[1],b.y12[2])), 4, 4, action=:stroke)
    #=
    if legend == :legend
        text(L"y^{12}", coordNO(Point(b.y12[1]-20, b.y12[2]+20)), halign=:left, valign=:baseline, angle=0)
    end
	box(coordNO(Point(b.y21[1],b.y21[2])), 4, 4, action=:stroke)
    if legend == :legend
        text(L"y^{21}", coordNO(Point(b.y21[1]+20, b.y21[2]-30)), halign=:center, valign=:baseline, angle=0)
    end
=#
    return nothing

end


# prologue --------------------------------------------------------------------

# supprime le fichier image si celui-ci existe
if isfile("figUFLPbox.png")
	rm("figUFLPbox.png")
end

# parametres de la zone de dessin
Drawing(CanvasSize, CanvasSize, "figUFLPbox.png")
background("white")


# trace les axes x et y -------------------------------------------------------
sethue("grey")
arrow( coordNO(Point(50, 30)), coordNO(Point(50, CanvasSize)) )
arrow( coordNO(Point(30, 50)), coordNO(Point(CanvasSize, 50)) )

# place les reperes sur les axes ----------------------------------------------
fontsize(12)
text(L"(50;50)", coordNO(Point(25, 30)), halign=:center, valign=:baseline, angle=0)
fontsize(16)
text(L"f^1(x)", coordNO(Point(CanvasSize-30, 30)), halign=:center, valign=:baseline, angle=0)
text(L"$f^2(x)$", coordNO(Point(30, CanvasSize-30)), halign=:center, valign=:baseline, angle=0)

if draw_paving
    for ib in eachindex(paving)
        #if !pruned[ib]
        displayBoxJ1(paving[ib], "B"*string(ib),"black", :nolegend)
        #end
    end
end


# visualise CR ----------------------------------------------------------------

if draw_allRunningCosts

    setopacity(1.0)
    sethue("purple")
    for b in eachindex(paving)
        circle(coordNO(Point(paving[b].CR[1], paving[b].CR[2])), 3, action = :fill)
    end

end


# visualise CR + CA12/CA21 pour une boite -------------------------------------

# indices de 2 boites particulieres
ib1 = 2# 25
ib2 = 1 

#for ib1=6:6
if draw_CRplusCA12andCA21

    sethue("red")
    arrow(coordNO(Point(0,0)),coordNO(Point(paving[ib1].CR[1], paving[ib1].CR[2])))#, action = :stroke) 
    arrow(coordNO(Point(paving[ib1].CR[1], paving[ib1].CR[2])),coordNO(Point(paving[ib1].CR[1]+paving[ib1].CA12[1], paving[ib1].CR[2]+paving[ib1].CA12[2])))#, action = :stroke) 
    arrow(coordNO(Point(paving[ib1].CR[1], paving[ib1].CR[2])),coordNO(Point(paving[ib1].CR[1]+paving[ib1].CA21[1], paving[ib1].CR[2]+paving[ib1].CA21[2])))#, action = :stroke) 
    displayBoxJ1(paving[ib1], "B"*string(ib1),"red", :nolegend)
    #displayBoxJ1unique(paving[ib1], "B"*string(ib1),"red", :legend)
    #displayBoxJ1(paving[ib2], "B"*string(ib2),"blue", :nolegend)

end
#end
#@assert false "stop"

# visualise les yNS calcules pour une boite -----------------------------------

if draw_yNS

    setopacity(1.0)
    sethue("red")
    for ip in eachindex(paving[ib1].listyNS)
        circle(coordNO(Point(paving[ib1].listyNS[ip][1], paving[ib1].listyNS[ip][2])), 5, action = :fill) 
    end

    sethue("blue")
    for ip in eachindex(paving[ib2].listyNS)
        circle(coordNO(Point(paving[ib2].listyNS[ip][1], paving[ib2].listyNS[ip][2])), 5, action = :fill) 
    end
end

# visualise pour une boite la decomposition en sous-boites realisee  ----------

if draw_decompositionBox

    # accede a une variable globale !!
    for ib in eachindex(decomposition_yNS)
        displayBoxJ1(decomposition_yNS[ib], "B"*string(ib),"red", :nolegend)
    end

end


# visualise Y_N ---------------------------------------------------------------

if draw_nondominatedPointsALGO

    setopacity(1.0)
    sethue("green4")
    for i in eachindex(all_YN)
        circle(coordNO(Point(all_YN[i][1], all_YN[i][2])), 2, action = :fill)
    end
    for i in eachindex(ND_YN)
        circle(coordNO(Point(ND_YN[i][1], ND_YN[i][2])), 8, action = :fill)
    end

end

if draw_nondominatedPoints

    setopacity(1.0)
    sethue("yellow")
    for i in eachindex(Y_N)
        circle(coordNO(Point(Y_N[i][1], Y_N[i][2])), 4, action = :fill)
    end

end




# visualise BoundSet_yI -------------------------------------------------------

if draw_EBpointsIdeaux

    setopacity(1.0)
    sethue("blue")
    #setdash("dash")
    for i in eachindex(ND_yI)
        circle(coordNO(Point(ND_yI[i][1], ND_yI[i][2])), 3, action = :fill)
    end
    poly([coordNO(Point(BoundSet_yI[i][1],BoundSet_yI[i][2])) for i in eachindex(BoundSet_yI)], :stroke)

end


# epilogue --------------------------------------------------------------------
finish()
preview()