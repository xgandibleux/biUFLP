# ==============================================================================
# display graphically the results

print("Loading and compiling the graphic package (Luxor, MathTeXEngine)...")
using Luxor               # dessin vectoriel
using MathTeXEngine       # affichage avec commandes latex 
println(" done!")

# ==============================================================================
# parameters to control the drawings

const draw_paving             = true

const draw_allRunningCosts    = true
const draw_CRplusCA12andCA21  = true
const draw_yNSbox             = true
const draw_yN                 = false


# ==============================================================================
"""
    coordNO(p::Point)

    transform the coordinates of a point to correspond to the reference (0,0) at top-left of luxor
"""
function coordNO(p::Point)

    x = p.x
    y = (CanvasSize) - p.y

    return Point(x,y)

end


# ==============================================================================
"""
    displayBoxJ1(b, name, color, legend)

    draw a box
    - legend gets :legend or :nolegend
"""
function displayBoxJ1(b::Box, name::String, color, legend)

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


# ==============================================================================
function drawGraphics(paving,YN)

    # prologue ----------------------------------------------------------------

    # erase the file .png if it already exists 
    if isfile("figUFLPbox.png")
	    rm("figUFLPbox.png")
    end

    # parameters of the drawing area
    Drawing(CanvasSize, CanvasSize, "figUFLPbox.png")
    background("white")

    # draw the axes x et y 
    sethue("grey")
    arrow( coordNO(Point(50, 30)), coordNO(Point(50, CanvasSize)) )
    arrow( coordNO(Point(30, 50)), coordNO(Point(CanvasSize, 50)) )

    # add the names of axes  
    fontsize(12)
    text(L"(50;50)", coordNO(Point(25, 30)), halign=:center, valign=:baseline, angle=0)
    fontsize(16)
    text(L"f^1(x)", coordNO(Point(CanvasSize-30, 30)), halign=:center, valign=:baseline, angle=0)
    text(L"$f^2(x)$", coordNO(Point(30, CanvasSize-30)), halign=:center, valign=:baseline, angle=0)

    # draw a paving -----------------------------------------------------------
    if draw_paving

        for ib in eachindex(paving)
            displayBoxJ1(paving[ib], "B"*string(ib),"black", :nolegend)
        end

    end


    # draw all CR -------------------------------------------------------------

    if draw_allRunningCosts

        setopacity(1.0)
        sethue("purple")
        for b in eachindex(paving)
            circle(coordNO(Point(paving[b].CR[1], paving[b].CR[2])), 3, action = :fill)
        end

    end


    # draw CR + CA12/CA21 for a box -------------------------------------------

    # index of a box
    ib1 = 3

    if draw_CRplusCA12andCA21

        sethue("red")
        arrow(coordNO(Point(0,0)),coordNO(Point(paving[ib1].CR[1], paving[ib1].CR[2])))#, action = :stroke) 
        arrow(coordNO(Point(paving[ib1].CR[1], paving[ib1].CR[2])),coordNO(Point(paving[ib1].CR[1]+paving[ib1].CA12[1], paving[ib1].CR[2]+paving[ib1].CA12[2])))#, action = :stroke) 
        arrow(coordNO(Point(paving[ib1].CR[1], paving[ib1].CR[2])),coordNO(Point(paving[ib1].CR[1]+paving[ib1].CA21[1], paving[ib1].CR[2]+paving[ib1].CA21[2])))#, action = :stroke) 
        displayBoxJ1(paving[ib1], "B"*string(ib1),"red", :nolegend)

    end


    # draw yNS available for a box --------------------------------------------

    # index of a box
    ib = 3

    if draw_yNSbox

        #for ib in eachindex(paving)
            setopacity(1.0)
            sethue("red")
            for ip in eachindex(paving[ib].listyNS)
                circle(coordNO(Point(paving[ib].listyNS[ip][1], paving[ib].listyNS[ip][2])), 3, action = :fill) 
            end
        #end

    end


    # draw Y_N ----------------------------------------------------------------

    if draw_yN

        setopacity(1.0)
        sethue("green2")
        for i in eachindex(YN)
            circle(coordNO(Point(YN[i][1], YN[i][2])), 4, action = :fill)
        end

    end

end