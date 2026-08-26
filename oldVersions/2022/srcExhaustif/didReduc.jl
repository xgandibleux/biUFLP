# =============================================================================
# illustration didactique de situation
# - de dominance (D) ou
# - domainance faible (W) ou
# - aucune dominance immediate entre 2 boites 


# packages usites
print("Loading and compiling the graphic package...")
using Luxor               # dessin vectoriel
using MathTeXEngine       # affichage avec commandes latex 
println(" done! \nRunning... \n")


# -----------------------------------------------------------------------------
# size of the drawing area (square)
const CanvasSize = 500


# -----------------------------------------------------------------------------
# structure of a box
mutable struct Box
    yI  :: Tuple{Int64, Int64}       # ideal point
    yN  :: Tuple{Int64, Int64}       # nadir point
    y12 :: Tuple{Int64, Int64}       # optLex(1;2) point
    y21 :: Tuple{Int64, Int64}       # optLex(2;1) point
    #
     Box(yI   :: Tuple{Int64,Int64},  
         yN   :: Tuple{Int64,Int64}, 
         y12  :: Tuple{Int64,Int64}, 
         y21  :: Tuple{Int64,Int64}) = new(yI,yN,y12,y21) 
end 

# -----------------------------------------------------------------------------
# transforme la coordonnee du point pour correspondre a la reference (0,0) haut-gauche de luxor

function coordNO(p::Point)

    x = p.x
    y = (CanvasSize - 1) - p.y

    return Point(x,y)

end


# -----------------------------------------------------------------------------
# cree une boite

function createBox(x,y,h,w)

    yI = (x, y)     # ideal
    yN = (x+w, y+h) # nadir
    y12 = (x,y+h)   # lex 1|2
    y21 = (x+w,y)   # lex 2|1

    return Box(yI, yN, y12, y21)

end


# -----------------------------------------------------------------------------
# trace d'une boite 

function displayBox(x,y, h,w, name, color)

    # contour et nom de la boite
	sethue(color)
	setdash("dash")
	box( coordNO(Point(x, y)), coordNO(Point(x+w, y+h)), action=:stroke)
	text(name, coordNO(Point(x+w/2, y+h/2)), halign=:center, valign=:baseline, angle=0)
		
    # points caracteristiques (ideal nadir) de la boite
	setdash("solid")
	circle(coordNO(Point(x, y)), 3, action = :fill)
	text(L"y^I", coordNO(Point(x-10, y-20)), halign=:left, valign=:baseline, angle=0)
	circle(coordNO(Point(x+w, y+h)), 3, action = :fill)
	text(L"y^N", coordNO(Point(x+w+10, y+h+10)), halign=:center, valign=:baseline, angle=0)
	
    # points correspondants aux solutions lexicographiquement optimales
	box(coordNO(Point(x,y+h)), 4, 4, action=:stroke)
	text(L"y^{12}", coordNO(Point(x-10, y+h+10)), halign=:left, valign=:baseline, angle=0)
	box(coordNO(Point(x+w,y)), 4, 4, action=:stroke)
	text(L"y^{21}", coordNO(Point(x+w+10, y-20)), halign=:center, valign=:baseline, angle=0)

    return nothing

end


# =============================================================================
# Main

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
arrow( coordNO(Point(50, 30)), coordNO(Point(50, CanvasSize-30)) )
arrow( coordNO(Point(30, 50)), coordNO(Point(CanvasSize-30, 50)) )

# place les reperes sur les axes ----------------------------------------------
fontsize(16)
text(L"(0;0)", coordNO(Point(30, 30)), halign=:center, valign=:baseline, angle=0)
text(L"f^1(x)", coordNO(Point(CanvasSize-50, 30)), halign=:center, valign=:baseline, angle=0)
text(L"$f^2(x)$", coordNO(Point(30, CanvasSize-50)), halign=:center, valign=:baseline, angle=0)

# boite de reference ----------------------------------------------------------
h0=150; w0=150      # hauteur et largeur
x0=150; y0=150      # point inferieur gauche (ideal) de la boite 

b0 = createBox(x0,y0,h0,w0)
displayBox(x0,y0,h0,w0,"B0","black")
box(coordNO(Point(x0+25,y0+100)), 4, 4, action=:stroke)
box(coordNO(Point(x0+50,y0+50)), 4, 4, action=:stroke)
box(coordNO(Point(x0+100,y0+25)), 4, 4, action=:stroke)

# boite candidate -------------------------------------------------------------

# taille et positionnement aleatoire
h1=115#rand(50:150)              # hauteur et largeur
w1=98#rand(50:150)              
x1=191#rand(50:CanvasSize-w1-20)  # point inferieur gauche (ideal) de la boite 
y1=195#rand(50:CanvasSize-h1-20)

b1 = createBox(x1,y1,h1,w1)
displayBox(x1,y1,h1,w1,"B1","blue")
box(coordNO(Point(x1+h1/3,y1+w1/3)), 4, 4, action=:stroke)
#box(coordNO(Point(x1+h1/6,y1+w1/0.85)), 4, 4, action=:stroke)

sethue("red")
setopacity(0.15)


#If one rectangle is on left side of other
if (b0.y12[1] >= b1.y21[1] || b1.y12[1] >= b0.y21[1])
    println("Pas de chevauchement")
else
    # If one rectangle is above other
    if (b0.y21[2] >= b1.y12[2] || b1.y21[2] >= b0.y12[2])
        println("Pas de chevauchement")
    else
        println("Chevauchement")
    end
end


# analyse ---------------------------------------------------------------------

sethue("red")
setopacity(0.15)


if     (Point(x0+25,y0+100).x < b1.yI[1]) &&  (Point(x0+25,y0+100).y < b1.yI[2])
    box( coordNO(Point(x0+25,y0+100)), coordNO(Point(CanvasSize, CanvasSize)), action=:fill)

elseif (Point(x0+50,y0+50).x < b1.yI[1]) &&  (Point(x0+50,y0+50).y < b1.yI[2])
    box( coordNO(Point(x0+50,y0+50)), coordNO(Point(CanvasSize, CanvasSize)), action=:fill)

elseif (Point(x0+100,y0+25).x < b1.yI[1]) &&  (Point(x0+100,y0+25).y < b1.yI[2])
    box( coordNO(Point(x0+100,y0+25)), coordNO(Point(CanvasSize, CanvasSize)), action=:fill)
end

#= exemple redA redB
h1=151#rand(50:150)              # hauteur et largeur
w1=99#rand(50:150)              
x1=317#rand(50:CanvasSize-w1-20)  # point inferieur gauche (ideal) de la boite 
y1=83#rand(50:CanvasSize-h1-20)
box( coordNO(Point(x0+100,y0+25)), coordNO(Point(CanvasSize, CanvasSize)), action=:fill)

sethue("blue")
setopacity(0.5)
box(coordNO(Point(x1+h1/6,y1+w1/0.85)), 4, 4, action=:stroke)
box( coordNO(Point(x1+h1/6,y1+w1/0.85)), coordNO(Point(b1.y21[1],b1.y21[2])), action=:fill)
=#

#= exemple red0 red1 red2
h1=148#rand(50:150)              # hauteur et largeur
w1=118#rand(50:150)              
x1=190#rand(50:CanvasSize-w1-20)  # point inferieur gauche (ideal) de la boite 
y1=237#rand(50:CanvasSize-h1-20)
box(coordNO(Point(x1+h1/3,y1+w1/3)), 4, 4, action=:stroke)

box( coordNO(Point(x0+25,y0+100)), coordNO(Point(CanvasSize, CanvasSize)), action=:fill)
sethue("blue")
setopacity(0.5)
setdash("dash")
box( coordNO(Point(x1+h1/3,y1+w1/3)), coordNO(Point(b1.y21[1],b1.y21[2])), action=:fill)

sethue("red")
setopacity(0.15)
box( coordNO(Point(x0+50,y0+50)), coordNO(Point(CanvasSize, CanvasSize)), action=:fill)
=#



#=
if (b0.y12[1] < b1.y12[1]) &&  (b0.y21[1] > b1.y12[1])
    sethue("green")
    box( coordNO(Point(b1.y12[1],b0.y21[2])), coordNO(Point(b1.y21[1], b1.y21[2])), action=:fill)
end
=#

#=
if (b0.y21[2] < b1.y12[2]) &&  (b0.yN[2] < b1.yI[2])
    box( coordNO(Point(b1.y12[1],b0.y21[2])), coordNO(Point(b1.y21[1], b1.y21[2])), action=:fill)
elseif (b0.y12[1] < b1.y21[1]) &&  (b0.y12[2] > b1.y12[2])
    sethue("green")
    box( coordNO(Point(b1.y12[1],b1.y12[2])), coordNO(Point(b0.y12[1], b0.y12[2])), action=:fill)
elseif (b1.y12[1] < b0.y21[1]) &&  (b0.y12[2] > b1.y21[2])
    sethue("green")
    box( coordNO(Point(b1.y12[1],b1.y12[2])), coordNO(Point(b0.y12[1], b0.y12[2])), action=:fill)
end
=#

#=
# conclusions pouvant etre immediatement deduites 
if (b0.yN[1] < b1.yI[1]) &&  (b0.yN[2] < b1.yI[2])
    box( coordNO(Point(b0.yN[1],b0.yN[2])), coordNO(Point(CanvasSize, CanvasSize)), action=:fill)
    println("B0 D B1  car  b0.yN=$(b0.yN) D b1.yI=$(b1.yI)")
    
elseif (b0.yN[1] == b1.yI[1]) &&  (b0.yN[2] == b1.yI[2])
    box( coordNO(Point(b0.yN[1],b0.yN[2])), coordNO(Point(CanvasSize, CanvasSize)), action=:fill)
    println("B0 W B1  car  b0.yN=$(b0.yN) = b1.yI=$(b1.yI)")  
    
elseif (b0.y12[1] < b1.yI[1]) &&  (b0.y12[2] < b1.yI[2])
    box( coordNO(Point(b0.y12[1],b0.y12[2])), coordNO(Point(CanvasSize, CanvasSize)), action=:fill)
    println("B0 D B1  car  b0.y12=$(b0.y12) D b1.yI=$(b1.yI)")
    
elseif (b0.y12[1] == b1.yI[1]) &&  (b0.y12[2] < b1.yI[2])
    box( coordNO(Point(b0.y12[1],b0.y12[2])), coordNO(Point(CanvasSize, CanvasSize)), action=:fill)
    println("B0 W B1  car  b0.y12=$(b0.y12) W b1.yI=$(b1.yI)") 
    
elseif (b0.y21[1] < b1.yI[1]) &&  (b0.y21[2] < b1.yI[2])
    box( coordNO(Point(b0.y21[1],b0.y21[2])), coordNO(Point(CanvasSize, CanvasSize)), action=:fill)    
    println("B0 D B1  car  b0.y21=$(b0.y21) D b1.yI=$(b1.yI)")

elseif (b0.y21[1] < b1.yI[1]) &&  (b0.y21[2] == b1.yI[2])
    box( coordNO(Point(b0.y21[1],b0.y21[2])), coordNO(Point(CanvasSize, CanvasSize)), action=:fill)    
    println("B0 W B1  car  b0.y21=$(b0.y21) W b1.yI=$(b1.yI)")
        
else
    println("pas de D|W immediate : b0=($(b0.yI);$(b0.yN))  b1=($(b1.yI);$(b1.yN))")
end
=#

# epilogue --------------------------------------------------------------------
finish()
preview()
