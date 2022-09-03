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
	text(name, coordNO(Point(x+w/2, 5+y+h/2)), halign=:center, valign=:baseline, angle=0)
		
    # points caracteristiques (ideal nadir) de la boite
	setdash("solid")
	circle(coordNO(Point(x, y)), 3, action = :fill)
	text(L"y^I", coordNO(Point(x-20, y-30)), halign=:left, valign=:baseline, angle=0)
	circle(coordNO(Point(x+w, y+h)), 3, action = :fill)
	text(L"y^N", coordNO(Point(x+w+20, y+h+20)), halign=:center, valign=:baseline, angle=0)
	
    # points correspondants aux solutions lexicographiquement optimales
	box(coordNO(Point(x,y+h)), 4, 4, action=:stroke)
	text(L"y^{12}", coordNO(Point(x-20, y+h+20)), halign=:left, valign=:baseline, angle=0)
	box(coordNO(Point(x+w,y)), 4, 4, action=:stroke)
	text(L"y^{21}", coordNO(Point(x+w+20, y-30)), halign=:center, valign=:baseline, angle=0)

    return nothing

end

# -----------------------------------------------------------------------------
# trace d'une boite 

function displayCR(x,y, h,w, name, color)

    # contour et nom de la boite
	sethue(color)
	setdash("dash")
	box( coordNO(Point(x, y)), coordNO(Point(x+w, y+h)), action=:stroke)
    text(name, coordNO(Point(x+w/2, 5+y+h/2)), halign=:center, valign=:baseline, angle=0)
		
    # points caracteristiques (ideal nadir) de la boite
	setdash("solid")
	circle(coordNO(Point(x+w,y)), 3, action = :fill)
	text(L"y^{R}", coordNO(Point(x+w+20, y-30)), halign=:center, valign=:baseline, angle=0)

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
h=0; w=0      # hauteur et largeur
x=200; y=200      # point inferieur gauche (ideal) de la boite 

b0 = createBox(x,y,h,w)
displayBox(x,y,h,w,"B0","black")
displayCR(x-100,y-100,h,w,"B0","black")


# boite candidate -------------------------------------------------------------

# taille et positionnement aleatoire
h=0#rand(50:100)              # hauteur et largeur
w=0#rand(50:100)              
x=rand(50:CanvasSize-w-20)  # point inferieur gauche (ideal) de la boite 
y=rand(50:CanvasSize-h-20)

b1 = createBox(x,y,h,w)
displayCR(x-100,y-100,h,w,"B1","blue")
#displayBox(x,y,h,w,"B1","blue")


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


# epilogue --------------------------------------------------------------------
finish()
preview()
