print("Loading and compiling the graphic package (PyPlot)...")
using PyPlot         
println(" done!")


function setupGraphic(fname, data, paving)

    # Set the axes orthonormed ------------------------------------------------
    y_f1max = paving[1].y12
    y_f2max = paving[1].y21
    for i = 2:length(paving)
        if paving[i].y12[2] > y_f2max[2]
            y_f2max = paving[i].y12
        end
        if paving[i].y21[1] > y_f1max[1]
            y_f1max = paving[i].y21
        end        
    end
    vmax = max(y_f1max[1], y_f1max[2], y_f2max[1], y_f2max[2])

    # Initialize the graphic --------------------------------------------------
    figure("Bi-objective binary UFLP",figsize=(6,6)) # Create a new figure
    title("fname=" * fname * " | nI=" * string(data.nI) * " nJ=" * string(data.nJ))
#    title("fname=" * fname * " | nI=" * string(data.nI) * " nJ=" * string(data.nJ))

    xlabel(L"$f^1(x,s),<$")
    ylabel(L"$f^2(x,s),<$")

    xlim(0.0,vmax+0.15*vmax)
    ylim(0.0,vmax+0.15*vmax)
    return nothing
end


function displayBoxJ1(b, name, thecolor, thelegend, thelinestyle)

    if thelinestyle=="-"
      plot([b.y12[1],b.yN[1],b.y21[1],b.yI[1],b.y12[1]], [b.y12[2],b.yN[2],b.y21[2],b.yI[2],b.y12[2]], lw=0.5, c=thecolor, linestyle=thelinestyle)
    end
    
    if thelinestyle=="--"
      plot([b.y12[1],b.yN[1],b.y21[1],b.yI[1],b.y12[1]], [b.y12[2],b.yN[2],b.y21[2],b.yI[2],b.y12[2]], lw=0.5, c=thecolor, linestyle=thelinestyle)
    end
    
      
    scatter([b.y12[1],b.y21[1]], [b.y12[2],b.y21[2]],s=10,c=thecolor,marker="o")
    #scatter([b.y12[1],b.y21[1]], [b.y12[2],b.y21[2]],s=15,c="white",marker="D")    
    scatter([b.yN[1],b.yI[1]], [b.yN[2],b.yI[2]],s=10,c=thecolor,marker="o")
    scatter([b.yN[1],b.yI[1]], [b.yN[2],b.yI[2]],s=5,c="white",marker="o")    
    #scatter([b.yN[1]],[b.yN[2]],s=25,c="black",marker="s")
    #scatter([b.yN[1]],[b.yN[2]],s=15,c="white",marker="s")    
    #scatter([b.yI[1]],[b.yI[2]],s=45,c="black",marker="*")
    #scatter([b.yI[1]],[b.yI[2]],s=25,c="white",marker="*")    
    scatter([b.y12[1],b.y21[1]], [b.y12[2],b.y21[2]],s=10,c=thecolor,marker="o")    
    
    if thelegend == :legendOn
        #text( b.yN[1],  b.yN[2],  L"y^{N}",  ha="left",   va="bottom", size=12)
        #text( b.y12[1], b.y12[2], L"y^{12}", ha="center", va="bottom", size=12)
        #text( b.y21[1], b.y21[2], L"y^{21}", ha="left",   va="top", size=12)
        #text( b.yI[1],  b.yI[2],  L"y^{I}",  ha="right",  va="top", size=12)
        text( (b.y12[1]+b.y21[1])/2, (b.y12[2]+b.y21[2])/2 -15, L"\ \{"*replace(string(b.J1)[begin+1:end-1], " " => "")*"}", ha="center", va="top", size=8, color=thecolor)
    end

    return nothing
end

function displayBoxJ1bis(b, name, thecolor, thelegend, thelinestyle)

    if thelinestyle=="-"
      plot([b.y12[1],b.yN[1],b.y21[1],b.yI[1],b.y12[1]], [b.y12[2],b.yN[2],b.y21[2],b.yI[2],b.y12[2]], lw=0.5, c=thecolor, linestyle=thelinestyle)
    end
    
    if thelinestyle=="--"
      plot([b.y12[1],b.yN[1],b.y21[1],b.yI[1],b.y12[1]], [b.y12[2],b.yN[2],b.y21[2],b.yI[2],b.y12[2]], lw=0.5, c=thecolor, linestyle=thelinestyle)
    end
    
      
    scatter([b.y12[1],b.y21[1]], [b.y12[2],b.y21[2]],s=10,c=thecolor,marker="o")
    #scatter([b.y12[1],b.y21[1]], [b.y12[2],b.y21[2]],s=15,c="white",marker="D")    
    scatter([b.yN[1],b.yI[1]], [b.yN[2],b.yI[2]],s=10,c=thecolor,marker="o")
    scatter([b.yN[1],b.yI[1]], [b.yN[2],b.yI[2]],s=5,c="white",marker="o")    
    #scatter([b.yN[1]],[b.yN[2]],s=25,c="black",marker="s")
    #scatter([b.yN[1]],[b.yN[2]],s=15,c="white",marker="s")    
    #scatter([b.yI[1]],[b.yI[2]],s=45,c="black",marker="*")
    #scatter([b.yI[1]],[b.yI[2]],s=25,c="white",marker="*")    
    scatter([b.y12[1],b.y21[1]], [b.y12[2],b.y21[2]],s=10,c=thecolor,marker="o")    
    
    if thelegend == :legendOn
        #text( b.yN[1],  b.yN[2],  L"y^{N}",  ha="left",   va="bottom", size=12)
        #text( b.y12[1], b.y12[2], L"y^{12}", ha="center", va="bottom", size=12)
        #text( b.y21[1], b.y21[2], L"y^{21}", ha="left",   va="top", size=12)
        #text( b.yI[1],  b.yI[2],  L"y^{I}",  ha="right",  va="top", size=12)
        text( (b.y12[1]+b.y21[1])/2+10, (b.y12[2]+b.y21[2])/2, L"\ \{"*replace(string(b.J1)[begin+1:end-1], " " => "")*"}", ha="left", va="center", size=8, color=thecolor)
    end

    return nothing
end

function displayAllRunningCosts(paving)

    X=(Int64)[]
    Y=(Int64)[]
    for i in reverse(eachindex(paving))
        i <=length(paving)-5 && break
        #if #(i == length(paving)-1) || (i == length(paving)-3) || (i == length(paving)-5) || 
  #      if    ((i == length(paving)-9+1) || (i == length(paving)-14+1) || (i == length(paving)-7+1))
            push!(X,paving[i].CR[1]); push!(Y,paving[i].CR[2])
           displayCRplusCA12andCA21(paving[i])
   #     end
    end
   
    scatter(X, Y, s=10, c="red", marker="o")
    #text( X, Y,   L"y^{R}",  ha="left",  va="top", size=12)
    scatter(X, Y, s=5, c="white", marker="o")
    return nothing
end


function displayCRplusCA12andCA21(b::Box)

    thecolor = "red"; thelegend = :legendOn; thelinestyle = "-"
    displayBoxJ1bis(b, 0, thecolor, thelegend, thelinestyle)

    annotate("", xytext=(0,0), xy=(b.CR[1],b.CR[2]), arrowprops=Dict("arrowstyle"=>"->","linestyle"=>"--","color"=>"red"))
    annotate("", xytext=(b.CR[1],b.CR[2]), xy=(b.y12[1],b.y12[2]), arrowprops=Dict("arrowstyle"=>"->","linestyle"=>"--","color"=>"red"))
    annotate("", xytext=(b.CR[1],b.CR[2]), xy=(b.y21[1],b.y21[2]), arrowprops=Dict("arrowstyle"=>"->","linestyle"=>"--","color"=>"red"))
  
   # text( (b.y12[1]+b.y21[1])/2, (b.y12[2]+b.y21[2])/2 +5, L"\ \{"*replace(string(b.J1)[begin+1:end-1], " " => "")*"}", ha="center", va="bottom", size=8,color="red")

    return nothing
end


function displayBoxedyNS(b::Box)

    X=(Int64)[]
    Y=(Int64)[]
    for i in eachindex(b.listyNS)
        push!(X,b.listyNS[i][1]); push!(Y,b.listyNS[i][2])
    end

    scatter(X, Y, s=30, c="red", marker=".")
    return nothing
end


function displayYN(Y_N,thecolor,themarkerstyle,themarkersize)

    X=(Int64)[]
    Y=(Int64)[]
    for i in eachindex(Y_N)
        push!(X,Y_N[i][1]); push!(Y,Y_N[i][2])
    end
    scatter(X, Y, s=themarkersize, c=thecolor, marker=themarkerstyle)
    #text( X[2], Y[2],  L"\ y \in Y_{N}",  ha="left",  va="bottom", size=10)
    return nothing
end

function displayPaving(paving)

    thecolor = "grey"; thelegend = :legendOn; thelinestyle = "-"
    for i in reverse(eachindex(paving))
        #i <=length(paving)-15 && break
       # if (i == length(paving)-1+1) || (i == length(paving)-3+1) || (i == length(paving)-5+1) # || ((i == length(paving)-9+1) || (i == length(paving)-14+1)|| (i == length(paving)-7+1))
            displayBoxJ1(paving[i], i, thecolor, thelegend, thelinestyle)
       #end
        
    end
    return nothing
end

#=
setupGraphic(fname, data, paving)

ib=1
thecolor = "red"; thelegend = :legendOn; thelinestyle = "-"
displayBoxJ1(paving[ib], ib, thecolor, thelegend, thelinestyle)



ib=2
thecolor = "black"; thelegend = :legendOn; thelinestyle = "--"
displayBoxJ1(paving[ib], ib, thecolor, thelegend, thelinestyle)

if draw_allRunningCosts
    displayAllRunningCosts(paving)
end

if draw_CRplusCA12andCA21
    ib=3
    displayCRplusCA12andCA21(paving[ib])
end

if draw_yNS
    ib=3
    displayBoxedyNS(paving[ib])
end


if draw_nondominatedPointsALGO
    thecolor = "orange"; themarkerstyle = "*"; themarkersize = 15
    displayYN(all_YN,thecolor,themarkerstyle,themarkersize)
    thecolor = "green"; themarkerstyle = "o"; themarkersize = 15
    displayYN(ND_YN,thecolor,themarkerstyle,themarkersize)

end

if draw_nondominatedPoints
    thecolor = "yellow"; themarkerstyle = "+"; themarkersize = 25
    displayYN(Y_N,thecolor,themarkerstyle,themarkersize)
end

draw_paving = true
if draw_paving
    displayPaving(paving)
end
=#