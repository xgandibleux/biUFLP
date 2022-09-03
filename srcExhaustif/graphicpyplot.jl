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
    figure("Bi-objective discrete UFLP",figsize=(6,6)) # Create a new figure
    title("fname=" * fname * " | nI=" * string(data.nI) * " nJ=" * string(data.nJ))

    xlabel(L"$f^1(x,s),<$")
    ylabel(L"$f^2(x,s),<$")

    xlim(0.0,vmax+0)
    ylim(0.0,vmax+0)
    return nothing
end


function displayBoxJ1(b, name, thecolor, thelegend, thelinestyle)

    if thelinestyle=="-"
      plot([b.y12[1],b.yN[1],b.y21[1],b.yI[1],b.y12[1]], [b.y12[2],b.yN[2],b.y21[2],b.yI[2],b.y12[2]], lw=0.5, c=thecolor, linestyle=thelinestyle)
    end
    
    if thelinestyle=="--"
      plot([b.y12[1],b.yN[1],b.y21[1],b.yI[1],b.y12[1]], [b.y12[2],b.yN[2],b.y21[2],b.yI[2],b.y12[2]], lw=0.5, c=thecolor, linestyle=thelinestyle)
    end
    
      
    scatter([b.y12[1],b.y21[1]], [b.y12[2],b.y21[2]],s=2,c="black",marker="o")
    scatter([b.yN[1],b.yI[1]], [b.yN[2],b.yI[2]],s=2,c="black",marker="D")
    
    if thelegend == :legendOn
        text( b.yN[1],  b.yN[2],  L"z^{N}",  ha="left",   va="bottom")
        text( b.y12[1], b.y12[2], L"z^{12}", ha="center", va="bottom")
        text( b.y21[1], b.y21[2], L"z^{21}", ha="left",   va="top")
        text( b.yI[1],  b.yI[2],  L"z^{I}",  ha="right",  va="top")
        text( (b.y12[1]+b.y21[1])/2, (b.y12[2]+b.y21[2])/2, L"B"*string(name), ha="center", va="center")
    end

    return nothing
end

function displayAllRunningCosts(paving)

    X=(Int64)[]
    Y=(Int64)[]
    for i in eachindex(paving)
        push!(X,paving[i].CR[1]); push!(Y,paving[i].CR[2])
    end

    scatter(X, Y, s=2, c="purple", marker="^")
    return nothing
end


function displayCRplusCA12andCA21(b::Box)

    annotate("", xytext=(0,0), xy=(b.CR[1],b.CR[2]), arrowprops=Dict("arrowstyle"=>"->","linestyle"=>"--","color"=>"red"))
    annotate("", xytext=(b.CR[1],b.CR[2]), xy=(b.y12[1],b.y12[2]), arrowprops=Dict("arrowstyle"=>"->","linestyle"=>"--","color"=>"red"))
    annotate("", xytext=(b.CR[1],b.CR[2]), xy=(b.y21[1],b.y21[2]), arrowprops=Dict("arrowstyle"=>"->","linestyle"=>"--","color"=>"red"))
  
    thecolor = "red"; thelegend = :legendOn; thelinestyle = "-"
    displayBoxJ1(paving[ib], ib, thecolor, thelegend, thelinestyle)
    return nothing
end


function displayBoxedyNS(b::Box)

    X=(Int64)[]
    Y=(Int64)[]
    for i in eachindex(b.listyNS)
        push!(X,b.listyNS[i][1]); push!(Y,b.listyNS[i][2])
    end

    scatter(X, Y, s=2, c="red", marker=".")
    return nothing
end

function displayYN(Y_N,thecolor,themarkerstyle,themarkersize)

    X=(Int64)[]
    Y=(Int64)[]
    for i in eachindex(Y_N)
        push!(X,Y_N[i][1]); push!(Y,Y_N[i][2])
    end
    scatter(X, Y, s=themarkersize, c=thecolor, marker=themarkerstyle)
    return nothing
end

function displayPaving(paving)

    thecolor = "grey"; thelegend = :legendOff; thelinestyle = "-"
    for ib in eachindex(paving)
        displayBoxJ1(paving[ib], ib, thecolor, thelegend, thelinestyle)
    end
    return nothing
end

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