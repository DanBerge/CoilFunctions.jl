# module Drawing

function drawwire(p::Point, awg::AWG, fillfactor=1.0;u=u"mm")
    f = x->ustrip(u(x))
    setcolor("green")
    circle(p,f(awg.d_insul)/2,:fill)
    setcolor("blue")
    circle(p,f(awg.d_bare)/2,:fill)
    # setcolor("black")
    # circle(p,f(awg.d_insul)/sqrt(fillfactor)/2, :stroke)
    # @info f(awg.d_insul)
end

function drawcoil(awg::AWG,winding::WindingLayers,maxturns::Integer,fillfactor::Real=1.0)
    turns = zero(maxturns)
    d = awg.d_insul/sqrt(fillfactor) |> u"mm" |> ustrip
    for (layer,y) = enumerate(range(d/2;step=d*√(3)/2,length=winding.layers))
        (turnlayers,xoffset) = isodd(layer) ? (winding.odd_layer,d/2) : (winding.even_layer,d)
        for x in range(xoffset;step=d,length=turnlayers)
            turns >= maxturns && break
            drawwire(Point(x,y),awg,fillfactor)
            turns += 1
        end
    end
end
drawcoil(a,fillfactor) = drawcoil(CoilFunctions.AWG_Chart[a.gauge], a.winding, a.turns, fillfactor)

function drawcoil(coil::CoilGeometry)
    y=(coil.od - coil.id)/2 |> u"mm" |> ustrip
    x=ustrip(u"mm"(coil.len))
    gsave()
    setcolor("black")
    box(Point(0,0),Point(x,y), :stroke)
    grestore()
end

function foo()
    let x=800,y=300,pad=20
        # Drawing(x+2pad, y+2pad, "luxor-drawing-(timestamp).svg")
        Drawing(x+2pad, y+2pad, :png)
        origin(pad,y+pad)
        s = ustrip(x/coil.len)
        scale(s,-s)
        background("white")
        sethue("black")
        drawcoil(coil)
        drawcoil(a,a.fillfactor)
        # drawcoil(merge(a,(turns=100,)),a.fillfactor)
        finish()
        preview()
    end
end

# end
