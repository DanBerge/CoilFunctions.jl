module Drawing

using CoilFunctions, Unitful, Luxor

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

function illustratecoil(coil::CoilGeometry, Ω::Unitful.ElectricalResistance, fillfactor=1.0; filename=:png, pad=20,
    resolution=(800,600), minturns=0, maxturns=typemax(Int64))

    x,y = resolution
    a=optimalcoil(coil,Ω,fillfactor)
    Drawing(x+2pad, y+2pad, filename)
    origin(pad,y+pad)
    s = ustrip(x/coil.len)
    scale(s,-s)
    background("white")
    sethue("black")
    drawcoil(coil)
    turns = clamp(a.turns,minturns, maxturns)
    drawcoil(merge(a,(turns=turns,)),a.fillfactor)
    finish()
    @info a
    preview()
end

function illustratecoil(coil::CoilGeometry, v::Unitful.Voltage, p::Unitful.Power, fillfactor=1.0; kwargs...)
    illustratecoil(coil, v^2/p |> u"Ω", fillfactor;kwargs...)
end

function illustratecoil(s::Symbol, v::Unitful.Voltage, fillfactor=1.0; kwargs...)
    coil=CoilFunctions.coil[s]
    Ω=CoilFunctions.Ω(s,v)
    illustratecoil(coil,Ω,fillfactor;kwargs...)
end

end
