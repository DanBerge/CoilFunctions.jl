
"""
isvalidcoil(a::CoilGeometry)

Test whether the coil has valid geometry.
"""
isvalidcoil(a::CoilGeometry) = a.od > a.id

"""
Calculate the maximum number of turns of wire for a given `CoilGeometry`, coil resistance `Ω`, and `fillfactor`.
Optional keyword `AWG_Chart::AbstractDict{AWG}` to specify wire options.

    optimalcoil(coil::CoilGeometry,Ω::Unitful.ElectricalResistance,fillfactor::Real=1.0; AWG_Chart=AWG_Chart

"""
function optimalcoil(coil::CoilGeometry,Ω::Unitful.ElectricalResistance,fillfactor::Real=1.0; AWG_Chart=AWG_Chart)
    @assert 0 < fillfactor <= 1
    y=[(gauge=gauge,estimatetruefill(coil,awg,Ω,fillfactor)...,fillfactor=fillfactor) for (gauge,awg) in AWG_Chart]
    filter!(x->x.Ω >= Ω,y)
    ind=last(findmax(map(x->x.turns,y)))
    y[ind]
end

estimatetruefill(coil::CoilGeometry, awg::AWG, Ω::Unitful.ElectricalResistance,fillfactor::Number=1.0;
    maxturns::Number=Inf) = estimatetruefill(coil,awg;Ω=Ω, fillfactor=fillfactor, maxturns=maxturns)

"""
Estimate the resulting winding given a `coil` and `gauge` of wire, and `fillfactor`.
Optionally limit resistance `Ω`, `maxturns`

    estimatetruefill(coil::CoilGeometry, awg::AWG; Ω::Unitful.ElectricalResistance=Inf*u"Ω",fillfactor::Number=1.0, maxturns::Number=Inf)
"""
function estimatetruefill(coil::CoilGeometry, awg::AWG; Ω::Unitful.ElectricalResistance=Inf*u"Ω",
    fillfactor::Number=1.0, maxturns::Number=Inf)

    f=ideal_fill(coil,awg,fillfactor)
    if f.turns==0
        # return (turns=0,Ω=zero(Ω),winding=f,coil=coil)
    end
    maxturns=min(maxturns,f.turns) |> Int
    wire_diam = awg.d_insul / sqrt(fillfactor)
    d0 = coil.id + wire_diam
    dstep = √3*wire_diam
    state = (layers=0, turns=0, Ω=0.0u"Ω")
    for (layer,d) in Iterators.take(enumerate(Iterators.countfrom(d0, dstep)), f.layers)
        layerturns = isodd(layer) ? f.odd_layer : f.even_layer
        layerturns = Int(min(maxturns-state.turns,layerturns)) # turns due to turn limit)

        turnΩ = awg.Ω*π*d |> unit(Ω)
        Ωturns = isinf(Ω) ? layerturns : (Ω-state.Ω)/turnΩ |> ceil |> Int # turns due to resistance limit

        layerturns = min(layerturns, Ωturns)

        #update State
        state = (layers=layer,
            turns=state.turns+layerturns,
            Ω=state.Ω+layerturns*turnΩ)

        if (state.turns >= maxturns) || (state.Ω ≳ Ω)
            break
        end
    end
    dfinal = coil.id+begin
        if state.layers <= 1
            2wire_diam
        else
            2wire_diam + sqrt(3)*wire_diam*(state.layers-1)
        end
    end
    coil_section=(coil.od-coil.id)/2*coil.len
    relativefill= (state.turns * π*(wire_diam/2)^2)/coil_section |> upreferred
    truefill= (state.turns * π*(awg.d_insul/2)^2)/coil_section |> upreferred
    (turns=state.turns,Ω=state.Ω,winding=f,
        coil=CoilGeometry(coil.id,unit(coil.od)(dfinal),coil.len),
        relativefill=relativefill,truefill=truefill,
        coil_section=coil_section,
        wire_cross=pi*(wire_diam/2)^2,)
end

function optimalgeometry(coil::CoilGeometry, Ω::Unitful.ElectricalResistance, fillfactor::Real=1.0;
    rng=0u"mm":0.01u"mm":2u"mm")

    f = insul -> let
        newcoil = CoilGeometry(coil; od=coil.od-2insul)
        len=gaugebreakpoint(newcoil,Ω,newcoil.len,fillfactor)
        err = 2abs(coil.len-len) + insul
        (insul=insul,len=len,error=err)
    end
    y=tmap(f,rng)
    y[last(findmin(map(x->x.error,y)))]
end
