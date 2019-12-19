module CoilFunctions

# import CSV
using Unitful, DataStructures, MappedArrays
using ThreadTools: tmap
using Transducers
import Unitful: 𝐈,𝐋,𝐌,𝐓
using Roots: find_zero

export ideal_fill, AWG_Chart, optimalcoil, estimatetruefill

@derived_dimension ResistanceLength dimension(u"Ω/m") #𝐈^-2*𝐋*𝐌*𝐓^-3

include("types.jl")
include("defaults.jl")

const AWG_Chart=let
    src=readlines((@__DIR__) * "\\AWG Chart.csv")
    data=map(src[2:end]) do x
        x=parse.(Float64,split(x,','))
        AWG(x[1],x[2]u"inch",x[5]u"inch",x[3]*u"Ω/1000ft",x[7]u"ft/lb")
    end
    AWG_Chart=Dict(x.gauge=>x for x in data)
    AWG_Chart=OrderedDict(sort(AWG_Chart, by=first))
    #Filter out half gauages
    filter!(x->isinteger(x[1]),AWG_Chart)
end

const _default_awg=AWG_Chart[30.0]

≲(a,b) = (a ≈ b) || (a < b)
≳(a,b) = (a ≈ b) || (a > b)

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

"""
Calculate the orthocyclic coil winding for a give wire diameter, returning a `WindingLayers`.

ideal_fill(w, h, d)

    w => Horizontal width of coil for layers
    h => Vertical height to fill with layers
    d => Diameter of Wire

    ←           w             →
    ###########################
    #           ↑             #
    #                         #
    #           h             #
    #                         #
    #           ↓             #
    ###########################
                ↑
              (id)
"""
function ideal_fill(w,h,d)

    # layers=length(d/2:d*sqrt(3)/2:h-d/2)
    #TODO Improve this to the rounding isn't necessary
    layers = ideal_layers(h,d)
    odd_layer= Int(w ÷ d)
    even_layer= Int((w-d/2) ÷ d)

    odd_turns = length(1:2:layers)*odd_layer
    even_turns = length(2:2:layers)*even_layer

    WindingLayers(layers, odd_layer, even_layer, odd_turns, even_turns, odd_turns+even_turns)
end

layers(h,d) = h>d ? Int(1+(h-d)÷(d*√3/2)) : 0
eventurns(w,d) = Int(w÷d)
oddturns(w,d) = Int((w-d/2)÷d)

"""
Calculate the number of layers of wire.  Given the `h` height of the layter, and the `d` diameter of the wire.

    ideal_layers(h, d)

    ###########################
    #           ↑             #
    #                         #
    #           h             #
    #                         #
    #           ↓             #
    ###########################
                ↑
              (id)
"""
function ideal_layers(h, d)
    layers = h>d ? 1 + (h-d)/(d*√3/2) : 0
    if round(layers) ≈ layers
        return Int(round(layers))
    else
        Int(floor(layers))
    end
end


"""
ideal_fill(coil::CoilGeometry, awg::AWG, fill::Real=1.0)

    coil => CoilGeometry
    awg => AWG of the wire
    fillfactor => fillfactor of the wire
"""
function ideal_fill(coil::CoilGeometry, awg::AWG, fillfactor::T=1.0) where {T<:Real}
    zero(T) < fillfactor <= one(T) || throw(DomainError(fillfactor))
    w=coil.len * sqrt(fillfactor)
    h=(coil.od/2-coil.id/2)*sqrt(fillfactor)
    ideal_fill(w,h,awg.d_insul)
end

"""
isvalidcoil(a::CoilGeometry)

Test whether the coil has valid geometry.
"""
isvalidcoil(a::CoilGeometry) = a.od > a.id

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

"""

    minimumcoillength(coil, Ω, minturns, fillfactor=1.0)

Calculate the minimum length coil necessary to achieve the minimum turns of wire


"""
function minimumcoillength(coil::CoilGeometry, Ω::Unitful.ElectricalResistance,
    minturns::Integer, fillfactor::Real=1.0)

    f=let Ω=u"Ω"(Ω)  #Simplify for type inference, necessary on <=1.3?
        x -> begin
            _coil = CoilGeometry(coil.id,coil.od,x)
            # _Ω = upreferred(Ω)
            optimalcoil(_coil, Ω, fillfactor)
        end
    end

    rng=5u"mm":0.01u"mm":300u"mm"
    m=mappedarray(x->f(x).turns,rng)
    i=searchsortedfirst(m,minturns)
    f(rng[i])
end

"""

    gaugebreakpoint(coil::CoilGeometry, Ω::Unitful.ElectricalResistance, len::Unitful.Length, fillfactor::Real=1.0)

Calculate length of coil with the nearest wire gauge break point, such that the geometry has the highest fill factor.

"""
function gaugebreakpoint(coil::CoilGeometry, Ω::Unitful.ElectricalResistance,
    len::Unitful.Length, fillfactor::Real=1.0)

    f=let Ω=u"Ω"(Ω)  #Simplify for type inference, necessary on <=1.3?
        x -> begin
            _coil = CoilGeometry(coil.id,coil.od,x)
            # _Ω = upreferred(Ω)
            optimalcoil(_coil, Ω, fillfactor)
        end
    end
    gauge = f(len).gauge

    rng=5u"mm":0.01u"mm":300u"mm"
    m=mappedarray(x->f(x).gauge,rng)
    a=rng[searchsortedfirst(m,gauge;rev=true)]
    b=rng[searchsortedlast(m,gauge;rev=true)]
    abs(len-a) < abs(len-b) ? a : b
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



function crossprocess_wind(coil::CoilGeometry, Ω::Unitful.ElectricalResistance,
    processfill::Number,targetfill::Number)

    targetwind=optimalcoil(coil,Ω,targetfill)
    turns  =targetwind.turns
    gauge  = AWG_Chart[targetwind.gauge]
    targetΩ=targetwind.Ω

    f=id -> begin
        newcoil = CoilGeometry(id, 2coil.od,coil.len)
        estimatetruefill(newcoil,gauge;Ω=targetΩ,maxturns=turns,fillfactor=processfill).Ω-targetΩ
    end
    id=find_zero(f,coil.id)
    y=estimatetruefill(CoilGeometry(coil;id=id),gauge;Ω=targetΩ,maxturns=turns,fillfactor=processfill)
    merge((gauge=targetwind.gauge,),y)
end

function minimumturns(coil::CoilGeometry, power::Unitful.Power, fillfactor::Real=1.0; op=minimum,
    vref::Unitful.Voltage=24u"V", rng::AbstractVector{<:Unitful.Voltage} = 12.0u"V":1u"V":208u"V")

    tmap(rng) do v
        Ω = v^2/power |> u"Ω"
        y = optimalcoil(coil, Ω, fillfactor)
        # @info (v=v,turns=y.turns,relative = y.turns * vref/v)
        y.turns * vref/v |> floor |> Int # Equivalent number of turns
    end |> op

    # g = v -> begin
    #     Ω = v^2/power |> u"Ω"
    #     y = optimalcoil(coil, Ω, fillfactor)
    #     # @info (v=v,turns=y.turns,relative = y.turns * vref/v, id=Threads.threadid())
    #     y.turns * vref/v |> floor |> Int
    # end
    # reduce(op, Map(g), rng)
end

function enclosewinding(layers, oddlayer, evenlayer, id::Unitful.Length, d::Unitful.Length)
    od = id + 2d + (layers-1)*d*√3
    len = oddlayers > evenlayer ? d*oddlayer : d/2 + d*evenlayer
    f = unit(id)
    CoilGeometry(id,f(od),f(len))
end

function enclosewinding(coil::CoilGeometry, winding::WindingLayers, d::Unitful.Length)
    enclosewinding(winding.layers, winding.odd_layers, winding.even_layers, coil.od, d)
end
enclosewinding(coil::CoilGeometry, winding::WindingLayers, gauge::AWG) = enclosewinding(coil,winding,gauge.d_insul)

function istrianglesorted(a::AbstractVector)
    i = argmax(a)
    issorted(@view a[1:i]) && issorted(@view a[i:end];lt= > )
end

end # module
