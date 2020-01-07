
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

"""
Calculate the minimum length coil necessary to achieve the minimum turns of wire

    minimumcoillength(coil, Ω, minturns, fillfactor=1.0)
    minimumcoillength(coil, v, p, minturns, fillfactor=1.0)

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

minimumcoillength(coil::CoilGeometry, v::Unitful.Voltage, p::Unitful.Power,
    minturns::Integer, fillfactor::Real=1.0) = minimumcoillength(coil, u"Ω"(v^2/p), minturns, fillfactor)

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

"""

Calculate the coil wire length and weight from an `optimalcoil` calculation.

"""
function weight(a)
    coil=a.coil
    turns=a.turns
    Ω=a.Ω
    gauge=AWG_Chart[a.gauge]
    feet=Ω/gauge.Ω
    return (length=feet,weight=feet/gauge.ft_lb)
end
