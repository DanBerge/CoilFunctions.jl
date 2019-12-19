using Plots, ThreadTools, Transducers

let s=:SBE5
    rng=(12:1:204).*u"V"
    coil=CoilFunctions.coil[s]
    p=CoilFunctions.power[s]
    f = v -> let
        Ω=CoilFunctions.Ω(s,v)
        turns=optimalcoil(coil,Ω,0.9).turns
        (v=v,turns=turns*24u"V"/v |> upreferred)
    end
    plot(;xlabel="Voltage (V)",ylabel="Relative Turns",title="Relative Turns $s")
    plot!(x->ustrip(x.v),y->y.turns,f.(rng);legend=:none)
end

let s=:SBE9,reqturns=520
    rng=(12:1:204).*u"V"
    coil=CoilFunctions.coil[s]
    p=CoilFunctions.power[s]
    f = v -> let
        Ω=CoilFunctions.Ω(s,v)
        minturns=reqturns*v/24u"V" |> ceil |> Int
        y = CoilFunctions.minimumcoillength(coil,Ω,minturns,0.9)
        (v=v,len=y.coil.len)
    end
    # g = minturns -> let
    #     maximum(x->x.len, f.(rng) )
    # end
    plot(;xlabel="Voltage (V)",ylabel="Relative Turns",title="Relative Turns $s")
    plot(x->ustrip(x.v),y->ustrip(y.len),f.(rng))
end

let s=:SBE7,reqturns=520
    vrng=(12:1:204).*u"V"
    coil=CoilFunctions.coil[s]
    p=CoilFunctions.power[s]
    f = (v,minturns) -> let
        Ω=v^2/p |> u"Ω"
        relative_turns=minturns*v/24u"V" |> ceil |> Int
        y = CoilFunctions.minimumcoillength(coil,Ω,relative_turns,0.9)
        (v=v,len=y.coil.len)
    end
    g = minturns -> let
        (minturns=minturns,len=maximum(x->x.len, f.(vrng,minturns)),)
    end
    # g(500)
    y=tmap(g,400:10:600)
    plot(;xlabel="Turns",ylabel="Minimum Coil Length (mm)",title="Relative Coil Length $s")
    plot!(x->x.minturns,y->ustrip(y.len),y)
end

plt=let
    vrng=(12:1:204).*u"V"
    plt=plot(;xlabel="Turns",ylabel="Minimum Coil Length (mm)",title="Relative Coil Length",
        legend=:topleft,ylims=(15,35))
    for s in (:SBE4,:SBE5,:SBE7,:SBE9)
        g = minturns -> let
            coil=CoilFunctions.coil[s]
            p=CoilFunctions.power[s]
            f = (v,minturns) -> let
                Ω=v^2/p |> u"Ω"
                relative_turns=minturns*v/24u"V" |> ceil |> Int
                y = CoilFunctions.minimumcoillength(coil,Ω,relative_turns,0.9)
                (v=v,len=y.coil.len)
            end
            (minturns=minturns,len=maximum(x->x.len, f.(vrng,minturns)),)
        end
        y=tmap(g,400:10:500)
        plot!(x->x.minturns,y->ustrip(y.len),y;label="$s")
    end
    plt
end

let s=:SBE4
    coil=CoilFunctions.coil[s]
    p=CoilFunctions.power[s]
    vrng = 12u"V":2u"V":208u"V"
    f = v -> let
        Ω = v^2/p |> u"Ω"
        merge((v=v,),CoilFunctions.optimalgeometry(coil,Ω,0.9))
    end
    y=tmap(f,vrng)
    plot(x->ustrip(x.v),y->ustrip(y.len-26u"mm"),y;label="Coil Length Variation")
    plot!(x->ustrip(x.v),y->ustrip(y.insul),y;label="Insulation")
end

let s=:SBE4
    coil=CoilGeometry(CoilFunctions.coil[s];len=24u"mm")
    p=CoilFunctions.power[s]
    vrng = 12u"V":2u"V":208u"V"
    f = v -> let
        Ω = v^2/p |> u"Ω"
        y=CoilFunctions.optimalcoil(coil,Ω,0.9)
        (v=v,insul=(y.coil.od-coil.od)/-2)
    end
    y=tmap(f,vrng)
    plot!(x->ustrip(x.v),y->ustrip(y.insul),y)
end

g=(s) -> let
    p=CoilFunctions.power[s]
    # vrng = 12u"V":1u"V":208u"V"
    vrng = 140u"V":0.1u"V":160u"V"
    rng = 26u"mm" .+ (-3u"mm":0.01u"mm":3u"mm")
    coil=CoilFunctions.coil[s]
    f = (v,len) -> let
        newcoil=CoilGeometry(coil;len=len)
        Ω = v^2/p |> u"Ω"
        y=CoilFunctions.optimalcoil(newcoil,Ω,0.9)
        (v=v,insul=(coil.od-y.coil.od)/2,len=len)
    end
    a=((v,len) for v in vrng, len in rng)
    g = x -> f(x[1],x[2])
    x=tmap(g,a)
    y=minimum(x->x.insul,x;dims=2)
    plot(ustrip.(vrng),ustrip.(y))
end



coil,p  = CoilGeometry(62.5mm,87.5mm,25mm,2mm), 30W
coil,p  = CoilGeometry(82.5mm,110mm,30mm,2mm), 44W
coil,p = CoilGeometry(100mm,140mm,30mm,2mm), 50W
coil,p  = CoilGeometry(140mm,190mm,30mm,2mm), 86W

map([12,24,48,104,180,207].*u"V") do x
    Ω=(x)^2/p |> u"Ω"
    I=(x)/Ω |> u"A"
    opt=CoilFunctions.optimalcoil(coil,Ω,0.9)
    # @info opt
    (V=x,Ω=Ω,gauge=Int(opt.gauge),turns=opt.turns,ampturns=opt.turns*I)
end

a=map((10:1:100).*W) do P
    P
    y=map([24,104]) do x
        Ω=(x*V)^2/P |> u"Ω"
        I=(x*V)/Ω |> u"A"
        opt=CoilFunctions.optimalcoil(coil,Ω)
        (V=x*V,Ω=Ω,gauge=opt.gauge,turns=opt.turns,ampturns=opt.turns*I)
    end
    ustrip(P),y[1].ampturns / y[2].ampturns
end

Plots.plot(x->first(x),x->last(x),a;ylabel="Ratio",xlabel="Power (Watts)")

b=let#
    h=30mm:0.1mm:60mm
    AWG=filter(x->first(x)==32,CoilFunctions.AWG_Chart)
    Ω=(104V)^2/30W |> u"Ω";
    # I=104V /
    a=map(h) do x
        a=CoilGeometry(coil;len=x)
        opt=CoilFunctions.optimalcoil(a,Ω;AWG_Chart=AWG)
    end
    Plots.plot(ustrip.(h),map(x->x.turns,a))

    h=20mm:0.1mm:60mm
    a=map(h) do x
        a=CoilGeometry(coil;len=x)
        opt=CoilFunctions.optimalcoil(a,Ω)
    end
    Plots.plot!(ustrip.(h),map(x->x.turns,a))
end

a = let
    id=100.0u"mm":1u"mm":105u"mm"
    od=140.0u"mm"
    v=24.0u"V"
    p=50.0u"W"
    Ω = v^2/p |> u"Ω"
    rng=(20:0.1:30)*u"mm"
    awg_chart = values(AWG_Chart)
    f = (id,len,gauge) -> begin
        pad=1.5u"mm"
        coil=CoilGeometry(id+2pad,od-2pad,len)
        CoilFunctions.estimatetruefill(coil,gauge,Ω,0.9)
    end
    y=((id,len,gauge) for len in rng for id in id for gauge in awg_chart) |> collect
    # map(x->f(x...),y)
    m = Map(x->f(x...)) |> Filter(x-> Ω/1.1 <= x.Ω <= Ω/0.9)
    tcollect(m,y)

end
