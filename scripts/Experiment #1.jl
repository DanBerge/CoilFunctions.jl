# id = range(80u"mm";stop=140u"mm",step=0.1u"mm")
# od = range(81u"mm";stop=141u"mm",step=0.1u"mm")
id = coil.id:0.1mm:coil.od
od = id

const awg=CoilFunctions.AWG_Chart
resistance = x->x.Ω
turns = x->x.turns
ampturns=x->x.ampturns
a=let #awg=awg[23]
    f = (id,od,len,awg) -> begin
        coil = CoilGeometry(id,od,len)
        y=CoilFunctions.estimatetruefill(coil,awg)
        V=104u"V"
        ampturns=V/resistance(y)*turns(y) |> u"A"
        merge(y,(awg=awg,ampturns=ampturns))
    end
    [f(id,od,len,awg[i]) for id in id for od in od for len in 20u"mm":0.1u"mm":40u"mm" for i in 31:34 if id<od]
end
b=filter(x->isapprox(resistance(x),360.5u"Ω";rtol=0.01),a)
b=filter(x->isapprox(2574,turns(x);rtol=0.01),b)

let
    id = x->x.coil.id |> ustrip
    od = x->x.coil.od |> ustrip
    turns = x->x.turns
    plt=plot(xlabel="ID (mm)",ylabel="OD (mm)")
    for g in unique(gauge,b)
        y=filter(x->gauge(x)==gauge(g),b)
        sort!(y;by=id)
        plot!(id,od,y,label="$(gauge(g)) AWG")
    end
    plt
end

let
    id = x->x.coil.id |> ustrip
    od = x->x.coil.od |> ustrip
    turns = x->x.turns
    plt1=scatter()

end

let
    a=0:0.1:2
    f = (a,b) ->isapprox(a,b,atol(0.11))
    isequal = (a,b) -> isapprox(a,b,atol=0.11)
    unique(a)
end

function approxunique(a;rtol=0.001)
    y=[first(a)]
    f = (a,b) -> isapprox(a,b,rtol=rtol)
    id = x->x.coil.id
    od = x->x.coil.od
    turns = x->x.turns
    for a in a
        if  all(y) do x
                !f(id(a),id(x)) && !f(od(a),od(x)) && !f(turns(a),turns(x))
            end
            push!(y,a)
        end
    end
    y
end
