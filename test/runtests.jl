using CoilFunctions
using CoilFunctions: ideal_fill, AWG_Chart
using Unitful
using Test

@testset "CoilGeometry" begin
    @test CoilGeometry(10u"mm",20u"mm",30u"mm") == CoilGeometry(10.0u"mm",20.0u"mm",30.0u"mm")
    @test CoilGeometry(10u"mm",20u"mm",30u"mm", 1.5u"mm") == CoilGeometry(13u"mm",17.0u"mm",27u"mm")
    @test CoilGeometry(10u"mm",20u"mm",30u"mm") ≈ CoilGeometry(prevfloat(10.0)u"mm",nextfloat(20.0)u"mm",30.0u"mm")
end

@testset "Ideal Fill" begin
    let
        fillfactor=0.9
        awg=AWG_Chart[30]
        a=CoilGeometry(10u"mm",10u"mm"+10u"mm",30.0u"mm")
        b=CoilGeometry(10u"mm",10u"mm"+10u"mm"*sqrt(fillfactor),30.0u"mm"*sqrt(fillfactor))
        @test ideal_fill(a,awg) == CoilFunctions.WindingLayers(19, 101, 101, 1010, 909, 1919)
        @test ideal_fill(a,awg,fillfactor) == CoilFunctions.WindingLayers(18, 96, 96, 864, 864, 1728)
        @test ideal_fill(a,awg,fillfactor) == ideal_fill(b,awg)

        b=CoilGeometry(10u"mm",10u"mm"+2.01awg.d_insul,10.01awg.d_insul)
        @test ideal_fill(b,awg) == WindingLayers(1, 10, 9, 10, 0, 10)
    end
end

@testset "Winding Associativity" begin
    let
        g = (layers, odd, id, gauge) -> begin
            f = (layers, odd, even, id, gauge) -> begin
                coil = enclosewinding(layers, odd, even, id, gauge)
                winding = ideal_fill(coil, gauge)
                winding.layers == layers && winding.odd_layer == odd && winding.even_layer == even
            end
            all(f(layers,odd,odd,id,gauge) for layers in layers for odd in odd for id in id) &&
            all(f(layers,odd,odd-1,id,gauge) for layers in layers for odd in odd for id in id)
        end
        @test all(x->g(1:100,1:200,range(10u"mm";stop=120u"mm",length=17),x),values(AWG_Chart))
    end
end

@testset "Layer Associativity" begin
    let
        f = (n,gauge) -> begin
            d=gauge.d_insul
            w=CoilFunctions.oddwidth(n,d) |> u"mm"
            _n =CoilFunctions.oddturns(w,d)
            # @info (n,_n,d,w)
            n == _n
        end
        @test all(f(n,gauge) for n in 1:1000 for (key,gauge) in CoilFunctions.AWG_Chart)
    end
    let
        f = (n,gauge) -> begin
            d=gauge.d_insul
            w=CoilFunctions.evenwidth(n,d) |> u"mm"
            _n =CoilFunctions.eventurns(w,d)
            # @info (n,_n,d,w)
            n == _n
        end
        @test all(f(n,gauge) for n in 1:1000 for (key,gauge) in CoilFunctions.AWG_Chart)
    end
    let
        f = (n,gauge) -> begin
            d=gauge.d_insul
            h=CoilFunctions.layerheight(n,d) |> u"mm"
            _n =CoilFunctions.ideal_layers(h,d)
            # @info (n=n,n2=_n,d=d,h=h)
            n == _n
        end
        @test all(f(n,gauge) for n in 1:100 for (key,gauge) in CoilFunctions.AWG_Chart)
    end
end

@testset "Basic Calculations" begin
    @test CoilFunctions.resistance(1000u"ft"/π,CoilFunctions.AWG_Chart[30]) ≈ 102.7u"Ω"
end

@testset "Assertion Errors" begin
    @test_throws AssertionError enclosewinding(10,10,8,10u"mm",0.2u"mm")
    #TODO enclosedwinding from coil
end
