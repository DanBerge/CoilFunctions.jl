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
