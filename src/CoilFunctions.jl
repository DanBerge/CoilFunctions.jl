module CoilFunctions

# import CSV
using Unitful, DataStructures, MappedArrays
using ThreadTools: tmap
using Transducers
import Unitful: 𝐈,𝐋,𝐌,𝐓
using Roots: find_zero

export ideal_fill, AWG_Chart, optimalcoil, estimatetruefill, enclosewinding, resistance

@derived_dimension ResistanceLength dimension(u"Ω/m") #𝐈^-2*𝐋*𝐌*𝐓^-3

include("types.jl")
include("defaults.jl")
include("winding.jl")
include("coil.jl")
include("optimization.jl")

include("drawing.jl")
# using .Drawing

const AWG_Chart=let
    src=readlines(joinpath((@__DIR__),"AWG Chart.csv"))
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


function istrianglesorted(a::AbstractVector)
    i = argmax(a)
    issorted(@view a[1:i]) && issorted(@view a[i:end];lt= > )
end

end # module
