# baremodule Defaults
# using CoilFunctions, Unitful
# # using Unitful
#
# # abstract type Product end
# # struct SBE4 <: Product end
# # struct SBE5 <: Product end
# # struct SBE7 <: Product end
# # struct SBE9 <: Product end
# #
# # CoilGeometry(::Type{SBE4}) = CoilGeometry(66.5u"mm",83.5u"mm",21u"mm")
# # CoilGeometry(::Type{SBE5}) = CoilGeometry(86.5u"mm",106u"mm",26u"mm")
# # CoilGeometry(::Type{SBE7}) = CoilGeometry(104u"mm",136u"mm",26u"mm")
# # CoilGeometry(::Type{SBE9}) = CoilGeometry(144u"mm",186u"mm",26u"mm")
# # CoilGeometry(::SBE5) = CoilGeometry(86.5u"mm",106u"mm",26u"mm")
# #
# # power(::Type{SBE4})=30.0u"W"
# # power(::Type{SBE5})=44.0u"W"
# # power(::Type{SBE7})=50.0u"W"
# # power(::Type{SBE9})=86.0u"W"
#
# # coil=Dict(:SBE4 => CoilGeometry(66.5u"mm",83.5u"mm",21u"mm"),
# #     :SBE5 => CoilGeometry(86.5u"mm",106u"mm",26u"mm"),
# #     :SBE7 => CoilGeometry(104u"mm",136u"mm",26u"mm"),
# #     :SBE9 => CoilGeometry(144u"mm",186u"mm",26u"mm"))
#
# end

const coil=Dict(:SBE4 => CoilGeometry(66.5u"mm",83.5u"mm",21u"mm"),
    :SBE5 => CoilGeometry(86.5u"mm",106u"mm",26u"mm"),
    :SBE7 => CoilGeometry(104.0u"mm",136u"mm",26u"mm"),
    :SBE9 => CoilGeometry(144.0u"mm",186u"mm",26u"mm"))

const coilcavity=Dict(:SBE4 => CoilGeometry(62.5u"mm",83.5u"mm",25u"mm"),
    :SBE5 => CoilGeometry(82.5u"mm",110u"mm",30u"mm"),
    :SBE7 => CoilGeometry(100.0u"mm",140u"mm",30u"mm"),
    :SBE9 => CoilGeometry(140.0u"mm",190u"mm",30u"mm"))


const power=Dict(:SBE4 => 30.0u"W", :SBE5 => 44.0u"W", :SBE7=>50.0u"W",:SBE9=>86.0u"W")
Ω(s::Symbol,v::Unitful.Voltage) = v^2/power[s] |> u"Ω"
