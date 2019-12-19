export AWG, CoilGeometry, WindingLayers

import Base: ==, ≈
const Length = Unitful.Length

struct AWG{T<:Real, L<:Unitful.Length, R<:ResistanceLength, Q<:Unitful.Quantity}
    gauge::T
    d_bare::L
    d_insul::L
    Ω::R
    ft_lb::Q
end

struct CoilGeometry{T<:Unitful.Length}
    id::T
    od::T
    len::T
end
CoilGeometry{T}(x::CoilGeometry) where {T<:Length} = CoilGeometry{T}(x.id,x.od,x.len)
Base.length(::CoilGeometry)=1
==(a::CoilGeometry, b::CoilGeometry) = a.id == b.id && a.od == b.od && a.len == b.len
≈(a::CoilGeometry, b::CoilGeometry) = a.id ≈ b.id && a.od ≈ b.od && a.len ≈ b.len
Base.:convert(::Type{T}, x::CoilGeometry) where {T<:CoilGeometry} = T(x)
Base.:promote_rule(::Type{CoilGeometry{T}},::Type{CoilGeometry{S}}) where {T,S} = CoilGeometry{promote_type(T,S)}

function CoilGeometry(id::Length,od::Length,len::Length, gap::Length=zero(id))
    CoilGeometry(promote(id+2gap, od-2gap, len-2gap)...)
end
function CoilGeometry(coil::CoilGeometry;id=coil.id,od=coil.od,len=coil.len)
    CoilGeometry(id,od,len)
end

struct WindingLayers
    layers::Int
    odd_layer::Int
    even_layer::Int
    odd_turns::Int
    even_turns::Int
    turns::Int
end

function WindingLayers(layers::T, odd_layer::T, even_layer::T) where {T<:Integer}
    div,rem = divrem(layers,2)
    odd_turns = odd_layer * (div+rem)
    even_turns = even_layer * div
    WindingLayers(layers, odd_layer, even_layer, odd_turns, even_turns, odd_turns+even_turns)
end
