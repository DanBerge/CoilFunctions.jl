
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
