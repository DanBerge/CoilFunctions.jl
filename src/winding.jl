
oddturns(w,d) = Int(floor(round(w/d;digits=6)))
eventurns(w,d) = Int(floor(round((w-d/2)/d;digits=6)))

oddwidth(n::Integer,d::Unitful.Length) = n*d
evenwidth(n::Integer,d::Unitful.Length) = (n+0.5)*d
function layerheight(n::Integer,d::Unitful.Length)
    @assert n>0
    y = d
    if n>1
        y += (n-1)*(d/2)*√3
    end
    y
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
    odd_layer= oddturns(w,d)
    even_layer= eventurns(w,d)

    odd_turns = length(1:2:layers)*odd_layer
    even_turns = length(2:2:layers)*even_layer

    WindingLayers(layers, odd_layer, even_layer, odd_turns, even_turns, odd_turns+even_turns)
end

"""
ideal_fill(coil::CoilGeometry, awg::AWG, [fill::Real=1.0])

    coil => CoilGeometry
    awg => AWG of the wire
    fillfactor => fillfactor of the wire
"""
function ideal_fill(coil::CoilGeometry, awg::AWG, fillfactor::T=1.0) where {T<:Real}
    zero(T) < fillfactor <= one(T) || throw(DomainError(fillfactor))
    w=coil.len
    h=(coil.od-coil.id)/2
    ideal_fill(w,h,awg.d_insul/sqrt(fillfactor))
end

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
ideal_layers(h,d) = h≳d ? 1+Int(round((h-d)/(d*√3/2))) : 0
ideal_layers(coil::CoilGeometry,d::Unitful.Length) = ideal_layers((coil.od-coil.id)/2,d)

"""
Construct the minimum size coil to fit the specified winding

    enclosewinding(layers, oddlayer, evenlayer, id::Unitful.Length, d::Unitful.Length)
    enclosewinding(layers, oddlayer, evenlayer, id::Unitful.Length, gauge::AWG)
"""
function enclosewinding(layers, oddlayer, evenlayer, id::Unitful.Length, d::Unitful.Length)
    @assert layers>0
    @assert evenlayer == oddlayer || evenlayer == oddlayer-1
    # od = id + 2d + (layers-1)*d*√3
    # len = oddlayer > evenlayer ? d*oddlayer : d/2 + d*evenlayer
    w=oddlayer>evenlayer ? oddwidth(oddlayer,d) : evenwidth(evenlayer,d)
    od=2layerheight(layers,d) + id
    f = unit(id)
    CoilGeometry(id,f(od),f(w))
end
enclosewinding(layers, oddlayer, evenlayer, id::Unitful.Length, gauge::AWG) =
    enclosewinding(layers, oddlayer, evenlayer, id, gauge.d_insul)

"""
Construct the minimum size `coil` to fit the `winding` within the limits of the existing coil.

enclosewinding(coil::CoilGeometry, winding::WindingLayers, d::Unitful.Length)
enclosewinding(coil::CoilGeometry, winding::WindingLayers, gauge::AWG)
"""
function enclosewinding(coil::CoilGeometry, winding::WindingLayers, d::Unitful.Length)
#TODO Set check to make sure the calculated winding is within the limits of the existing coil
    new_coil = enclosewinding(winding.layers, winding.odd_layers, winding.even_layers, coil.od, d)
    @assert new_coil.id >= coil.id && new_coil.od <= coil.od && new_coil.len <= coil.len
    @assert isvalidcoil(new_coil)
    new_coil
end
enclosewinding(coil::CoilGeometry, winding::WindingLayers, gauge::AWG) = enclosewinding(coil,winding,gauge.d_insul)
