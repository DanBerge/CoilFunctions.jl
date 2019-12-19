function plot_circle!(plt, x, y, r; points=100, color=:black)
    rng = range(0; stop=2pi, length=points)
    plot!(plt, θ -> x + r*cos(θ), θ -> y + r*sin(θ), rng; color=color)
end

function plot_winding(winding, awg, maxturns, fill=1.0)

    coil=winding.coil

    width = (coil.od-coil.id)/2 |> ustrip
    h = ustrip(coil.len)
    plt=plot([0, h,h, 0, 0], [0, 0, width, width, 0]; aspect_ratio=:equal,legend=:none)

    d = awg.d_insul |>  unit(winding.coil.id) |> ustrip
    r = d/2
    d /= sqrt(fill)
    turns = 0
    foreach(enumerate(range(d/2;step=sqrt(3)/2*d,length=winding.winding.layers))) do (layer,y)
        x_offset = isodd(layer) ? d/2 : d
        layerturns = isodd(layer) ? winding.winding.odd_layer : winding.winding.even_layer
        foreach(range(x_offset;step=d,length=layerturns)) do x
            if turns < maxturns
                plot_circle!(plt,x,y,r)
                turns += 1
            else

            end
        end
    end
    plt
end
