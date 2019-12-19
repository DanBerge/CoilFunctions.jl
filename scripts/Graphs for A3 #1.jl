fdir = "J:\\ENG\\PAC_Group\\CurrentProjects\\Electric ServoMotor Brake\\Plan & Proof\\Plan\\2019-9-13 Coil Winding"
_savefig = (plt,fname)-> begin
    fdir = "J:\\ENG\\PAC_Group\\CurrentProjects\\Electric ServoMotor Brake\\Plan & Proof\\Plan\\2019-9-13 Coil Winding"
    savefig(plt,joinpath(fdir,fname*".pdf"))
    savefig(plt,joinpath(fdir,fname*".png"))
end

plt=let
    step=1u"mm"
    id=70u"mm":step:140u"mm"
    od=80u"mm":step:150u"mm"
    len=26u"mm"
    power=50.0u"W"
    global f = x -> begin
        coil=CoilGeometry(x.id,x.od,x.len)
        if CoilFunctions.isvalidcoil(coil)
            CoilFunctions.minimumturns(coil,power,1.0)
        else
            NaN
        end
    end
    @time y=tmap(f,((id=id,od=od,len=len) for id in id, od in od))
    contour(ustrip.(id),ustrip.(od),y';fill=true,
        title="SBE7 24V Minimum Equivalent Turns",xlabel="ID (mm)",ylabel="OD (mm)",dpi=300)
    plot!([0, 100,100],[140,140,0];color=:red,linestyle=:dash,
        ylims=ustrip.(extrema(od)),xlims=ustrip.(extrema(id)),
        legend=:bottomright,label="SBE7 Reference")
end
_savefig(plt,"SBE7 24V Minimum Equivalent Turns")

plt=let
    step=1u"mm"
    id=70u"mm":step:140u"mm"
    od=80u"mm":step:150u"mm"
    len=26u"mm"
    power=50.0u"W"
    global f = x -> begin
        coil=CoilGeometry(x.id,x.od,x.len)
        if CoilFunctions.isvalidcoil(coil)
            CoilFunctions.minimumturns(coil,power,1.0;op=extrema)
        else
            (NaN, NaN)
        end
    end
    @time y=tmap(f,((id=id,od=od,len=len) for id in id, od in od))
    contour(ustrip.(id),ustrip.(od),map(x->x[2]-x[1],permutedims(y,(2,1)));fill=true,
        title="SBE7 24V Equivalent Range of Turn Count",xlabel="ID (mm)",ylabel="OD (mm)",dpi=300)
    plot!([0, 100,100],[140,140,0];color=:red,linestyle=:dash,
        ylims=ustrip.(extrema(od)),xlims=ustrip.(extrema(id)),
        legend=:bottomright,label="SBE7 Reference")
end
_savefig(plt,"SBE7 24V Equivalent Range of Turn Count")

plt=let
    step=0.1u"mm"
    id=70u"mm":step:140u"mm"
    od=80u"mm":step:150u"mm"
    len=26u"mm"
    power=50.0u"W"
    v = 24.0u"V"
    Ω = v^2/power |> u"Ω"
    global f = x -> begin
        coil=CoilGeometry(x.id,x.od,x.len)
        if CoilFunctions.isvalidcoil(coil)
            optimalcoil(coil,Ω,1.0).turns
        else
            NaN
        end
    end
    @time y=tmap(f,((id=id,od=od,len=len) for id in id, od in od))
    contour(ustrip.(id),ustrip.(od),y';fill=true,
        title="SBE7 24V Coil Minimum Turns",xlabel="ID (mm)",ylabel="OD (mm)",dpi=300)
    plot!([0, 100,100],[140,140,0];color=:red,linestyle=:dash,
        ylims=ustrip.(extrema(od)),xlims=ustrip.(extrema(id)),
        legend=:bottomright,label="SBE7 Reference")
end
_savefig(plt,"SBE7 24V Coil Minimum Turns")

plt=let
    step=0.1u"mm"
    id=80u"mm":step:120u"mm"
    od=130u"mm":step:150u"mm"
    len=26u"mm"
    power=50.0u"W"
    v = 24.0u"V"
    Ω = v^2/power |> u"Ω"
    global f = x -> begin
        coil=CoilGeometry(x.id,x.od,x.len)
        if CoilFunctions.isvalidcoil(coil)
            optimalcoil(coil,Ω,1.0).gauge
        else
            NaN
        end
    end
    @time y=tmap(f,((id=id,od=od,len=len) for id in id, od in od))
    contour(ustrip.(id),ustrip.(od),y';fill=true,
        title="SBE7 24V Coil Wire Gauge",xlabel="ID (mm)",ylabel="OD (mm)",dpi=300,
        levels=length(unique(y)))
end
_savefig(plt,"SBE7 24V Coil Wire Gauge")
