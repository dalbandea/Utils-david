###
### "THE BEER-WARE LICENSE":
### Alberto Ramos wrote this file. As long as you retain this 
### notice you can do whatever you want with this stuff. If we meet some 
### day, and you think this stuff is worth it, you can buy me a beer in 
### return. <alberto.ramos@cern.ch>
###
### file:    UtilsPlot.jl
### created: Wed Feb  3 11:45:01 2021
###                               


plot_data_yerr(xv, yv, dtlb; kwargs...) = Plots.plot!(xv, ADerrors.value.(yv);
                                                      label=dtlb,
                                                      seriestype=:scatter,
                                                      yerror=ADerrors.err.(yv),
                                                      kwargs...)
plot_data_xyerr(xv, yv, dtlb) = Plots.plot!(ADerrors.value.(xv),
                                            ADerrors.value.(yv),
                                            yerror=ADerrors.err.(yv),
                                            xerror=ADerrors.err.(xv),
                                            label=dtlb, seriestype=:scatter)
function plot_func(f, prm, xlim; ns=20, lbl="Fit")

    xv = range(xlim[1], stop=xlim[2], length=ns)
    fv = Vector{ADerrors.uwreal}(undef, ns)
    for i in 1:ns
        fv[i] = f(xv[i], prm)
    end
    ADerrors.uwerr.(fv)

    return Plots.plot(xv, ADerrors.value.(fv), ribbons=ADerrors.err.(fv), label=lbl, reuse=false) # David changed this from Plots.plot!(...) to Plots.plot(..., reuse=false)
end

function plot_func!(pl, f, prm, xlim; ns=20, lbl="Fit", kwargs...)

    xv = range(xlim[1], stop=xlim[2], length=ns)
    fv = Vector{ADerrors.uwreal}(undef, ns)
    for i in 1:ns
        fv[i] = f(xv[i], prm)
    end
    ADerrors.uwerr.(fv)

    return Plots.plot!(pl, xv, ADerrors.value.(fv); ribbons=ADerrors.err.(fv), label=lbl, kwargs...) # David changed this from Plots.plot!(...) to Plots.plot(..., reuse=false)
end

function plot_data(xv, yv, dtlb)

    ADerrors.uwerr.(yv)
    if eltype(xv) == ADerrors.uwreal
        ADerrors.uwerr.(xv)
        return plot_data_xyerr(xv,yv, dtlb)
    else
        return plot_data_yerr(xv,yv, dtlb)
    end
end
    
function plot_fit(xv, yv, f, prm; dtlb="Y", flbl="Fit", ns=20, pborder=20)

    ADerrors.uwerr.(yv)
    if eltype(xv) == ADerrors.uwreal
        ADerrors.uwerr.(xv)
        plot_data_xyerr(xv,yv, dtlb)
        xmin = minimum(ADerrors.value.(xv))
        xmax = maximum(ADerrors.value.(xv))
    else
        plot_data_yerr(xv,yv, dtlb)
        xmin = minimum(xv)
        xmax = maximum(xv)
    end
    xlims = (xmin, xmax)
    
    return plot_func(f, prm, xlims, ns=ns, lbl=flbl)
end


function plot_fit!(pl, xv, yv, f, prm; dtlb="Y", flbl="Fit", ns=20, pborder=20, kwargs...)

    ADerrors.uwerr.(yv)
    if eltype(xv) == ADerrors.uwreal
        ADerrors.uwerr.(xv)
        plot_data_xyerr(xv,yv, dtlb)
        xmin = minimum(ADerrors.value.(xv))
        xmax = maximum(ADerrors.value.(xv))
    else
        plot_data_yerr(xv,yv, dtlb; kwargs...)
        xmin = minimum(xv)
        xmax = maximum(xv)
    end
    xlims = (xmin, xmax)
    
    return plot_func!(pl, f, prm, xlims; ns=ns, lbl=flbl, kwargs...)
end
