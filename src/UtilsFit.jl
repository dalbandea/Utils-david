###
### "THE BEER-WARE LICENSE":
### Alberto Ramos wrote this file. As long as you retain this 
### notice you can do whatever you want with this stuff. If we meet some 
### day, and you think this stuff is worth it, you can buy me a beer in 
### return. <alberto.ramos@cern.ch>
###
### file:    UtilsFit.jl
### created: Tue Feb  2 16:00:40 2021
###

function fit_defs_yerr(f, x, dy)
    
    function lmfit(prm, y)
        
	nof = round(Int64, length(y))
	res = Vector{eltype(prm)}(undef,nof)
	for i in 1:nof
	    res[i] = (y[i] - f(x[i], prm)) / dy[i]
	end
        
	return res
    end
    chisq(prm,data) = sum(lmfit(prm, data) .^ 2)
    
    return lmfit, chisq
end

function fit_data_yerr(f, xv, yv, prms0, wpm = missing) # Argument prms0 and wpm added by David
    
    if ismissing(wpm)
        ADerrors.uwerr.(yv)
    else
        for uwvalue in yv
            ADerrors.uwerr(uwvalue, wpm)
        end
    end
    lm, csq = fit_defs_yerr(f, xv, ADerrors.err.(yv))
    
    prm0 = zeros(nparameters(f))
    
    # Control flow added by David
    if length(prms0) == nparameters(f)
	    prm0 = copy(prms0)
    else
	    println("Initial conditions do not match number of parameters")
    end

    fit  = LeastSquaresOptim.optimize(xx -> lm(xx, ADerrors.value.(yv)), prm0,
		    LeastSquaresOptim.LevenbergMarquardt(), autodiff = :forward)
    
    if ismissing(wpm)
        fitp, csqexp = ADerrors.fit_error(csq, fit.minimizer, yv)
    else
        fitp, csqexp = ADerrors.fit_error(csq, fit.minimizer, yv, wpm)
    end
    
    return fitp, csqexp, fit.ssr    
end

function fit_defs_xyerr(f, dx, dy)
    
    function lmfit(prm, xy)
        
	nof = round(Int64, length(xy)/2)
	res = Vector{eltype(prm)}(undef,2*nof)
        nfp = nparameters(f)
	for i in 1:nof
	    res[i] = (xy[i] - prm[nfp+i]) / dx[i]
	end
	for i in 1:nof
	    res[i+nof] = (xy[i+nof] - f(prm[i+nfp], view(prm,1:nfp))) / dy[i]
	end

	return res
    end
    chisq(prm,data) = sum(lmfit(prm, data) .^ 2)
    
    return lmfit, chisq
end

function fit_data_xyerr(f, xv, yv)
    
    ADerrors.uwerr.(yv)
    ADerrors.uwerr.(xv)
    lm, csq = fit_defs_xyerr(f, ADerrors.err.(xv), ADerrors.err.(yv))
    
    prm0 = zeros(nparameters(f)+length(xv))
    fit  = LeastSquaresOptim.optimize(xx -> lm(xx, ADerrors.value.([xv;yv])), prm0,
		    LeastSquaresOptim.LevenbergMarquardt(), autodiff = :forward)
    
    fitp, csqexp = ADerrors.fit_error(csq, fit.minimizer, [xv;yv])
    
    return fitp, csqexp, fit.ssr
end

function fit_data(f, xv, yv, prms0; wpm = missing) # Argument prms0 and wpm added by David
    if (isa(xv, Vector{ADerrors.uwreal}))
        return fit_data_xyerr(f, xv, yv)
    else
        return fit_data_yerr(f, xv, yv, prms0, wpm)
    end
end

function fit_defs_global_yerr(f, fc, x, dy, asq)
    
    function lmfit(prm, y)

        nfp = nparameters(f)
        nfc = nparameters(fc)
	nof = round(Int64, length(y))
	res = Vector{eltype(prm)}(undef,nof)
	for i in 1:nof
	    res[i] = (y[i] - (f(x[i], view(prm,1:nfp)) + fc(x[i], view(prm,nfp+1:nfp+nfc))*asq[i])) / dy[i]
	end
        
	return res
    end
    chisq(prm,data) = sum(lmfit(prm, data) .^ 2)
    
    return lmfit, chisq
end

function fit_data_global_yerr(f, fc, xv, yv, asq)
    
    ADerrors.uwerr.(yv)
    lm, csq = fit_defs_global_yerr(f, fc, xv, ADerrors.err.(yv), asq)
    
    prm0 = zeros(nparameters(f)+nparameters(fc))
    fit  = LeastSquaresOptim.optimize(xx -> lm(xx, ADerrors.value.(yv)), prm0,
		    LeastSquaresOptim.LevenbergMarquardt(), autodiff = :forward)
    
    fitp, csqexp = ADerrors.fit_error(csq, fit.minimizer, yv)
    
    return fitp, csqexp, fit.ssr    
end

function fit_defs_global_xyerr(f, fc, dx, dy, asq)
    
    function lmfit(prm, xy)
        
	nof = round(Int64, length(xy)/2)
	res = Vector{eltype(prm)}(undef,2*nof)
        nfp = nparameters(f)
        nfc = nparameters(fc)
	for i in 1:nof
	    res[i] = (xy[i] - prm[nfp+nfc+i]) / dx[i]
	end
	for i in 1:nof
	    res[i+nof] = (xy[i+nof] -
                          (f(prm[i+nfp+nfc], view(prm,1:nfp)) +
                           fc(prm[i+nfp+nfc], view(prm,nfp+1:nfp+nfc)) * asq[i])
                          ) / dy[i]
	end

	return res
    end
    chisq(prm,data) = sum(lmfit(prm, data) .^ 2)
    
    return lmfit, chisq
end

function fit_data_global_xyerr(f, fc, xv, yv, asq)
    
    ADerrors.uwerr.(yv)
    ADerrors.uwerr.(xv)
    lm, csq = fit_defs_global_xyerr(f, fc, ADerrors.err.(xv), ADerrors.err.(yv), asq)
    
    prm0 = zeros(nparameters(f)+nparameters(fc)+length(xv))
    for i in 1:length(xv)
        prm0[i+nparameters(f)+nparameters(fc)] = ADerrors.value(xv[i])
    end
    fit  = LeastSquaresOptim.optimize(xx -> lm(xx, ADerrors.value.([xv;yv])), prm0,
		    LeastSquaresOptim.LevenbergMarquardt(), autodiff = :forward)
    
    fitp, csqexp = ADerrors.fit_error(csq, fit.minimizer, [xv;yv])
    
    return fitp, csqexp, fit.ssr
end

function fit_data(f, fc, xv, yv, asq)
    if (isa(xv, Vector{ADerrors.uwreal}))
        return fit_data_global_xyerr(f, fc, xv, yv, asq)
    else
        return fit_data_global_yerr(f, fc, xv, yv, asq)
    end
end
