###
### "THE BEER-WARE LICENSE":
### Alberto Ramos wrote this file. As long as you retain this 
### notice you can do whatever you want with this stuff. If we meet some 
### day, and you think this stuff is worth it, you can buy me a beer in 
### return. <alberto.ramos@cern.ch>
###
### file:    UtilsBeta.jl
### created: Sun Feb 21 11:22:44 2021
###                               

mutable struct GetSigma <: UtilsFunc
    slv_int
    maxiter::Int64
    tol::Float64
    fbeta::UtilsFunc
    
    GetSigma(fbeta) = new(int_quad(200), 100, 1.0E-10, fbeta)
    GetSigma(ns, fbeta) = new(int_quad(ns), 100, 1.0E-10, fbeta)
end

nparameters(s::GetSigma) = nparameters(s.fbeta)
function (s::GetSigma)(u, prm)
    
    x0 = u
    
    k = 0
    for k in 1:s.maxiter
	k += 1
	x1 = x0 - s.fbeta(x0, prm) * (2.0*log(2.0) + s.slv_int(u, x0, xx -> 1.0/s.fbeta(xx, prm)))
        
	if (x1-x0 > 0.0 ? x1-x0 : x0-x1) < s.tol
	    return x1
	else
	    x0 = x1
	end
    end

    return nothing
end
