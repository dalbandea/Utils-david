###
### "THE BEER-WARE LICENSE":
### Alberto Ramos wrote this file. As long as you retain this 
### notice you can do whatever you want with this stuff. If we meet some 
### day, and you think this stuff is worth it, you can buy me a beer in 
### return. <alberto.ramos@cern.ch>
###
### file:    UtilsInt.jl
### created: Sun Feb 21 11:26:20 2021
###                               

function int_quad(ns)

    nd, w = FastGaussQuadrature.gausslegendre(ns)
    function do_int(a, b, f)
	r1 = (b-a)/2.0
	r2 = (b+a)/2.0
        
	sum = r1 * w[1] * f(r1*nd[1]+r2) 
	for i in 2:ns
	    sum += r1 * w[i] * f(r1*nd[i]+r2) 
	end
        
	return sum 
    end
    return do_int
end

