###
### "THE BEER-WARE LICENSE":
### Alberto Ramos wrote this file. As long as you retain this 
### notice you can do whatever you want with this stuff. If we meet some 
### day, and you think this stuff is worth it, you can buy me a beer in 
### return. <alberto.ramos@cern.ch>
###
### file:    UtilsFunctions.jl
### created: Tue Feb  2 12:04:49 2021
###                               


abstract type UtilsFunc end
abstract type MixedUtilsFunc <: UtilsFunc end

mutable struct FuncSeries{T1,T2,T3} <: UtilsFunc
    ip::Vector{T1}
    fc::Vector{T2}
    fp::Vector{T3}
end

FuncSeries(ip;coeffs=Vector{Float64}(),powers=Vector{Int64}()) = FuncSeries(ip,coeffs,powers)

mutable struct FuncPade <: UtilsFunc
    nn::Int64
    nd::Int64
end

nparameters(s::FuncPade)   = s.nn+s.nd+1
nparameters(s::FuncSeries) = length(s.ip)

function (s::FuncSeries)(x, p)

    if (length(p) == 0) && (length(s.fc) == 0)
        return 0.0
    end
    
    if (length(p) > 0)
        res = p[1]*x^(s.ip[1])
        for k in 2:length(p)
            res = res + p[k]*x^(s.ip[k])
        end
    else
        res = 0.0
    end

    if (length(s.fc) > 0)
        res = res + s.fc[1]*x^(s.fp[1])
        for k in 2:length(s.fc)
            res = res + s.fc[k]*x^(s.fp[k])
        end
    end
    
    return res
end

function (s::FuncPade)(x, prm)

    if (s.nn == 0) && (s.nd == 0)
        return 0.0
    end

    if (s.nn >= 0)
        vn = prm[1]
        for k in 2:s.nn+1
            vn = vn + prm[k]*x^(k-1)
        end
    else
        vn = 0.0
    end

    if (s.nd > 0)
        vd = 1.0 + prm[s.nn+2]*x
        for k in 2:s.nd
            vd = vd + prm[nn+k+1]*x^k
        end
    end

    return vn/vd
end


#
# Simple correlator functions
#
# Correlator(T/a, sign)
#
# Is a function of the form
#  - s=0:  A * e^{-mx_0}
#  - s=1:  A * (e^{-mx_0} + e^{-m(T/a-x_0)})
#  - s=-1: A * (e^{-mx_0} - e^{-m(T/a-x_0)})
#
mutable struct Correlator <: UtilsFunc
    T::Int64
    s::Int64
end

nparameters(s::Correlator) = 2
(s::Correlator)(x, p) = p[1] * (exp(-p[2]*x) + s.s*exp(-p[2]*(s.T-x)))


mutable struct ConstantFunc <: UtilsFunc
end

nparameters(s::ConstantFunc) = 1

function (s::ConstantFunc)(x,p)
	return p[1]
end


mutable struct Einstein <: UtilsFunc
end

nparameters(s::Einstein) = 1

function (s::Einstein)(x,p)
	return p[1]^2 + x[1]^2
end

@Base.kwdef mutable struct SymCorrelator <: UtilsFunc
    T::Int64
    s::Int64
    c::Int64 = 0
end

# Symmetric Correlator
nparameters(s::SymCorrelator) = abs(s.s)*2 + s.c
# (s::SymCorrelator)(x, p) = p[1] * cosh( p[2]*(x-s.T/2) )
# (s::SymCorrelator)(x, p) = p[1] * sinh( p[2]*(x-s.T/2) )
function (s::SymCorrelator)(x,p)
	if(s.s == 1)
		if(s.c == 0)
			return p[1] * cosh( p[2]*(x-s.T/2) )
		elseif(s.c == 1)
			return p[1] * cosh( p[2]*(x-s.T/2) ) + p[3]
		elseif(s.c == 2)
			return p[1] * cosh( p[2]*(x-s.T/2) ) + p[3] * exp(-p[4]*x)
		end
	elseif(s.s == -1)
		return p[1] * sinh( p[2]*(x-s.T/2) )
	elseif(s.s == -2)
		return p[1] * sinh( p[2]*(x-s.T/2) ) + p[3] * exp(-p[4]*x)
	end

end


# Mixed Correlator, with input from two different correlators
@Base.kwdef mutable struct MixedCorrelator <: MixedUtilsFunc
    T::Int64
    s::Int64
    c::Int64 = 0
end


function (s::MixedCorrelator)(x,p)
	if(s.s == 1)
		return x[2] * p[1] * cosh( p[2]*(x[1] - s.T/2) ) + p[3] * (1 .- x[2] ) * sinh( p[2]*(x[1] - s.T/2) )
	elseif(s.s == -1)
		return x[2] * p[1] * sinh( p[2]*(x[1] - s.T/2) ) + p[3] * (1 .- x[2] ) * sinh( p[2]*(x[1] - s.T/2) )
	elseif(s.s == 2)
		return x[2] * (p[1] * sinh( p[2]*(x[1] - s.T/2) ) + p[4]*exp(-p[5]*x[1])) + (1 .- x[2] ) * (p[3]  * sinh( p[2]*(x[1] - s.T/2) ) + p[6]*exp(-p[5]*x[1]) )
	end
end

nparameters(s::MixedCorrelator) = 3*abs(s.s)


# Projected Symmetric Correlator
@Base.kwdef mutable struct ProjSymCorrelator <: UtilsFunc
    T::Int64
    s::Int64
end

function (s::ProjSymCorrelator)(x,p)
	if(s.s == 1)
		return p[1] * cosh( x[2]*(x[1]-s.T/2) )
	elseif(s.s == -2)
		return p[1] * sinh( x[2]*(x[1]-s.T/2) ) + p[2] * exp(-x[3]*x[1])
	end
end

nparameters(s::ProjSymCorrelator) = abs(s.s)


# Projected Mixed Correlator, with input from two different correlators
@Base.kwdef mutable struct ProjMixedCorrelator <: MixedUtilsFunc
    T::Int64
    s::Int64
    c::Int64 = 0
end


# x = [[x,M1,M2],α]
function (s::ProjMixedCorrelator)(x,p)
	if(s.s == 2)
		return x[2] * (p[1] * sinh( x[1][2]*(x[1][1] - s.T/2) ) + p[3]*exp(-x[1][3]*x[1][1])) + (1 .- x[2]) * (p[2]  * sinh( x[1][2]*(x[1][1] - s.T/2) ) + p[4]*exp(-x[1][3]*x[1][1]) )
	end
end

nparameters(s::ProjMixedCorrelator) = 2*abs(s.s)


mutable struct AutCorrScaling <: UtilsFunc end
nparameters(s::AutCorrScaling) = 3
(s::AutCorrScaling)(x, p) = p[1] * x^p[2] + p[3]


mutable struct AutCorrScalingwoOffset <: UtilsFunc end
nparameters(s::AutCorrScalingwoOffset) = 2
(s::AutCorrScalingwoOffset)(x, p) = p[1] * x^p[2]
