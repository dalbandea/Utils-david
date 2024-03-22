###
### "THE BEER-WARE LICENSE":
### Alberto Ramos wrote this file. As long as you retain this 
### notice you can do whatever you want with this stuff. If we meet some 
### day, and you think this stuff is worth it, you can buy me a beer in 
### return. <alberto.ramos@cern.ch>
###
### file:    UtilsMSbar.jl
### created: Tue Feb  9 12:03:50 2021
###                               

beta_b0(Nf::Int64) = (11.0-2.0*Nf/3.0)/(4.0*pi)^2
beta_b1(Nf::Int64) = (102.0 - (8.0/3.0 + 10.0)*Nf) / (4.0*pi)^4
beta_b2MS(Nf::Int64) = 1.0/4.0^3 *
    ( 2857.0/2.0 - 5033.0/18.0 * Nf + 325.0/54.0 * Nf^2 ) /
    (4.0*pi^2)^3
beta_b3MS(Nf::Int64) = 1.0/4.0^4 * ( 
    149753.0/6.0 + 3564.0*1.2020569031595942853997 - 
    ( 1078361.0/162.0 + 6508.0/27.0*1.2020569031595942853997 ) * Nf + 
    ( 50065.0/162.0 + 6472.0/81.0*1.2020569031595942853997   ) * Nf^2 + 
    ( 1093.0/729.0 ) * Nf^3 ) /
    (4.0*pi^2)^4
beta_b4MS(Nf::Int64) = (524.56 - 181.8*Nf + 17.16*Nf^2 - 0.22586*Nf^3 - 0.0017993*Nf^4) /
    (4.0*pi^2)^5


function betaMS(nf::Int64; nl::Int64 = 5)
    
    res = Vector{Float64}(undef, nl)
    res[1] = beta_b0(nf)
    if nl==1
        return res
    end
    
    res[2] = beta_b1(nf)
    if nl==2
        return res
    end

    res[3] = beta_b2MS(nf)
    if nl==3
        return res
    end

    res[4] = beta_b3MS(nf)
    if nl==4
        return res
    end

    res[5] = beta_b4MS(nf)
    if nl==5
        return res
    end

    return res
end
