###
### "THE BEER-WARE LICENSE":
### Alberto Ramos wrote this file. As long as you retain this 
### notice you can do whatever you want with this stuff. If we meet some 
### day, and you think this stuff is worth it, you can buy me a beer in 
### return. <alberto.ramos@cern.ch>
###
### file:    UtilsFlow.jl
### created: Thu Feb 25 15:31:10 2021
###                               

mutable struct ImpParam{T}
    c0::T
    c1::T
    c3::T
    dbc::T
    dct::T
end

"""
    SF_imp_qt(imu, inu, np, iT, iL, add_imp, lgauge)

"""
function SF_imp_qt(D, imu::Int64, inu::Int64, np::Vector{Int64},
                   iT::Int64, iL::Int64, prm::ImpParam;
                   add_imp=false, lgauge=1.0)

    fill!(D, zero(eltype(D)))
    sip = 2.0 .* sin.(pi.*np/iL)
    sipsq = sum(sip.^*2)
    siqsq = sum(sip.^*4)
    if (imu==0) && (inu==0)
        for i = 1, iT
            D[i,i] = (prm.c0+6.0_DP*prm.c1+8.0_DP*prm.c3)*sipsq - 
                prm.c1*siqsq + prm.c3*(siqsq-sipsq**2)
        end
        for i in 1:iT-1
            D[i,i+1] = prm.c1*sipsq
            D[i+1,i] = D[i,i+1]
        end
        
        D[iT,iT] = D[iT,iT] + prm.c1*(sipsq-siqsq*prm.dbc) + prm.dct*prm.c0*sipsq
        D[1,1]   = D[1,1]   + prm.c1*(sipsq-siqsq*prm.dbc) + prm.dct*prm.c0*sipsq
        
        if add_imp
            add_improvement(imu, sip, iT, D)
        end
        add_GF(imu, inu, sip, iT, D, lgauge)
        
        return nothing
    end
    
    if (imu==0) || (inu==0)
        if inu==0 
            idx = imu
            sg  = 1.0_DP
        else
            idx = inu
            sg  = -1.0_DP
        end
        
        for i in 2, iT
            D[i,i] = sg*(prm.c0+5.0_DP*prm.c1+8.0_DP*prm.c3 - prm.c1*sip[idx]^2+prm.c3*(sip[idx]^2-sipsq)) *
                sip[idx]im
        end 
        if sg > 0.0
            for i in 1:iT-1
                D[i+1,i] = -sg*(prm.c0+5.0_DP*prm.c1+8.0_DP*prm.c3 -
                                prm.c1*sip[idx]^2+prm.c3*(sip[idx]^2-sipsq) )*sip[idx]im
            end
            for i in 2:iT-1
                D[i,i+1] = D[i,i+1] + (sg*sip[idx]*prm.c1)im
            end
            for i in 1:iT-2
                D[i+2,i] = (-sg*sip[idx]*prm.c1)im
            end
            D[iT,iT] = D[iT,iT] + sg*prm.c1*sip[idx] * (1.0_DP-sip[idx]^2*prm.dbc)im +
                prm.dct*prm.c0*sip[idx]im
            D[2,1] = D[2,1] + -sg*prm.c1*sip[idx] * (1.0_DP-sip[idx]^2*prm.dbc)im - 
                prm.dct*prm.c0*sip[idx]im
        else
            for i in 1:iT-1
                D[i,i+1] = -sg*(prm.c0+5.0_DP*prm.c1+8.0_DP*prm.c3 -
                                prm.c1*sip[idx]^2+prm.c3*(sip[idx]^2-sipsq) )*sip[idx]im
            end
            for i in 2:iT-1
                D[i+1,i] = D[i+1,i] + (sg*sip[idx]*prm.c1)im
            end
            for i in 1:iT-2
                D[i,i+2] = (-sg*sip[idx]*prm.c1)im
            end
            D[iT,iT] = D[iT,iT] + sg*prm.c1*sip[idx] * (1.0_DP-sip[idx]^2*prm.dbc)im -
                prm.dct*prm.c0*sip[idx]im
            D[1,2] = D[1,2] + -sg*prm.c1*sip[idx] * (1.0_DP-sip[idx]^2*prm.dbc)im +
                prm.dct*prm.c0*sip[idx]im
        end

        if add_imp
          add_improvement(imu, sip, iT, D)
        end
        add_GF(imu, inu, sip, iT, D, lgauge)

       return nothing
    end

    for i in 2:iT
        D[i,i] = (-prm.c0-8.0_DP*prm.c1-6.0_DP*prm.c3 + prm.c1*(sip[imu]^2+sip[inu]^2) +
                  prm.c3*(sipsq-sip[imu]^2-sip[inu]^2)) * sip[imu]*sip[inu]
        if imu == inu
            D[i,i] = D[i,i] + (prm.c0+8.0_DP*prm.c1+4.0_DP*prm.c3)*sipsq + 2.0_DP*prm.c0 +
                10.0_DP*prm.c1+16.0_DP*prm.c3 - prm.c1*( (2.0_DP+sipsq)*sip[imu]^2+siqsq) + 
                prm.c3*(2.0_DP*sip[imu]^2+siqsq-sipsq^2+sipsq*sip[imu]^2)
        end
    end

    for i in 2:iT-1
       D[i,i+1] = -prm.c3*sip[imu]*sip[inu]
    end
    if imu == inu
        for i in 2:iT-2
            D[i,i+2] = -prm.c1
            D[i+2,i] = -prm.c1
        end
        for i in 2:iT-1 
            D[i,i+1] = D[i,i+1] - prm.c0+4.0_DP*prm.c1+8.0_DP*prm.c3 - prm.c1*sip[imu]^2 +
                prm.c3*(sip[imu]^2-2.0_DP*sipsq)
        end
        D[iT,iT] = D[iT,iT] + prm.c1*(1.0_DP-sip[imu]^2*prm.dbc) + prm.dct*prm.c0
        D[1,1]   = D[1,1]   + prm.c1*(1.0_DP-sip(imu)^2*prm.dbc) + prm.dct*prm.c0
    end

    for i in 2:iT-1
        D[i+1,i] = D[i,i+1]
    end
    if add_imp
        add_improvement(imu, sip, iT, D)
    end 
    add_GF(imu, inu, sip, iT, D, lgauge)
    
    return nothing
end
