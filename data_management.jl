######################################################################################
# This file was created by Alessandro Conigli 
# Here we define routines for reading and computing the correlators G33 G88 and G08
# used in the HVP analyses. We also define routines for the corresponding renormalization 
# constants Z3 Z8 and Z08
########################################################################################

function get_Z3(ens::EnsInfo; impr_set::String="1")
    
    ml = 0.5 * (1 / ens.kappa_l - 1 / ens.kappa_crit )
    ms = 0.5 * (1 / ens.kappa_s - 1 / ens.kappa_crit )
    trmq = 2*ml + ms
    if impr_set == "1"
        Z3 = ZV(ens.beta) * (1. + bv_bar(ens.beta)*trmq + bv(ens.beta)*ml)
    elseif  impr_set == "2"
        Z3 = ZV_set2(ens.beta) * (1. + bv_bar_set2(ens.beta)*trmq + bv_set2(ens.beta)*ml)
    end
    return Z3
end

function get_Z8(ens::EnsInfo; impr_set::String="1")

    ml = 0.5 * (1 / ens.kappa_l - 1 / ens.kappa_crit )
    ms = 0.5 * (1 / ens.kappa_s - 1 / ens.kappa_crit )
    trmq = 2*ml + ms
    if impr_set == "1"
        Z8 = ZV(ens.beta) * (1. + bv_bar(ens.beta)*trmq + bv(ens.beta)*(ml + 2*ms)/3)
    elseif impr_set == "2"
        Z8 = ZV_set2(ens.beta) * (1. + bv_bar_set2(ens.beta)*trmq + bv_set2(ens.beta)*(ml + 2*ms)/3)
    end
    return Z8
end

function get_Z08(ens::EnsInfo; impr_set::String="1")

    ml = 0.5 * (1 / ens.kappa_l - 1 / ens.kappa_crit )
    ms = 0.5 * (1 / ens.kappa_s - 1 / ens.kappa_crit )
    if impr_set == "1"
        Z08 = 1/3 * ZV(ens.beta) * bv(ens.beta) * (2/sqrt(3)) * (ml - ms)
    elseif impr_set == "2"
        Z08 = 1/3 * ZV_set2(ens.beta) * bv_set2(ens.beta) * (2/sqrt(3)) * (ml - ms)
    end
    return Z08
end

function renormalize!(cc::Corr, Z::uwreal)
    cc.obs[:] .*= Z
    return nothing
end
function renormalize!(g3l::Vector{uwreal}, Z::uwreal)
    g3l[:] .*= Z
    return nothing
end

@doc raw"""
    corr33(path_data::String, ens::EnsInfo; path_rw::Union{Nothing,String}=nothing, L::Int64=1, frw_bcwd::Bool=true, impr::Bool=true, cons::Bool=true, std::Bool=false)

This function return the G33 correlator given path_data and the EnsInfo data type of the corresponding ensemble.  

Optional flags:    

    - path_rw  : if provided, correlators are reweighted.  
    - L        : correlators are normalised with the volume L^3. L=1 by default.    
    - frw_bcwd : if true, the forward backward symmetrization is performed. 
    - impr     : if true, the correlator is improved.  
    - impr_set : either "1" or "2", to select set1 or set2 of improvement coefficients accordingly.  
    - cons     : if true, the function returns both the local-local and the local-conserved G88 correlators, else only the local-local is returned.  
    - std      : if true, standard symmetric derivatives are used for the vector-tensor correlator. If false, improved derivatives are used. The latter are thought to improve the short distance cutoff effects.

Examples:
```@example
g33_ll, g33_lc = corr33(path, ensinfo, path_rw=path_rw, frw_bcwd=true, impr=true, cons=true, std=false) 
```
"""
function corr33(path_data::String, ens::EnsInfo; sector::String="light", path_rw::Union{Nothing,String}=nothing, L::Int64=1, frw_bcwd::Bool=true, impr::Bool=true, impr_set::String="1", cons::Bool=true, std::Bool=false)
    
    Gamma_l = ["V1V1", "V2V2", "V3V3", "V1T10", "V2T20", "V3T30"]
    Gamma_c = ["V1V1c", "V2V2c", "V3V3c", "V1cT10", "V2cT20", "V3cT30"]

    v1v1 = get_corr(path_data, ens, sector, Gamma_l[1], path_rw=path_rw, frw_bcwd=false, L=L)
    v2v2 = get_corr(path_data, ens, sector, Gamma_l[2], path_rw=path_rw, frw_bcwd=false, L=L)
    v3v3 = get_corr(path_data, ens, sector, Gamma_l[3], path_rw=path_rw, frw_bcwd=false, L=L)
    
    v1t10 = get_corr(path_data, ens, sector, Gamma_l[4], path_rw=path_rw, frw_bcwd=false, L=L)
    v2t20 = get_corr(path_data, ens, sector, Gamma_l[5], path_rw=path_rw, frw_bcwd=false, L=L)
    v3t30 = get_corr(path_data, ens, sector, Gamma_l[6], path_rw=path_rw, frw_bcwd=false, L=L)

    if impr
        beta = ens.beta
        if impr_set == "1"
            cv_l = cv_loc(beta)
        elseif impr_set =="2"
            cv_l = cv_loc_set2(beta)
        end
        improve_corr_vkvk!(v1v1, v1t10, 2*cv_l, std=std)
        improve_corr_vkvk!(v2v2, v2t20, 2*cv_l, std=std)
        improve_corr_vkvk!(v3v3, v3t30, 2*cv_l, std=std)
    end
    g33 = Corr(0.5 .* (v1v1.obs .+ v2v2.obs .+ v3v3.obs)./3, v1v1.id, "G33ll")
    if frw_bcwd
        frwd_bckwrd_symm!(g33)
    end

    if cons
        v1v1_c = get_corr(path_data, ens, sector, Gamma_c[1], path_rw=path_rw, frw_bcwd=false, L=L)
        v2v2_c = get_corr(path_data, ens, sector, Gamma_c[2], path_rw=path_rw, frw_bcwd=false, L=L)
        v3v3_c = get_corr(path_data, ens, sector, Gamma_c[3], path_rw=path_rw, frw_bcwd=false, L=L)
    
        v1t10_c = get_corr(path_data, ens, sector, Gamma_c[4], path_rw=path_rw, frw_bcwd=false, L=L)
        v2t20_c = get_corr(path_data, ens, sector, Gamma_c[5], path_rw=path_rw, frw_bcwd=false, L=L)    
        v3t30_c = get_corr(path_data, ens, sector, Gamma_c[6], path_rw=path_rw, frw_bcwd=false, L=L)

        if impr
            beta = ens.beta
            if impr_set == "1"
                cv_l = cv_loc(beta) 
                cv_c = cv_cons(beta)
            elseif impr_set == "2"
                cv_l = cv_loc_set2(beta)
                cv_c = cv_cons(beta)
            end
            improve_corr_vkvk_cons!(v1v1_c, v1t10, v1t10_c, cv_l, cv_c, std=std)
            improve_corr_vkvk_cons!(v2v2_c, v2t20, v2t20_c, cv_l, cv_c, std=std)
            improve_corr_vkvk_cons!(v3v3_c, v3t30, v3t30_c, cv_l, cv_c, std=std)
        end
        g33_c = Corr(0.5 .* (v1v1_c.obs .+ v2v2_c.obs .+ v3v3_c.obs)./3, v1v1_c.id, "G33lc")
        if frw_bcwd
            frwd_bckwrd_symm!(g33_c)
        end
    end

    !cons ? (return g33) : (return g33, g33_c)
end

@doc raw"""
    corr88_conn(path_data::String, ens::EnsInfo, g33_ll::Corr; g33_lc::Union{Nothing, Corr}=nothing, path_rw::Union{Nothing,String}=nothing, L::Int64=1, frw_bcwd::Bool=true, impr::Bool=true, cons::Bool=true, std::Bool=true)

This function return the G88 connected correlator given path_data and the EnsInfo data type for the ensemble of interest.  

Optional flags:  

    - g33_lc   : either nothing or a local-conserved G33 correlator. This is used to build the local-conserevd G88. It is required only if cons is set to true.
    - path_rw  : if provided, correlators are reweighted.  
    - L        : correlators are normalised with the volume L^3. L=1 by default.    
    - frw_bcwd : if true, the forward backward symmetrization is performed. 
    - impr     : if true, the correlator is improved.  
    - impr_set : either "1" or "2", to select set1 or set2 of improvement coefficients accordingly.  
    - cons     : if true, the function returns both the local-local and the local-conserved G88 correlators, else only the local-local is returned.  
    - std      : if true, standard symmetric derivatives are used for the vector-tensor correlator. If false, improved derivatives are used. The latter are thought to improve the short distance cutoff effects.

Examples:
```@example
corr88_conn(path, ensinfo, path_rw=path_rw, frw_bcwd=true, impr=true, std=true) 
```
"""
function corr88_conn(path_data::String, ens::EnsInfo, g33_ll::Corr; g33_lc::Union{Nothing, Corr}=nothing, path_rw::Union{Nothing,String}=nothing, L::Int64=1, frw_bcwd::Bool=true, impr::Bool=true, impr_set::String="1", cons::Bool=true, std::Bool=true)

    Gamma_l = ["V1V1", "V2V2", "V3V3", "V1T10", "V2T20", "V3T30"]
    Gamma_c = ["V1V1c", "V2V2c", "V3V3c", "V1cT10", "V2cT20", "V3cT30"]


    v1v1 = get_corr(path_data, ens, "strange", Gamma_l[1], path_rw=path_rw, frw_bcwd=false, L=L)
    v2v2 = get_corr(path_data, ens, "strange", Gamma_l[2], path_rw=path_rw, frw_bcwd=false, L=L)
    v3v3 = get_corr(path_data, ens, "strange", Gamma_l[3], path_rw=path_rw, frw_bcwd=false, L=L)
    
    v1t10 = get_corr(path_data, ens, "strange", Gamma_l[4], path_rw=path_rw, frw_bcwd=false, L=L)
    v2t20 = get_corr(path_data, ens, "strange", Gamma_l[5], path_rw=path_rw, frw_bcwd=false, L=L)
    v3t30 = get_corr(path_data, ens, "strange", Gamma_l[6], path_rw=path_rw, frw_bcwd=false, L=L)

    if impr
        beta = ens.beta
        if impr_set == "1"
            cv_l = cv_loc(beta)
        elseif impr_set =="2"
            cv_l = cv_loc_set2(beta)
        end 
        improve_corr_vkvk!(v1v1, v1t10, 2*cv_l, std=std)
        improve_corr_vkvk!(v2v2, v2t20, 2*cv_l, std=std)
        improve_corr_vkvk!(v3v3, v3t30, 2*cv_l, std=std)
    end

    g88_ll = Corr(1/6 .* (2 * g33_ll.obs + 2 * (v1v1.obs + v2v2.obs + v3v3.obs) /3 ), v1v1.id, "G88ll") 
    if frw_bcwd
        frwd_bckwrd_symm!(g88_ll)
    end

    if cons
        if isnothing(g33_lc)
            error("G33_lc is required to compute G88_lc")
        end
        v1v1_c = get_corr(path_data, ens, "strange", Gamma_c[1], path_rw=path_rw, frw_bcwd=false, L=L)
        v2v2_c = get_corr(path_data, ens, "strange", Gamma_c[2], path_rw=path_rw, frw_bcwd=false, L=L)
        v3v3_c = get_corr(path_data, ens, "strange", Gamma_c[3], path_rw=path_rw, frw_bcwd=false, L=L)
    
        v1t10_c = get_corr(path_data, ens, "strange", Gamma_c[4], path_rw=path_rw, frw_bcwd=false, L=L)
        v2t20_c = get_corr(path_data, ens, "strange", Gamma_c[5], path_rw=path_rw, frw_bcwd=false, L=L)
        v3t30_c = get_corr(path_data, ens, "strange", Gamma_c[6], path_rw=path_rw, frw_bcwd=false, L=L)
        
        if impr
            beta = ens.beta
            if impr_set == "1"
                cv_l = cv_loc(beta) 
                cv_c = cv_cons(beta)
            elseif impr_set == "2"
                cv_l = cv_loc_set2(beta)
                cv_c = cv_cons_set2(beta)
            end
            improve_corr_vkvk_cons!(v1v1_c, v1t10, v1t10_c, cv_l, cv_c, std=std)
            improve_corr_vkvk_cons!(v2v2_c, v2t20, v2t20_c, cv_l, cv_c, std=std)
            improve_corr_vkvk_cons!(v3v3_c, v3t30, v3t30_c, cv_l, cv_c, std=std)
        end
        g88_lc = Corr(1/6 .* (2 * g33_lc.obs + 2 * (v1v1_c.obs + v2v2_c.obs + v3v3_c.obs) /3 ), v1v1.id, "G88lc") 
        if frw_bcwd
            frwd_bckwrd_symm!(g88_lc)
        end
    end 

    !cons ? (return g88_ll) : (return g88_ll, g88_lc) 
end

function corr88_disc(path_data::String, ens::EnsInfo; path_rw::Union{Nothing,String}=nothing, L::Int64=1, frw_bcwd::Bool=true, impr::Bool=true, impr_set::String="1", cons::Bool=true, std::Bool=true)

    Gamma_l = ["V1V1", "V2V2", "V3V3", "V1T10", "V2T20", "V3T30"]
    Gamma_c = ["V1V1c", "V2V2c", "V3V3c", "V1cT10", "V2cT20", "V3cT30"]

    v1v1 = get_corr(path_data, ens, "88", Gamma_l[1], path_rw=path_rw, frw_bcwd=false, L=L)
    v2v2 = get_corr(path_data, ens, "88", Gamma_l[2], path_rw=path_rw, frw_bcwd=false, L=L)
    v3v3 = get_corr(path_data, ens, "88", Gamma_l[3], path_rw=path_rw, frw_bcwd=false, L=L)
    
    v1t10 = get_corr(path_data, ens, "88", Gamma_l[4], path_rw=path_rw, frw_bcwd=false, L=L)
    v2t20 = get_corr(path_data, ens, "88", Gamma_l[5], path_rw=path_rw, frw_bcwd=false, L=L)
    v3t30 = get_corr(path_data, ens, "88", Gamma_l[6], path_rw=path_rw, frw_bcwd=false, L=L)

    if impr
        beta = ens.beta
        if impr_set == "1"
            cv_l = cv_loc(beta)
        elseif impr_set =="2"
            cv_l = cv_loc_set2(beta)
        end 
        improve_corr_vkvk!(v1v1, v1t10, 2*cv_l, std=std)
        improve_corr_vkvk!(v2v2, v2t20, 2*cv_l, std=std)
        improve_corr_vkvk!(v3v3, v3t30, 2*cv_l, std=std)
    end

    g88_ll = Corr(1/6 .* (2 * (v1v1.obs + v2v2.obs + v3v3.obs) /3 ), v1v1.id, "G88ll") 
    if frw_bcwd
        frwd_bckwrd_symm!(g88_ll)
    end

    if cons
        v1v1_c = get_corr(path_data, ens, "88", Gamma_c[1], path_rw=path_rw, frw_bcwd=false, L=L)
        v2v2_c = get_corr(path_data, ens, "88", Gamma_c[2], path_rw=path_rw, frw_bcwd=false, L=L)
        v3v3_c = get_corr(path_data, ens, "88", Gamma_c[3], path_rw=path_rw, frw_bcwd=false, L=L)
    
        v1t10_c = get_corr(path_data, ens, "88", Gamma_c[4], path_rw=path_rw, frw_bcwd=false, L=L)
        v2t20_c = get_corr(path_data, ens, "88", Gamma_c[5], path_rw=path_rw, frw_bcwd=false, L=L)
        v3t30_c = get_corr(path_data, ens, "88", Gamma_c[6], path_rw=path_rw, frw_bcwd=false, L=L)
        
        if impr
            beta = ens.beta
            if impr_set == "1"
                cv_l = cv_loc(beta) 
                cv_c = cv_cons(beta)
            elseif impr_set == "2"
                cv_l = cv_loc_set2(beta)
                cv_c = cv_cons_set2(beta)
            end
            improve_corr_vkvk_cons!(v1v1_c, v1t10, v1t10_c, cv_l, cv_c, std=std)
            improve_corr_vkvk_cons!(v2v2_c, v2t20, v2t20_c, cv_l, cv_c, std=std)
            improve_corr_vkvk_cons!(v3v3_c, v3t30, v3t30_c, cv_l, cv_c, std=std)
        end
        g88_lc = Corr(1/6 .* ( 2 * (v1v1_c.obs + v2v2_c.obs + v3v3_c.obs) /3 ), v1v1.id, "G88lc") 
        if frw_bcwd
            frwd_bckwrd_symm!(g88_lc)
        end
    end 

    !cons ? (return g88_ll) : (return g88_ll, g88_lc) 

end

function corr88(path_data::String, ens::EnsInfo, g33_ll::Corr; g33_lc::Union{Nothing, Corr}=nothing, path_rw::Union{Nothing,String}=nothing, L::Int64=1, frw_bcwd::Bool=true, impr::Bool=true, impr_set::String="1", cons::Bool=true, std::Bool=true)

    if cons
        ll88_conn, lc88_conn = corr88_conn(path_data, ens, g33_ll, g33_lc=g33_lc, path_rw=path_rw, L=L, frw_bcwd=frw_bcwd, impr=impr, impr_set=impr_set, cons=cons, std=std)
        ll88_disc, lc88_disc = corr88_disc(path_data, ens, path_rw=path_rw, L=L, frw_bcwd=frw_bcwd, impr=impr, impr_set=impr_set, cons=cons, std=std)

        T = length(ll88_conn.obs)
        println(T, " ", length(ll88_disc.obs))
        ll88_tot = Corr(ll88_conn.obs + ll88_disc.obs[1:T], ll88_conn.id, "G88ll")
        lc88_tot = Corr(lc88_conn.obs + lc88_disc.obs[1:T], lc88_conn.id, "G88lc")
        return ll88_tot, lc88_tot
    else
        ll88_conn = corr88_conn(path_data, ens, g33_ll, g33_lc=g33_lc, path_rw=path_rw, L=L, frw_bcwd=frw_bcwd, impr=impr, impr_set=impr_set, cons=cons, std=std)
        ll88_disc = corr88_disc(path_data, ens, path_rw=path_rw, L=L, frw_bcwd=frw_bcwd, impr=impr, impr_set=impr_set, cons=cons, std=std)
        T = length(ll88_conn.obs)
        println(T, " ", length(ll88_disc.obs))
        
        ll88_tot = Corr(ll88_conn.obs + ll88_disc.obs[1:T], ll88_conn.id, "G88ll")
        return ll88_tot
    end
end


@doc raw"""
corr08_conn(g33_ll::Corr, g88_ll::Corr; g33_lc::Union{Corr, Nothing}=nothing, g88_lc::Union{Corr, Nothing}=nothing)

This function return the local-loca G08 connected correlator from the G33 and G88 connected according to:

```math
G_{08}^{conn} = \frac{\sqrt{3}}{2}(G_{33} - G_{88}^{conn})
```

If the local-conserved G33_lc and G88_lc correlators are passed, this functions also returns the local-consereved G08 connected correlator.

Examples:
```@example
g08_ll, g08_lc = corr08_conn(g33, g88, g33_lc=g33_lc, g88_lc=g88_lc) 
```
"""
function corr08_conn(g33_ll::Corr, g88_ll::Corr; g33_lc::Union{Corr, Nothing}=nothing, g88_lc::Union{Corr, Nothing}=nothing)
    
    g08_ll =  Corr(sqrt(3)/2 * (g33_ll.obs - g88_ll.obs), g33_ll.id, "G08ll")
    if !isnothing(g33_lc) && !isnothing(g88_lc)
        g08_lc =  Corr(sqrt(3)/2 * (g33_lc.obs - g88_lc.obs), g33_lc.id, "G08lc")
        return g08_ll, g08_lc
    else
        return g08_ll
    end
end