using Revise
using HVPobs
using PyPlot, LaTeXStrings
using OrderedCollections
using TimerOutputs
using BDIO, ADerrors
import ADerrors: err

const path_data = "/Users/alessandroconigli/Lattice/data/HVP/2ptdata"
const path_rw   = "/Users/alessandroconigli/Lattice/data/HVP/rwf_deflated"

ensinfo = EnsInfo("H101")


## check bare vector current
v1v1 = get_corr(path_data, ensinfo, "light", "V1V1", path_rw=path_rw, frw_bcwd=false, L=1)
v2v2 = get_corr(path_data, ensinfo, "light", "V2V2", path_rw=path_rw, frw_bcwd=false, L=1)
v3v3 = get_corr(path_data, ensinfo, "light", "V3V3", path_rw=path_rw, frw_bcwd=false, L=1)


gvv_bare = Corr(abs.(v1v1.obs .+ v2v2.obs .+ v3v3.obs)./3, v1v1.id, "Gvvbare")
frwd_bckwrd_symm!(gvv_bare)
uwerr.(gvv_bare.obs)
for i in 1:9
    println(i, " ", value(gvv_bare.obs[i+1]), " ", err(gvv_bare.obs[i+1]) )
end
# agree!

## check bare improved derivative of the vector-tensor current

v1t10 = get_corr(path_data, ensinfo, "light", "V1T10", path_rw=path_rw, frw_bcwd=false, L=1)
v2t20 = get_corr(path_data, ensinfo, "light", "V2T20", path_rw=path_rw, frw_bcwd=false, L=1)
v3t30 = get_corr(path_data, ensinfo, "light", "V3T30", path_rw=path_rw, frw_bcwd=false, L=1)

# frwd_bckwrd_antisymm!(v1t10)
# frwd_bckwrd_antisymm!(v2t20)
# frwd_bckwrd_antisymm!(v3t30)
vt = Corr((v1t10.obs .+ v2t20.obs .+ v3t30.obs)./3, v1t10.id, "Gvt")
frwd_bckwrd_antisymm!(vt)
deriv_vt= HVPobs.Obs.improve_derivative(vt.obs, std=false)
uwerr.(deriv_vt)
for i in 1:9
    println(i, " ", value(deriv_vt[i]), " ", err(deriv_vt[i]) )
end
# agree! when following this order of doing things. First improving and then symmetrizing does not work. Also agree with standard central derivatives.

## check  and improved VV correlator 
v1v1 = get_corr(path_data, ensinfo, "light", "V1V1", path_rw=path_rw, frw_bcwd=false, L=1)
v2v2 = get_corr(path_data, ensinfo, "light", "V2V2", path_rw=path_rw, frw_bcwd=false, L=1)
v3v3 = get_corr(path_data, ensinfo, "light", "V3V3", path_rw=path_rw, frw_bcwd=false, L=1)

gvv = Corr(abs.(v1v1.obs .+ v2v2.obs .+ v3v3.obs)./3, v1v1.id, "Gvvbare")

v1t10 = get_corr(path_data, ensinfo, "light", "V1T10", path_rw=path_rw, frw_bcwd=false, L=1)
v2t20 = get_corr(path_data, ensinfo, "light", "V2T20", path_rw=path_rw, frw_bcwd=false, L=1)
v3t30 = get_corr(path_data, ensinfo, "light", "V3T30", path_rw=path_rw, frw_bcwd=false, L=1)

vt = Corr((v1t10.obs .+ v2t20.obs .+ v3t30.obs)./3, v1t10.id, "Gvt")

cv_l = cv_loc(ensinfo.beta)
improve_corr_vkvk!(gvv, vt, -2*cv_l, std=false)
frwd_bckwrd_symm!(gvv)

uwerr.(gvv.obs)
for i in 1:9
    println(i, " ", value(gvv.obs[i+1]), " ", err(gvv.obs[i+1]) )
end
## agree both with improved and standard derivative!

# new function corrConnected - already included in HPVobs.Automation
@doc raw"""
    corrConnected(path_data::String, ens::EnsInfo, sector::String; path_rw=Union{Nothing,String}=nothing, impr::Bool=true, impr_set::String="1", std::Bool=false, L::Int64=1)

Computes the Vector-Vector  local-local and local-conserved correlators for a given EnsInfo `ens` and sector `sector`.
The supported sectors are: light, strange, charm, charm_plus.
It returns the forward-backward symmetrised local-local and local-conserved correlators as Corr objects.
No charge factor is included in the computation. 
Optional flags:    

    - path_rw  : if provided, correlators are reweighted.  
    - impr     : if true, the correlator is improved.  
    - impr_set : either "1" or "2", to select set1 (Mainz) or set2 (ALPHA) improvement coefficients accordingly.  
    - std      : if true, standard symmetric derivatives are used for the vector-tensor correlator. If false, improved derivatives are used.
    - L        : correlators are normalised with the volume L^3. L=1 by default.    

Examples:
```@example
ens = EnsInfo("H101")
gvv_ll, gvv_lc = corrConnected(pathToData, ens, "ligth", path_rw=pathToRwf, impr=true, impr_set="1", std=false, L=1) 
```
"""
function corrConnected(path_data::String, ens::EnsInfo, sector::String; path_rw=Union{Nothing,String}=nothing, impr::Bool=true, impr_set::String="1", std::Bool=false, L::Int64=1)
    if sector ∉ ["light", "strange", "charm", "charm_plus"]
        error("Flavour sector '$(sector)' not recognised.")
    end
    # local-local VV
    Gamma_l = ["V1V1", "V2V2", "V3V3", "V1T10", "V2T20", "V3T30"]
    
    v1v1 = get_corr(path_data, ens, sector, Gamma_l[1], path_rw=path_rw, frw_bcwd=false, L=L)
    v2v2 = get_corr(path_data, ens, sector, Gamma_l[2], path_rw=path_rw, frw_bcwd=false, L=L)
    v3v3 = get_corr(path_data, ens, sector, Gamma_l[3], path_rw=path_rw, frw_bcwd=false, L=L)
    gvv_ll = Corr(abs.(v1v1.obs .+ v2v2.obs .+ v3v3.obs)./3, ens.id, "Gvv_ll_"*sector)
        
    # local-conserved VV
    Gamma_c = ["V1V1c", "V2V2c", "V3V3c", "V1cT10", "V2cT20", "V3cT30"]
    v1v1_c = get_corr(path_data, ens, sector, Gamma_c[1], path_rw=path_rw, frw_bcwd=false, L=L)
    v2v2_c = get_corr(path_data, ens, sector, Gamma_c[2], path_rw=path_rw, frw_bcwd=false, L=L)
    v3v3_c = get_corr(path_data, ens, sector, Gamma_c[3], path_rw=path_rw, frw_bcwd=false, L=L)
    gvv_lc = Corr(abs.(v1v1_c.obs .+ v2v2_c.obs .+ v3v3_c.obs)./3, ens.id, "Gvv_lc_"*sector)

    
    if impr
        # local-local VT
        v1t10 = get_corr(path_data, ens, sector, Gamma_l[4], path_rw=path_rw, frw_bcwd=false, L=L)
        v2t20 = get_corr(path_data, ens, sector, Gamma_l[5], path_rw=path_rw, frw_bcwd=false, L=L)
        v3t30 = get_corr(path_data, ens, sector, Gamma_l[6], path_rw=path_rw, frw_bcwd=false, L=L)
        gvt_ll = Corr((v1t10.obs .+ v2t20.obs .+ v3t30.obs)./3, ens.id, "Gvt_ll_"*sector)

        # local-consereved VT
        v1t10_c = get_corr(path_data, ens, sector, Gamma_c[4], path_rw=path_rw, frw_bcwd=false, L=L)
        v2t20_c = get_corr(path_data, ens, sector, Gamma_c[5], path_rw=path_rw, frw_bcwd=false, L=L)    
        v3t30_c = get_corr(path_data, ens, sector, Gamma_c[6], path_rw=path_rw, frw_bcwd=false, L=L)
        gvt_lc = Corr((v1t10_c.obs .+ v2t20_c.obs .+ v3t30_c.obs)./3, ens.id, "Gvt_lc_"*sector)

        beta = ens.beta
        if impr_set == "1"
            cv_l = cv_loc(beta)
            cv_c = cv_cons(beta)
        elseif impr_set =="2"
            cv_l = cv_loc_set2(beta)
            cv_c = cv_cons_set2(beta)
        end

        improve_corr_vkvk!(gvv_ll, gvt_ll, -2*cv_l, std=std)
        improve_corr_vkvk_cons!(gvv_lc, gvt_ll, gvt_lc, -cv_l, -cv_c, std=std)
    end
    frwd_bckwrd_symm!(gvv_ll)
    frwd_bckwrd_symm!(gvv_lc)

    return gvv_ll, gvv_lc
end

gvv_ll, gvv_lc = corrConnected(path_data, ensinfo, "light", path_rw=path_rw, impr=true, impr_set="1", std=true)

# testing with std deriv G33 for H101
aa = gvv_ll.obs ./2; uwerr.(aa)
bb = gvv_lc.obs ./2; uwerr.(bb)

##

@doc raw"""
    corrDisconnected(path_data::String, ens::EnsInfo, fl::String; path_rw::Union{Nothing,String}=nothing, impr::Bool=true, impr_set::String="1", std::Bool=true, L::Int64=1)

Computes the Vector-Vector  local-local, local-conserved and consereved-conserved disconnected correlators for a given EnsInfo `ens` and flavour `fl`.
The supported flavours are: [08, 0c, 80, 88, 8c, c0, c8, cc].
It returns the forward-backward symmetrised local-local, local-conserved, conserved-conserved correlators as Corr objects.
No charge factor is included in the computation. 

Optional flags:    

    - path_rw  : if !isnothing(path_rw), correlators are reweighted.  
    - impr     : if true, the correlator is improved.  
    - impr_set : either "1" or "2", to select set1 (Mainz) or set2 (ALPHA) improvement coefficients accordingly.  
    - std      : if true, standard symmetric derivatives are used for the vector-tensor correlator. If false, improved derivatives are used.
    - L        : correlators are normalised with the volume L^3. L=1 by default.    

Examples:
```@example
ens = EnsInfo("N200")
gvv_ll, gvv_lc , gvv_cc = corrDonnected(pathToData, ens, "cc", path_rw=pathToRwf, impr=true, impr_set="1", std=false, L=1) 
```
"""
function corrDisconnected(path_data::String, ens::EnsInfo, fl::String; path_rw::Union{Nothing,String}=nothing, impr::Bool=true, impr_set::String="1", std::Bool=true, L::Int64=1)
    if fl ∉ ["08", "0c", "80", "88", "8c", "c0", "c8", "cc"]
        error("Unrecognised flavour structure $(fl): choose from [08, 0c, 80, 88, 8c, c0, c8, cc]")
    end

    corr_dict = get_corr_disc(path_data, ens, fl, path_rw=path_rw, frw_bcwd=false, L=L)
    g_ll = Corr(corr_dict["VV"].obs, ens.id,   "G"*fl*"_ll_disc")
    g_lc = Corr(corr_dict["VVc"].obs, ens.id,  "G"*fl*"_lc_disc")
    g_cc = Corr(corr_dict["VcVc"].obs, ens.id, "G"*fl*"_cc_disc")

    if impr
        beta = ens.beta
        if impr_set == "1"
            cv_l = cv_loc(beta) 
            cv_c = cv_cons(beta)      
        elseif impr_set =="2"
            cv_l = cv_loc_set2(beta)
            cv_c = cv_cons_set2(beta)
        end
        improve_corr_vkvk!(g_ll, corr_dict["VT"], 2*cv_l, std=std)
        improve_corr_vkvk!(g_cc, corr_dict["VcT"], 2*cv_c, std=std)
        improve_corr_vkvk_cons!(g_lc, corr_dict["VT"], corr_dict["VcT"], cv_l, cv_c, std=std)
    end
   
    frwd_bckwrd_symm!(g_ll)
    frwd_bckwrd_symm!(g_lc)
    frwd_bckwrd_symm!(g_cc)

    return  g_ll, g_lc, g_cc
end
##
aa, bb, cc = corrDisconnected(path_data, EnsInfo("N200"), "cc", path_rw=path_rw)
