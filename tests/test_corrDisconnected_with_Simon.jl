using Revise
using HVPobs
using PyPlot, LaTeXStrings
using OrderedCollections
using TimerOutputs
using BDIO, ADerrors
import ADerrors: err

const path_data = "/Users/alessandroconigli/Lattice/data/HVP/2ptdata"
const path_rw   = "/Users/alessandroconigli/Lattice/data/HVP/rwf_deflated"

ensinfo = EnsInfo("E250")

cdata = get_data_disc(path_data, ensinfo.id, "88")
corr_dict = get_corr_disc(path_data, ensinfo, "88", path_rw=path_rw)

## unimproved and reweighted 88 -> agree with Simon both ll and lc
gll, glc, gcc = corrDisconnected(path_data, ensinfo, "88", path_rw=path_rw, impr=false)
## unimproved but reweighted 08 -> agree with Simon both ll and lc
gll, glc, gcc = corrDisconnected(path_data, ensinfo, "80", path_rw=path_rw, impr=false)
gll1, glc1, gcc1 = corrDisconnected(path_data, ensinfo, "08", path_rw=path_rw, impr=false)


## improved 88 and rweighted -> agree with Simon both ll and lc
gll, glc, gcc = corrDisconnected(path_data, ensinfo, "88", path_rw=path_rw, impr=true, impr_set="1", std=false)

## improved and reweighted 08 -> agree with Simon both ll and lc
gll, glc, gcc = corrDisconnected80(path_data, ensinfo, path_rw=path_rw, impr=true, impr_set="1", std=false)
gll, glc, gcc = corrDisconnected(path_data, ensinfo, "80", path_rw=path_rw, impr=true, impr_set="1", std=false)


## check for the total disconnected correlator

g88_ll_disc_aux, g88_lc_disc_aux, _ = corrDisconnected(path_data, ensinfo, "88", path_rw=path_rw, impr=true, impr_set="1", std=false)

g80_ll_disc_aux, g80_lc_disc_aux, _ = corrDisconnected80(path_data, ensinfo, path_rw=path_rw, impr=true, impr_set="1", std=false)


## VV total and renormalised
Z3 = get_Z3(ensinfo, impr_set="1")
Z8 = get_Z8(ensinfo, impr_set="1")
Z08 = get_Z08(ensinfo, impr_set="1")

VVtot = ( 1/3 * Z8^2 * g88_ll_disc_aux.obs .+ 1/3 * Z8*Z08 * 2*  g80_ll_disc_aux.obs)

VVc_tot =  ( 1/3 * Z8 * g88_lc_disc_aux.obs .+ 1/3 * Z08 * g80_lc_disc_aux.obs)


## renormalised 88 
V88_renorm = Z8 * glc.obs; uwerr.(V88_renorm); V88_renorm

## 
VVc_tot = VVc_tot
uwerr.(VVc_tot)
for t in 1:9
    println(t, " ", value(VVc_tot[t+1]), " ", err(VVc_tot[t+1]))
end

##
uwerr.(glc.obs)
for t in 1:9
    println(t, " ", value(glc.obs[t+1]), " ", err(glc.obs[t+1]))
end