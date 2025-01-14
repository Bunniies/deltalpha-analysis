######################################################################################
# This file was created by Alessandro Conigli - 03/2024
# Here we compute the correlators 
# The code works for a general kernel function KRNL to be used in the TMR integral.
# It also allows for two different windows (SD and ILD) to compute the isovector contribution.
# This splitting is required to correctly apply multiplicative tree-level improvement on the isovector piece.
# The isoscalar contribution is computed with the standard non-subtracted kernel 
# The observables are then stored in a BDIO file, foll owing the order 
########################################################################################

using Revise
using HVPobs
using PyPlot, LaTeXStrings
using OrderedCollections
using TimerOutputs
using BDIO, ADerrors
import ADerrors: err
using ALPHAio

include("../utils/const.jl")


const path_data = "/Users/alessandroconigli/Lattice/data/HVP/2ptdata"
const path_rw   = "/Users/alessandroconigli/Lattice/data/HVP/rwf_deflated"
const path_store_bdio = "/Users/alessandroconigli/MyDrive/postdoc-mainz/projects/deltalpha/corr/impr_deriv/"

const IMPR = true
const STD_DERIV = false
const RENORM = true
println("STD_DERIV: ", STD_DERIV)

# enslist = ["H101"]
#ensinfo = EnsInfo.(enslist)
#enslist = sort([ #"H101", "H102", "N101", "C101", "C102", "D150",
        #"B450", "N451", "D450", "D451", "D452"])
        # "N202", "N203", "N200", "D251", "D200", "D201", "E250"])
        # "N300", "J303", "J304", "E300",
        # "J500", "J501"])

# enslist = sort(["A653"])        
enslist = sort(["D450"])        
ensinfo = EnsInfo.(enslist)
##
#============ OBSERVABLE ALLOCATIONS ============#

# isovector
g33_ll = Vector{Corr}(undef, 1)
g33_lc = Vector{Corr}(undef, 1)

# isoscalar
g88_ll_conn = Vector{Corr}(undef, 1)
g88_lc_conn = Vector{Corr}(undef, 1)

g88_ll_disc = Vector{Corr}(undef, 1)
g88_lc_disc = Vector{Corr}(undef, 1)

gdelta_iso_ll = Vector{Corr}(undef, 1) # (G33 - G88) ll
gdelta_iso_lc = Vector{Corr}(undef, 1) # (G33 - G88) lc

# isoscalar singlet
g08_ll_conn = Vector{Corr}(undef, 1)
g08_lc_conn = Vector{Corr}(undef, 1)
g00_lc_conn = Vector{Corr}(undef, 1)

g08_ll_disc = Vector{Corr}(undef, 1)
g08_lc_disc = Vector{Corr}(undef, 1)

g08_ll = Vector{Corr}(undef, 1)
g08_lc = Vector{Corr}(undef, 1)

# charm
gcc_ll_conn = Vector{Corr}(undef, 1)
gcc_lc_conn = Vector{Corr}(undef, 1)
gcc_ll_conn_plus = Vector{Corr}(undef, 1)
gcc_lc_conn_plus = Vector{Corr}(undef, 1)

gcc_ll_disc = Vector{Corr}(undef, 1)
gcc_lc_disc = Vector{Corr}(undef, 1)
gcc_cc_disc = Vector{Corr}(undef, 1)

gc8_ll_disc = Vector{Corr}(undef, 1)
gc8_lc_disc = Vector{Corr}(undef, 1)
gc8_cc_disc = Vector{Corr}(undef, 1)

@time begin
    for (k, ens) in enumerate(ensinfo)
        @info("Reading data ensemble: $(ens.id)")

        for impr_set in ["2"] #["1", "2"]
            println("   - Impr Set: ", impr_set)

            println("        - G33 ll and lc correlator")
            gll_ll, gll_lc = corrConnected(path_data, ens, "light", path_rw=path_rw, impr=IMPR, impr_set=impr_set, std=STD_DERIV)
            
            g33_ll[1] = Corr(0.5*gll_ll.obs, ens.id, "G33_ll")
            g33_lc[1] = Corr(0.5*gll_lc.obs, ens.id, "G33_lc")

            println("        - G88 connected ll and lc correlator")  
            gss_ll, gss_lc = corrConnected(path_data, ens, "strange", path_rw=path_rw, impr=IMPR, impr_set=impr_set, std=STD_DERIV)
            
            g88_ll_conn[1]  = Corr(1/6 .*( gll_ll.obs + 2*gss_ll.obs), ens.id, "G88_ll_conn" )
            g88_lc_conn[1]  = Corr(1/6 .*( gll_lc.obs + 2*gss_lc.obs), ens.id, "G88_lc_conn" )
                        
            println("        - G08 connected ll and lc correlator")
            g08_ll_conn[1] = Corr(1/(2*sqrt(3)) .* (gll_ll.obs .- gss_ll.obs), ens.id,  "G08_ll_conn")
            g08_lc_conn[1] = Corr(1/(2*sqrt(3)) .* (gll_lc.obs .- gss_lc.obs), ens.id,  "G08_lc_conn")

            println("        - G08 connected ll and lc correlator")
            g00_lc_conn[1] = Corr(1/4 * (2*gll_lc.obs .+ gss_lc.obs), ens.id, "G00_lc_conn")

            println("        - Gcc connected ll and lc correlator")
            try
                gcc_ll_conn_plus[1], gcc_lc_conn_plus[1] = corrConnected(path_data, ens, "charm_plus", path_rw=path_rw, impr=IMPR, impr_set=impr_set, std=STD_DERIV) 
                gcc_ll_conn[1], gcc_lc_conn[1] = corrConnected(path_data, ens, "charm", path_rw=path_rw, impr=IMPR, impr_set=impr_set, std=STD_DERIV) 
            catch
                println("        - cc connected not found for ens $(ens.id)")
                T = HVPobs.Data.get_T(ens.id)
                gcc_ll_conn[1] = gcc_lc_conn[1] = Corr(fill(uwreal(0.0), T), ens.id, "Gcc") 
                gcc_ll_conn_plus[1] = gcc_lc_conn_plus[1] = Corr(fill(uwreal(0.0), T), ens.id, "Gcc") 
            end

            if ens.kappa_l != ens.kappa_s
                println("        - G88, G08, Gcc, Gc8 disconnected")
                #try    
                    g88_ll_disc[1], g88_lc_disc[1], _ = corrDisconnected(path_data, ens, "88", path_rw=path_rw, impr=IMPR, impr_set=impr_set, std=STD_DERIV)
                    g08_ll_disc[1], g08_lc_disc[1], _ = corrDisconnected(path_data, ens, "80", path_rw=path_rw, impr=IMPR, impr_set=impr_set, std=STD_DERIV)
                    
                    gcc_ll_disc[1], gcc_lc_disc[1], gcc_cc_disc[1] = corrDisconnected(path_data, ens, "cc", path_rw=path_rw, impr=false)
                    gc8_ll_disc[1], gc8_lc_disc[1], gc8_cc_disc[1] = corrDisconnected(path_data, ens, "c8", path_rw=path_rw, impr=false)
                #catch
                    println("disconnected failed")
                    T = HVPobs.Data.get_T(ens.id)
                    g88_ll_disc[1] = g88_lc_disc[1] = gcc_ll_disc[1] =  gcc_lc_disc[1] = g08_ll_disc[1] = g08_lc_disc[1] = gc8_ll_disc[1] = gc8_lc_disc[1] = gcc_cc_disc[1] =  gc8_cc_disc[1] = Corr(fill(uwreal(0.0), T), ens.id, "Gdisc")
                 end
            else
                T = HVPobs.Data.get_T(ens.id)
                g88_ll_disc[1] = g88_lc_disc[1] = gcc_ll_disc[1] =  gcc_lc_disc[1] = g08_ll_disc[1] = g08_lc_disc[1] = gc8_ll_disc[1] = gc8_lc_disc[1] = gcc_cc_disc[1] =  gc8_cc_disc[1] = Corr(fill(uwreal(0.0), T), ens.id, "Gdisc")
            end

            if RENORM
                # isovector
                Z3 = get_Z3(ens, impr_set=impr_set) 
                renormalize!(g33_ll[1], Z3^2)  
                renormalize!(g33_lc[1], Z3)    
                
                # isoscalaar
                Z8 = get_Z8(ens, impr_set=impr_set)
                Z08 = get_Z08(ens, impr_set=impr_set)
                g88_ll_conn[1].obs[:] = Z8^2 .* g88_ll_conn[1].obs[:] .+ 2 .* Z8 * Z08 .* g08_ll_conn[1].obs[:]
                g88_lc_conn[1].obs[:] = Z8 .* g88_lc_conn[1].obs[:] .+ Z08 .* g08_lc_conn[1].obs[:]
                
                g88_ll_disc[1].obs[:] = Z8^2 .* g88_ll_disc[1].obs[:] .+ 2 .* Z8 * Z08 .* g08_ll_disc[1].obs[:] 
                g88_lc_disc[1].obs[:] = Z8 .* g88_lc_disc[1].obs[:] .+   Z08 .* g08_lc_disc[1].obs[:] 

                # charm 
                Zcc = Zvc_l[ens.id]
                renormalize!(gcc_ll_conn[1], Zcc*Zcc)
                renormalize!(gcc_lc_conn[1], Zcc)
                renormalize!(gcc_ll_conn_plus[1], Zcc*Zcc)
                renormalize!(gcc_lc_conn_plus[1], Zcc)
                
                # isoscalar isosinglet
                g08_lc_conn[1].obs[:] = Z8 .* g08_lc_conn[1].obs[:] .+ Z08 .* g00_lc_conn[1].obs[:]
                g08_lc_disc[1].obs[:] = Z8 .* g08_lc_disc[1].obs[:] # here still missing the g00_lc_disc !

            end
            
            if ens.kappa_l != ens.kappa_s
                gdelta_iso_ll[1] = Corr(-g88_ll_conn[1].obs - g88_ll_disc[1].obs + g33_ll[1].obs, g33_ll[1].id, "deltaiso_ll")
                gdelta_iso_lc[1] = Corr(-g88_lc_conn[1].obs - g88_lc_disc[1].obs + g33_lc[1].obs, g33_lc[1].id, "deltaiso_lc")
                g08_lc[1] = Corr(g08_lc_conn[1].obs + g08_lc_disc[1].obs, g08_lc_conn[1].id, "G08lc")

            else
                gdelta_iso_ll[1] = Corr(-g88_ll_conn[1].obs + g33_ll[1].obs, g33_ll[1].id, "deltaiso_ll")
                gdelta_iso_lc[1] = Corr(-g88_lc_conn[1].obs + g33_lc[1].obs, g33_lc[1].id, "deltaiso_lc")
                g08_lc[1] = Corr(g08_lc_conn[1].obs , g08_lc_conn[1].id, "G08lc")
            end
            


            data_corr = Dict{String, Array{uwreal}}()
            
            data_corr["g33_ll"] = g33_ll[1].obs
            data_corr["g33_lc"] = g33_lc[1].obs

            data_corr["g88_ll"] = g88_ll_conn[1].obs .+ g88_ll_disc[1].obs
            data_corr["g88_lc"] = g88_lc_conn[1].obs .+ g88_lc_disc[1].obs

            data_corr["g3388_dlt_ll"] = gdelta_iso_ll[1].obs
            data_corr["g3388_dlt_lc"] = gdelta_iso_lc[1].obs 

            data_corr["gcc_ll_conn"] = gcc_ll_conn[1].obs
            data_corr["gcc_lc_conn"] = gcc_lc_conn[1].obs
            data_corr["gcc_ll_conn_plus"] = gcc_ll_conn_plus[1].obs
            data_corr["gcc_lc_conn_plus"] = gcc_lc_conn_plus[1].obs
            
            data_corr["g08_lc"] = g08_lc[1].obs

            data_corr["gcc_cc_disc"] = gcc_cc_disc[1].obs
            data_corr["gc8_cc_disc"] = gc8_cc_disc[1].obs
            
            io = IOBuffer()
            write(io, "$(ens.id)  HVP correlators, improvement set $(impr_set)")
            fb = ALPHAdobs_create(joinpath(path_store_bdio,"$(ens.id)_corr_set$(impr_set)"), io)

            extra = Dict{String, Any}("Ens" => ens.id, "Impr_Set" => impr_set)
            ALPHAdobs_write(fb, data_corr, extra=extra)
            
            ALPHAdobs_close(fb)

        end # end impr_set loop
    end # end ens loop

end # end time begin


## TEST WITH READING

fb = BDIO_open(joinpath(path_store_bdio,  "H102_corr_set2"), "r")

# res = Dict()
count=0
res = Dict()
while ALPHAdobs_next_p(fb)
    count+=1
    println("count:", count)
    d = ALPHAdobs_read_parameters(fb)
    nobs = d["nobs"]
    println("nobs: ", nobs)
    dims = d["dimensions"]
    println("dims: ", dims)
    sz = tuple(d["size"]...)
    println("size:", sz)
    ks = collect(d["keys"])
    # println(ks)
    # println(d["extra"])
    res = ALPHAdobs_read_next(fb, size=sz, keys=ks)
end
BDIO_close!(fb)