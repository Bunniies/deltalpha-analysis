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

include("const.jl")
include("data_management.jl")


const path_data = "/Users/alessandroconigli/Lattice/data/HVP/2ptdata"
const path_rw   = "/Users/alessandroconigli/Lattice/data/HVP/rwf_deflated"
const path_store_bdio = "/Users/alessandroconigli/MyDrive/postdoc-mainz/projects/deltalpha/corr"

const IMPR = true
const STD_DERIV = false
const RENORM = true


# enslist = ["H101"]
enslist = sort([ "H101", "H102", "C101", "C102", "D150",
         "B450", "N451", "D450", "D451", "D452",
         "N202", "N203", "N200", "D200", "D201", "E250",
         "N300", "J303", "E300",
         "J500", "J501"
         ])
ensinfo = EnsInfo.(enslist)

##
#============ OBSERVABLE ALLOCATIONS ============#

g33_ll = Vector{Corr}(undef, 1)
g33_lc = Vector{Corr}(undef, 1)

g88_ll_conn = Vector{Corr}(undef, 1)
g88_lc_conn = Vector{Corr}(undef, 1)

g88_ll_disc = Vector{Corr}(undef, 1)
g88_lc_disc = Vector{Corr}(undef, 1)

gdelta_iso_ll = Vector{Corr}(undef, 1) # (G33 - G88) ll
gdelta_iso_lc = Vector{Corr}(undef, 1) # (G33 - G88) lc

g08_ll_conn = Vector{Corr}(undef, 1)
g08_lc_conn = Vector{Corr}(undef, 1)

g08_ll_disc = Vector{Corr}(undef, 1)
g08_lc_disc = Vector{Corr}(undef, 1)

g08_ll = Vector{Corr}(undef, 1)
g08_lc = Vector{Corr}(undef, 1)

gcc_ll_conn = Vector{Corr}(undef, 1)
gcc_lc_conn = Vector{Corr}(undef, 1)

gcc_ll_disc = Vector{Corr}(undef, 1)
gcc_lc_disc = Vector{Corr}(undef, 1)
gcc_cc_disc = Vector{Corr}(undef, 1)

gc8_ll_disc = Vector{Corr}(undef, 1)
gc8_lc_disc = Vector{Corr}(undef, 1)
gc8_cc_disc = Vector{Corr}(undef, 1)

@time begin
    for (k, ens) in enumerate(ensinfo)
        @info("Reading data ensemble: $(ens.id)")

        for impr_set in ["1", "2"]
            println("   - Impr Set: ", impr_set)
            println("        - G33 ll and lc correlator")
            g33_ll[1], g33_lc[1] = corr33(path_data, ens, sector="light", path_rw=path_rw, impr=IMPR, impr_set=impr_set, cons=true, frw_bcwd=true, std=STD_DERIV)

            println("        - G88 connected ll and lc correlator")        
            g88_ll_conn[1], g88_lc_conn[1] = corr88_conn(path_data, ens, g33_ll[1], g33_lc=g33_lc[1], path_rw=path_rw, impr=IMPR, impr_set=impr_set, cons=true, frw_bcwd=true, std=STD_DERIV)
            
            println("        - G08 connected ll and lc correlator")
            g08_ll_conn[1], g08_lc_conn[1] = corr08_conn(g33_ll[1], g88_ll_conn[1], g33_lc=g33_lc[1], g88_lc=g88_lc_conn[1])

            println("        - Gcc connected ll and lc correlator")
            try
                gcc_ll_conn[1], gcc_lc_conn[1] = corrcc_conn(path_data, ens, path_rw=path_rw, impr=IMPR, impr_set=impr_set, cons=true, frw_bcwd=true, std=STD_DERIV)
            catch
                println("        - cc connected not found for ens $(ens.id)")
                T = HVPobs.Data.get_T(ens.id)
                gcc_ll_conn[1] = gcc_lc_conn[1] = Corr(fill(uwreal(0.0), T), ens.id, "Gcc") 
            end

            println("        - G88, Gcc, G08 disconnected")
            if ens.kappa_l == ens.kappa_s
                println("Symmetric point ensemble") 
                T = HVPobs.Data.get_T(ens.id)
                g88_ll_disc[1] = g88_lc_disc[1] = Corr(fill(uwreal(0.0), T), ens.id, "G88disc")  
                gcc_ll_disc[1] = gcc_lc_disc[1] = gcc_cc_disc[1] = Corr(fill(uwreal(0.0), T), ens.id, "Gccdisc")
                g08_ll_disc[1] = g08_lc_disc[1] = Corr(fill(uwreal(0.0), T), ens.id, "G08disc")
                gc8_ll_disc[1] = gc8_lc_disc[1] = gc8_cc_disc[1] = Corr(fill(uwreal(0.0), T), ens.id, "Gc8disc")

            else
                try    
                    g88_ll_disc[1], g88_lc_disc[1] = corrDisconnected(path_data, ens, "88", path_rw=path_rw, impr=IMPR, impr_set=impr_set,  frw_bcwd=true, cc=false, std=STD_DERIV)
                    gcc_ll_disc[1], gcc_lc_disc[1], gcc_cc_disc[1] = corrDisconnected(path_data, ens, "cc", path_rw=path_rw, impr=IMPR, impr_set=impr_set, cc=true, frw_bcwd=true, std=STD_DERIV)
                    g08_ll_disc[1], g08_lc_disc[1] = corrDisconnected(path_data, ens, "08", path_rw=path_rw, impr=IMPR, impr_set=impr_set, frw_bcwd=true, cc=false, std=STD_DERIV)
                    gc8_ll_disc[1], gc8_lc_disc[1], gc8_cc_disc[1] = corrDisconnected(path_data, ens, "c8", path_rw=path_rw, impr=IMPR, impr_set=impr_set, cc=true, frw_bcwd=true, std=STD_DERIV)
                catch
                    println("disconnected failed")
                    T = HVPobs.Data.get_T(ens.id)
                    g88_ll_disc[1] = g88_lc_disc[1] = gcc_ll_disc[1] =  gcc_lc_disc[1] = g08_ll_disc[1] = g08_lc_disc[1] = gc8_ll_disc[1] = gc8_lc_disc[1] = gcc_cc_disc[1] =  gc8_cc_disc[1] = Corr(fill(uwreal(0.0), T), ens.id, "Gdisc")
                end
            end

            if RENORM

                Z3 = get_Z3(ens, impr_set=impr_set)
                renormalize!(g33_ll[1], Z3^2)
                renormalize!(g33_lc[1], Z3)
                
                Z8 = get_Z8(ens, impr_set=impr_set)
                renormalize!(g88_ll_conn[1], Z8^2)
                renormalize!(g88_lc_conn[1], Z8)
                renormalize!(g88_ll_disc[1], Z8^2)
                renormalize!(g88_lc_disc[1], Z8)
    
                Zcc = Zvc_l[ens.id]
                renormalize!(gcc_ll_conn[1], Zcc*Zcc)
                renormalize!(gcc_lc_conn[1], Zcc)
                renormalize!(gcc_ll_disc[1], Zcc*Zcc)
                renormalize!(gcc_lc_disc[1], Zcc)
    
                Z08 = get_Z08(ens, impr_set=impr_set)
                renormalize!(g08_ll_conn[1], Z8*Z08)
                renormalize!(g08_lc_conn[1], Z08)
                renormalize!(g08_ll_disc[1], Z8*Z08)
                renormalize!(g08_lc_disc[1], Z08)
            end
            
            gdelta_iso_ll[1] = Corr(-g88_ll_conn[1].obs - g88_ll_disc[1].obs + g33_ll[1].obs, g33_ll[1].id, "deltaiso_ll")
            gdelta_iso_lc[1] = Corr(-g88_lc_conn[1].obs - g88_lc_disc[1].obs + g33_lc[1].obs, g33_lc[1].id, "deltaiso_lc")

            g08_ll[1] = Corr(g08_ll_conn[1].obs + g08_ll_disc[1].obs, g08_ll_conn[1].id, "G08ll")
            g08_lc[1] = Corr(g08_lc_conn[1].obs + g08_lc_disc[1].obs, g08_lc_conn[1].id, "G08lc")


            data_corr = Dict{String, Array{uwreal}}(
                "g33_ll" => g33_ll[1].obs,
                "g33_lc" => g33_lc[1].obs
            )
            
            data_corr["g3388_dlt_ll"] = gdelta_iso_ll[1].obs
            data_corr["g3388_dlt_lc"] = gdelta_iso_lc[1].obs

            data_corr["gcc_ll_conn"] = gcc_ll_conn[1].obs
            data_corr["gcc_lc_conn"] = gcc_lc_conn[1].obs

            data_corr["gcc_ll_disc"] = gcc_ll_disc[1].obs
            data_corr["gcc_lc_disc"] = gcc_lc_disc[1].obs
            data_corr["gcc_cc_disc"] = gcc_cc_disc[1].obs

            data_corr["g08_ll"] = g08_ll[1].obs
            data_corr["g08_lc"] = g08_lc[1].obs
 
            #data_corr["gc8_ll_disc"] = gc8_ll_disc[1].obs
            #data_corr["gc8_lc_disc"] = gc8_lc_disc[1].obs
            #data_corr["gc8_cc_disc"] = gc8_cc_disc[1].obs
            
            io = IOBuffer()
            write(io, "$(ens.id)  HVP correlators, improvement set $(impr_set)")
            fb = ALPHAdobs_create(joinpath(path_store_bdio,"$(ens.id)_corr_set$(impr_set)"), io)

            extra = Dict{String, Any}("Ens" => "H101", "Impr_Set" => impr_set)
            ALPHAdobs_write(fb, data_corr, extra=extra)
            
            ALPHAdobs_close(fb)

        end # end impr_set loop
    end # end ens loop

end # end time begin


## TEST WITH READING

fb = BDIO_open(joinpath(path_store_bdio, "H101_corr_set1"), "r")

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