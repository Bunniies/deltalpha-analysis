

######################################################################################
# This file was created by Alessandro Conigli - 04/2025
# Here we compute the correlators 
# This code computes the reconstructed correlator from Nolan's energies and overlap factors
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
path_spec = "/Users/alessandroconigli/Lattice/data/HVP/corr_recostruction"

const IMPR = true
const STD_DERIV = false
const RENORM = true
println("STD_DERIV: ", STD_DERIV)

enslist = sort(["E250"])        
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

        for impr_set in ["1", "2"]
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
                println("           - cc connected not found for ens $(ens.id)")
                T = HVPobs.Data.get_T(ens.id)
                gcc_ll_conn[1] = gcc_lc_conn[1] = Corr(fill(uwreal(0.0), T), ens.id, "Gcc") 
                gcc_ll_conn_plus[1] = gcc_lc_conn_plus[1] = Corr(fill(uwreal(0.0), T), ens.id, "Gcc") 
            end

            if ens.kappa_l != ens.kappa_s
                println("        - G88, G08, Gcc, Gc8 disconnected")
                try    
                    g88_ll_disc[1], g88_lc_disc[1], _ = corrDisconnected(path_data, ens, "88", path_rw=path_rw, impr=IMPR, impr_set=impr_set, std=STD_DERIV)
                    g08_ll_disc[1], g08_lc_disc[1], _ = corrDisconnected(path_data, ens, "80", path_rw=path_rw, impr=IMPR, impr_set=impr_set, std=STD_DERIV)
                    
                    gcc_ll_disc[1], gcc_lc_disc[1], gcc_cc_disc[1] = corrDisconnected(path_data, ens, "cc", path_rw=path_rw, impr=false)
                    gc8_ll_disc[1], gc8_lc_disc[1], gc8_cc_disc[1] = corrDisconnected(path_data, ens, "c8", path_rw=path_rw, impr=false)
                catch
                    println("            - disconnected failed")
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



function reconstr_corr(
    ens::EnsInfo,
    E::Vector{uwreal},Z::Vector{uwreal},Zimpr::Vector{uwreal};
    nmax::Union{Nothing,Int64}=nothing,
    impr_set::String="1",
    IMPR::Bool=true,
    RENORM::Bool=true,
    FreeEN::Bool=false,
    total::Bool=false,
)
    length(Z) != length(Zimpr) ? error("Length of Z and Zimpr are not the same") : nothing
    # nmax chooses the tower of states to be added
    # if nmax is not defined, then the max number of available states are used
    # FreeEN ensures that the last energy level is not used when saturating the corr; usefull for impr BM
    if isnothing(nmax)
        FreeEN ? nmax = minimum([length(E)-1,length(Z)]) : nmax = minimum([length(E),length(Z)])
    else
        if FreeEN
            nmax > length(E)-1 ? error("Input E data not large enough to tower $nmax (-1) states") : nothing
        else
            nmax > length(E) ? error("Input E data not large enough to tower $nmax states") : nothing
        end
        nmax > length(Z) ? error("Input Z data not large enough to tower $nmax states") : nothing
    end

    t = collect(1:Int64(HVPobs.Data.get_T(ens.id))/2+1)
    corrVec_PiPi = [corr_n(1,t.-1,E,Z.^2 ,ens.L)]
    corrVec_Impr = [corr_n(1,t.-1,E,Zimpr.*Z,ens.L)]
    for n=2:nmax
        push!(corrVec_PiPi,corrVec_PiPi[n-1] .+ corr_n(n,t.-1,E,Z.^2 ,ens.L))
        push!(corrVec_Impr,corrVec_Impr[n-1] .+ corr_n(n,t.-1,E,Zimpr.*Z,ens.L))
    end

    # corrVec = deepcopy(corrVec_PiPi)
    if IMPR
        if impr_set == "1"
            cv_l = cv_loc(ens.beta)
        elseif impr_set == "1old"
            cv_l = cv_loc_old(ens.beta)
        elseif impr_set == "2"
            cv_l = cv_loc_set2(ens.beta)
        else
            error("Impr Set $impr_set not recoginsed, please choose between '1', '1old' and '2'.")
        end
        # [1:Int64(HVPobs.Data.get_T(ens.id)/2+1)]
        # [improve_corr_vkvk!(corrVec[n], corrVec_JPi[n], 2*cv_l, std=std, treelevel=true) for n in collect(1:nmax)]
        corrVec = [corrVec_PiPi[n] .+ (2*cv_l).*corrVec_Impr[n] for n=1:nmax]
    else
        corrVec = corrVec_PiPi
    end

    if RENORM
        Z3 = get_Z3(ens, impr_set=impr_set)
        [renormalize!(corr, Z3^2) for corr in corrVec]
        if !total
            [renormalize!(corr, Z3^2) for corr in corrVec_PiPi]
            [renormalize!(corr, Z3^2) for corr in corrVec_Impr]
        end
    end

    !total ? (return corrVec) : (return corrVec, corrVec_PiPi, corrVec_Impr)
end

reconstr_corr(ensid::String,E::Vector{uwreal},Z::Vector{uwreal},Zimpr::Vector{uwreal};nmax::Union{Nothing,Int64}=nothing,impr_set::String="1",IMPR::Bool=true,RENORM::Bool=true,FreeEN::Bool=false,total::Bool=false,) = reconstr_corr(EnsInfo(ensid),E::Vector{uwreal},Z::Vector{uwreal},Zimpr::Vector{uwreal};nmax=nmax,impr_set=impr_set,IMPR=IMPR,RENORM=RENORM,FreeEN=FreeEN,total=total,)

function findfirst_uninterrupted(bools)
    i = 0
    stop = false
    while !stop
        arg = findfirst(bools[i+1:end])
        i = !isnothing(arg) ? (i+arg) : nothing
        if isnothing(i) || isnothing(findfirst(.!bools[i+1:end]))  # stops if there are unterrupted trues or if it reached the end
            stop = true
        end
    end
    return i
end