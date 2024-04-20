######################################################################################
# This file was created by Alessandro Conigli - 03/2024
# Here we compute the finite-a HVP required in the analysis of the running of δα
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

#================ SET UP PLOT PARAMETERS ==================#

rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
rcParams["text.usetex"] =  true
rcParams["mathtext.fontset"]  = "cm"
rcParams["font.size"] =13
rcParams["axes.labelsize"] =22
rcParams["axes.titlesize"] = 18
plt.rc("text", usetex=true) # set to true if a LaTeX installation is present

#======== INCLUDE, PATHS AND CONSTANTS ==========#
# Madrid scale setting
# const t0sqrt_ph = uwreal([0.1443, 0.0007], "sqrtt0 [fm]") 

include("const.jl")
include("pi_HVP_types.jl")
include("tools.jl")
include("data_management.jl")
include("plot_utils.jl")

const path_data = "/Users/alessandroconigli/Lattice/data/HVP/2ptdata"
const path_rw   = "/Users/alessandroconigli/Lattice/data/HVP/rwf_deflated"
const path_ms   = "/Users/alessandroconigli/Lattice/data/HVP/ms_t0_dat"
const path_fvc  = "/Users/alessandroconigli/Lattice/data/HVP/FSE"

const path_plot = "/Users/alessandroconigli/MyDrive/postdoc-mainz/projects/deltalpha/plots"
const path_bdio = "/Users/alessandroconigli/MyDrive/postdoc-mainz/projects/deltalpha/data"

const IMPR      = true
const IMPR_SET  = "1" # either "1" or "2"
const RENORM    = true
const STD_DERIV = true
const KRNLsub = krnl_dα_qhalf_sub # subtracted kernel
const KRNL = krnl_dα_qhalf        # standard kernel 
const WindSD = Window("SD")
const WindILD = Window("ILD")

enslist = ["B450"]
ensinfo = EnsInfo.(enslist)

Qgev = [3., 5., 9.] # Q^2
Qmgev = 9.0 # Qm^2



#============ OBSERVABLE ALLOCATIONS ============#

g33_ll = Vector{Corr}(undef, length(ensinfo))
g33_lc = Vector{Corr}(undef, length(ensinfo))

g88_ll_conn = Vector{Corr}(undef, length(ensinfo))
g88_lc_conn = Vector{Corr}(undef, length(ensinfo))

g88_ll_disc = Vector{Corr}(undef, length(ensinfo))
g88_lc_disc = Vector{Corr}(undef, length(ensinfo))

gdelta_iso_ll = Vector{Corr}(undef, length(ensinfo)) # (G33 - G88) ll
gdelta_iso_lc = Vector{Corr}(undef, length(ensinfo)) # (G33 - G88) lc

g08_ll_conn = Vector{Corr}(undef, length(ensinfo))
g08_lc_conn = Vector{Corr}(undef, length(ensinfo))

g08_ll_disc = Vector{Corr}(undef, length(ensinfo))
g08_lc_disc = Vector{Corr}(undef, length(ensinfo))

g08_ll = Vector{Corr}(undef, length(ensinfo))
g08_lc = Vector{Corr}(undef, length(ensinfo))

gcc_ll_conn = Vector{Corr}(undef, length(ensinfo))
gcc_lc_conn = Vector{Corr}(undef, length(ensinfo))

gcc_ll_disc = Vector{Corr}(undef, length(ensinfo))
gcc_lc_disc = Vector{Corr}(undef, length(ensinfo))
gcc_cc_disc = Vector{Corr}(undef, length(ensinfo))

gc8_ll_disc = Vector{Corr}(undef, length(ensinfo))
gc8_lc_disc = Vector{Corr}(undef, length(ensinfo))
gc8_cc_disc = Vector{Corr}(undef, length(ensinfo))

obs = Vector(undef, length(ensinfo)) 



#================ READING DATA ====================#

@info("Reading data...")
@time begin

    for (k, ens) in enumerate(ensinfo)
        
        println("\n - Ensemble: $(ens.id)")

        println("        - G33 ll and lc correlator")
        g33_ll[k], g33_lc[k] = corr33(path_data, ens, sector="light", path_rw=path_rw, impr=IMPR, impr_set=IMPR_SET, cons=true, frw_bcwd=true, std=STD_DERIV)
                
        println("        - G88 connected ll and lc correlator")        
        g88_ll_conn[k], g88_lc_conn[k] = corr88_conn(path_data, ens, g33_ll[k], g33_lc=g33_lc[k], path_rw=path_rw, impr=IMPR, impr_set=IMPR_SET, cons=true, frw_bcwd=true, std=STD_DERIV)
       
        println("        - G08 connected ll and lc correlator")
        g08_ll_conn[k], g08_lc_conn[k] = corr08_conn(g33_ll[k], g88_ll_conn[k], g33_lc=g33_lc[k], g88_lc=g88_lc_conn[k])

        println("        - Gcc connected ll and lc correlator")
        try
            gcc_ll_conn[k], gcc_lc_conn[k] = corrcc_conn(path_data, ens, path_rw=path_rw, impr=IMPR, impr_set=IMPR_SET, cons=true, frw_bcwd=true, std=STD_DERIV)
        catch
            println("        - cc connected not found for ens $(ens.id)")
            T = HVPobs.Data.get_T(ens.id)
            gcc_ll_conn[k] = gcc_lc_conn[k] = Corr(fill(uwreal(0.0), T), ens.id, "Gcc") 
        end

        println("        - G88, Gcc, G08 disconnected")
        if ens.kappa_l == ens.kappa_s
            println("Symmetric point ensemble") 
            T = HVPobs.Data.get_T(ens.id)
            g88_ll_disc[k] = g88_lc_disc[k] = Corr(fill(uwreal(0.0), T), ens.id, "G88disc")  
            gcc_ll_disc[k] = gcc_lc_disc[k]  = gcc_cc_disc[k] = Corr(fill(uwreal(0.0), T), ens.id, "Gccdisc")
            g08_ll_disc[k] = g08_lc_disc[k] = Corr(fill(uwreal(0.0), T), ens.id, "G08disc")
            gc8_ll_disc[k] = gc8_lc_disc[k] = gc8_cc_disc[k] = Corr(fill(uwreal(0.0), T), ens.id, "Gc8disc")
                      
        else
            try    
                g88_ll_disc[k], g88_lc_disc[k] = corrDisconnected(path_data, ens, "88", path_rw=path_rw, impr=IMPR, impr_set=IMPR_SET,  frw_bcwd=true, cc=false, std=STD_DERIV)
                gcc_ll_disc[k], gcc_lc_disc[k], gcc_cc_disc[k] = corrDisconnected(path_data, ens, "cc", path_rw=path_rw, impr=IMPR, impr_set=IMPR_SET, cc=true, frw_bcwd=true, std=STD_DERIV)
                g08_ll_disc[k], g08_lc_disc[k] = corrDisconnected(path_data, ens, "08", path_rw=path_rw, impr=IMPR, impr_set=IMPR_SET, frw_bcwd=true, cc=false, std=STD_DERIV)
                gc8_ll_disc[k], gc8_lc_disc[k], gc8_cc_disc[k] = corrDisconnected(path_data, ens, "c8", path_rw=path_rw, impr=IMPR, impr_set=IMPR_SET, cc=true, frw_bcwd=true, std=STD_DERIV)
             catch
                println("disconnected failed")
                T = HVPobs.Data.get_T(ens.id)
                g88_ll_disc[k] = g88_lc_disc[k] = gcc_ll_disc[k] =  gcc_lc_disc[k] = g08_ll_disc[k] = g08_lc_disc[k] = gc8_ll_disc[k] = gc8_lc_disc[k] = gcc_cc_disc[k] =  gc8_cc_disc[k] = Corr(fill(uwreal(0.0), T), ens.id, "Gdisc")
              end
        end
        

        if RENORM

            Z3 = get_Z3(ens, impr_set=IMPR_SET)
            renormalize!(g33_ll[k], Z3^2)
            renormalize!(g33_lc[k], Z3)
            
            Z8 = get_Z8(ens, impr_set=IMPR_SET)
            renormalize!(g88_ll_conn[k], Z8^2)
            renormalize!(g88_lc_conn[k], Z8)
            renormalize!(g88_ll_disc[k], Z8^2)
            renormalize!(g88_lc_disc[k], Z8)

            Zcc = Zvc_l[ens.id]
            renormalize!(gcc_ll_conn[k], Zcc*Zcc)
            renormalize!(gcc_lc_conn[k], Zcc)
            renormalize!(gcc_ll_disc[k], Zcc*Zcc)
            renormalize!(gcc_lc_disc[k], Zcc)

            Z08 = get_Z08(ens, impr_set=IMPR_SET)
            renormalize!(g08_ll_conn[k], Z8*Z08)
            renormalize!(g08_lc_conn[k], Z08)
            renormalize!(g08_ll_disc[k], Z8*Z08)
            renormalize!(g08_lc_disc[k], Z08)
        end

        gdelta_iso_ll[k] = Corr(-g88_ll_conn[k].obs - g88_ll_disc[k].obs + g33_ll[k].obs, g33_ll[k].id, "deltaiso_ll")
        gdelta_iso_lc[k] = Corr(-g88_lc_conn[k].obs - g88_lc_disc[k].obs + g33_lc[k].obs, g33_lc[k].id, "deltaiso_lc")

        g08_ll[k] = Corr(g08_ll_conn[k].obs + g08_ll_disc[k].obs, g08_ll_conn[k].id, "G08ll")
        g08_lc[k] = Corr(g08_lc_conn[k].obs + g08_lc_disc[k].obs, g08_lc_conn[k].id, "G08lc")

        println("        - Gradient flow t0")
        obs[k] = OrderedDict()
        obs[k]["t0"] = get_t0(path_ms, ens, path_rw=path_rw, pl=false)
    end
end
##
yy =  gcc_cc_disc[1].obs; uwerr.(yy)
errorbar(collect(1:length(yy)), value.(yy), err.(yy), fmt="s")
ylim(-0.00001, 0.00001)
display(gcf())

close("all")
## test disco
corr_dict = get_corr_disc(path_data, ensinfo[1], "cc", frw_bcwd=false, path_rw=path_rw, L=1)

#=========== COMPUTE FVC  ==============#
##
for (k, ens) in enumerate(ensinfo)
    
    Qlat = Qgev .* value.(t0sqrt_ph.^2 ./ obs[k]["t0"]) ./hc^2 *1e6
    # Qlat = 1.0 .* 0.1467.^2 ./ obs[k]["t0"] ./hc^2 *1e6

    fvc_raw_pi = get_fvc(joinpath(path_fvc, "JKMPI" ), ens.id)
    fvc_pi = vcat(sum(fvc_raw_pi, dims=1)...)
    
    fvc_raw_k = get_fvc(joinpath(path_fvc, "JKMK" ), ens.id)
    fvc_k = vcat(sum(fvc_raw_k, dims=1)...)
    
    T = Float64.(collect(1:length(fvc_pi)))

    aux_SD_pi = Vector{uwreal}(undef, 0)   # sub kernel 
    aux_SD_k = Vector{uwreal}(undef, 0)    # sub kernel
    aux_ILD_pi = Vector{uwreal}(undef, 0)  # sub kernel
    aux_ILD_k = Vector{uwreal}(undef, 0)   # sub kernel
    aux_pi = Vector{uwreal}(undef, 0)      # standard kernel
    aux_k = Vector{uwreal}(undef, 0)       # standard kernel

    for (j,q) in enumerate(Qlat) 
        qmlat = Qmgev * value.(t0sqrt_ph^2 / obs[k]["t0"]) ./hc^2 *1e6
        krnlsub = KRNLsub.(T, q, qmlat)
        krnl = KRNL.(T, q)
        # pion FVC
        push!(aux_SD_pi, sum(fvc_pi .* krnlsub .* WindSD(T.*value(a(ens.beta))))) 
        push!(aux_ILD_pi, sum(fvc_pi .* krnlsub .* WindILD(T.*value(a(ens.beta)))))
        push!(aux_pi, sum(fvc_pi .* krnl))
        # kaon FVC
        push!(aux_SD_k, sum(fvc_k .* krnlsub .* WindSD(T.*value(a(ens.beta))))) 
        push!(aux_ILD_k, sum(fvc_k .* krnlsub .* WindILD(T.*value(a(ens.beta)))))
        push!(aux_k, sum(fvc_k .* krnl))

    end
    
    # pion FVC
    obs[k]["fvc_pi_SD"]  =  abs.(aux_SD_pi)
    obs[k]["fvc_pi_ILD"] =  abs.(aux_ILD_pi)
    obs[k]["fvc_pi"]     = abs.(aux_pi)
    # kaon FVC
    obs[k]["fvc_k_SD"]  =  abs.(aux_SD_k)
    obs[k]["fvc_k_ILD"] =  abs.(aux_ILD_k)
    obs[k]["fvc_k"]     = abs.(aux_k)

end

#=========== COMPUTE TMR ==============#
##
for (k, ens) in enumerate(ensinfo)
    println("- Ensemble: ", ens.id)

    Qlat = Qgev .* value.(t0sqrt_ph.^2) ./ obs[k]["t0"] ./hc^2 *1e6
    qmlat = Qmgev * value(t0sqrt_ph^2) / obs[k]["t0"] /hc^2 *1e6

    # Qlat = 1.0 .* 0.1467.^2 ./ obs[k]["t0"] ./hc^2 *1e6

    pi_33_ll_SD, pi_33_lc_SD, pi_33_ll_ILD, pi_33_lc_ILD = [Vector{uwreal}(undef, 0) for _ = 1:4]       # G33 ll lc SD ILD
    pi_deltaiso_ll, pi_deltaiso_lc = [Vector{uwreal}(undef, 0) for _ = 1:2]                             # (G33 - G88)
    pi_cc_ll_conn, pi_cc_lc_conn, pi_cc_ll_disc, pi_cc_lc_disc, pi_cc_cc_disc = [Vector{uwreal}(undef, 0) for _ = 1:5] # Gcc ll lc cc conn disc
    pi_08_ll, pi_08_lc = [Vector{uwreal}(undef, 0) for _ = 1:2]                                         # G08 ll lc
    pi_c8_ll_disc, pi_c8_lc_disc, pi_c8_cc_disc = [Vector{uwreal}(undef, 0) for _ = 1:3]                # Gc8 ll lc cc  disc

    

    for (j, q) in enumerate(Qlat)    
        println("   - G33 SD window")
        push!(pi_33_ll_SD, tmr_integrand(g33_ll[k], q, qmlat, KRNLsub, pl=false, t0ens=obs[k]["t0"], wind=WindSD))
        push!(pi_33_lc_SD, tmr_integrand(g33_lc[k], q, qmlat, KRNLsub, pl=false, t0ens=obs[k]["t0"], wind=WindSD))

        println("   - G33 ILD window")
        push!(pi_33_ll_ILD, tmr_integrand(g33_ll[k], q, qmlat, KRNLsub, pl=false, t0ens=obs[k]["t0"], wind=WindILD))
        push!(pi_33_lc_ILD, tmr_integrand(g33_lc[k], q, qmlat, KRNLsub, pl=false, t0ens=obs[k]["t0"], wind=WindILD))

        println("   - (G33 - G88) ")
        push!(pi_deltaiso_ll, tmr_integrand(gdelta_iso_ll[k], q, KRNL, pl=false, t0ens=obs[k]["t0"]))
        push!(pi_deltaiso_lc, tmr_integrand(gdelta_iso_lc[k], q, KRNL, pl=false, t0ens=obs[k]["t0"]))
        
        println("   - Gcc connected ")
        push!(pi_cc_ll_conn, tmr_integrand(gcc_ll_conn[k], q, qmlat, KRNLsub, pl=false, t0ens=obs[k]["t0"]))
        push!(pi_cc_lc_conn, tmr_integrand(gcc_lc_conn[k], q, qmlat, KRNLsub, pl=false, t0ens=obs[k]["t0"]))

        println("   - Gcc disconnected ")
        push!(pi_cc_ll_disc, tmr_integrand(gcc_ll_disc[k], q, KRNL, pl=false, t0ens=obs[k]["t0"]))
        push!(pi_cc_lc_disc, tmr_integrand(gcc_lc_disc[k], q, KRNL, pl=false, t0ens=obs[k]["t0"]))
        push!(pi_cc_cc_disc, tmr_integrand(gcc_cc_disc[k], q, KRNL, pl=true, t0ens=obs[k]["t0"]))

        println("   - Gc8 disconnected ")
        push!(pi_c8_ll_disc, tmr_integrand(gc8_ll_disc[k], q, KRNL, pl=false, t0ens=obs[k]["t0"]))
        push!(pi_c8_lc_disc, tmr_integrand(gc8_lc_disc[k], q, KRNL, pl=false, t0ens=obs[k]["t0"]))
        push!(pi_c8_cc_disc, tmr_integrand(gc8_cc_disc[k], q, KRNL, pl=false, t0ens=obs[k]["t0"]))

        println("   - G08")
        push!(pi_08_ll, tmr_integrand(g08_ll[k], q, qmlat, KRNLsub, pl=false, t0ens=obs[k]["t0"]))
        push!(pi_08_lc, tmr_integrand(g08_lc[k], q, qmlat, KRNLsub, pl=false, t0ens=obs[k]["t0"]))

        println("")

    end
    
    # G33 SD and ILD
    obs[k]["pi_33_ll_SD"] = pi_33_ll_SD 
    obs[k]["pi_33_lc_SD"] = pi_33_lc_SD 
    obs[k]["pi_33_ll_ILD"] = pi_33_ll_ILD
    obs[k]["pi_33_lc_ILD"] = pi_33_lc_ILD

    # (G33 - G88)
    obs[k]["pi_dltiso_ll"] = ens.kappa_l==ens.kappa_s ? 0.0 : pi_deltaiso_ll 
    obs[k]["pi_dltiso_lc"] = ens.kappa_l==ens.kappa_s ? 0.0 : pi_deltaiso_lc 

    # Gcc
    obs[k]["pi_cc_ll_conn"] = pi_cc_ll_conn
    obs[k]["pi_cc_lc_conn"] = pi_cc_lc_conn
    obs[k]["pi_cc_ll_disc"] = pi_cc_ll_disc
    obs[k]["pi_cc_lc_disc"] = pi_cc_lc_disc
    obs[k]["pi_cc_cc_disc"] = pi_cc_cc_disc

    # Gc8
    obs[k]["pi_c8_ll_disc"] = pi_c8_ll_disc
    obs[k]["pi_c8_lc_disc"] = pi_c8_lc_disc
    obs[k]["pi_c8_cc_disc"] = pi_c8_cc_disc

    # G08
    obs[k]["pi_08_ll"] = pi_08_ll #ens.kappa_l==ens.kappa_s ? 0.0 : sqrt(3)/2 * (pi_33_ll_SD - pi_deltaiso_ll) 
    obs[k]["pi_08_lc"] = pi_08_lc #ens.kappa_l==ens.kappa_s ? 0.0 : sqrt(3)/2 * (pi_33_lc_SD - pi_deltaiso_lc) 
end
##
#========== ADD FVC ==========#
for k in eachindex(ensinfo)
    if ensinfo[k].kappa_l == ensinfo[k].kappa_s
        println("- Symmetric point ensemble")
        # isovector component
        obs[k]["pi_33_ll_SD"]  .+= 1.5 * obs[k]["fvc_pi_SD"]
        obs[k]["pi_33_lc_SD"]  .+= 1.5 * obs[k]["fvc_pi_SD"]
        obs[k]["pi_33_ll_ILD"] .+= 1.5 * obs[k]["fvc_pi_ILD"]
        obs[k]["pi_33_lc_ILD"] .+= 1.5 * obs[k]["fvc_pi_ILD"]  

    else
        # isovector component
        obs[k]["pi_33_ll_SD"]  .+= obs[k]["fvc_pi_SD"] .+ obs[k]["fvc_k_SD"]
        obs[k]["pi_33_lc_SD"]  .+= obs[k]["fvc_pi_SD"] .+ obs[k]["fvc_k_SD"]
        obs[k]["pi_33_ll_ILD"] .+= obs[k]["fvc_pi_ILD"] .+ obs[k]["fvc_k_ILD"]
        obs[k]["pi_33_lc_ILD"] .+= obs[k]["fvc_pi_ILD"] .+ obs[k]["fvc_k_ILD"]
        # G33 - G88 deltaiso
        pi_fvc = obs[k]["fvc_pi"] 
        k_fvc =  obs[k]["fvc_k"]
        obs[k]["pi_dltiso_ll"] .+= (pi_fvc + k_fvc - 2/9*k_fvc) # the -2/9 terms comes from the strange-conn piece in G88
        obs[k]["pi_dltiso_lc"] .+= (pi_fvc + k_fvc - 2/9*k_fvc)
    end
end
##
# comparison with older BDIO
path_bdio_old = "/Users/alessandroconigli/MyDrive/postdoc-mainz/projects/deltalpha/data/old_kernel/sub_kernel/"
bb = read_BDIO(joinpath(path_bdio_old,"H102", "H102_PI_hvp_set1.bdio"), "dalpha", "pi_33_ll"); uwerr.(bb); bb
cc = read_BDIO(joinpath(path_bdio_old,"H102", "H102_PI_hvp_set1.bdio"), "dalpha", "pi_33_ll"); uwerr.(cc); cc

old_res = bb #.+ cc ; uwerr.(old_res); old_res
aa = obs_SD[2]["pi_33_ll"] + obs_ILD[2]["pi_33_ll"] ; uwerr.(aa); aa


#============ PLOT TMR =============#
##
for (k,ens) in enumerate(ensinfo)
    Qlat = Qgev[2] * value(t0sqrt_ph^2) / obs[k]["t0"] /hc^2 *1e6
    qmlat = Qmgev * value(t0sqrt_ph^2) / obs[k]["t0"] /hc^2 *1e6

    pi33, data_33 = tmr_integrand(g33_lc[k], Qlat, KRNL, pl=false, t0ens=obs[k]["t0"], data=true)
    pi_88, data_88 = tmr_integrand(gdelta_iso_lc[k], Qlat, KRNL, pl=false, t0ens=obs[k]["t0"], data=true)
    pi_cc_conn, data_cc_conn = tmr_integrand(gcc_lc_conn[k], Qlat, KRNL, pl=false, t0ens=obs[k]["t0"], data=true)
    pi_cc_disc, data_cc_disc = tmr_integrand(gcc_cc_disc[k], Qlat, KRNL, pl=false, t0ens=obs[k]["t0"], data=true)
    pi_c8_disc, data_c8_disc = tmr_integrand(gc8_cc_disc[k], Qlat, KRNL, pl=false, t0ens=obs[k]["t0"], data=true)
    pi_08, data_08 = tmr_integrand(g08_lc[k], Qlat, KRNL, pl=false, t0ens=obs[k]["t0"], data=true)

    data_88 .*=  10
    data_cc_conn .*= (4/9)
    data_cc_disc .*= (4/9) * 1000
    data_c8_disc .*= (2/(3*sqrt(3))) * 1000
    data_08 .*= 1000

    uwerr.(data_33)
    uwerr.(data_88)
    uwerr.(data_cc_conn)
    uwerr.(data_cc_disc)
    uwerr.(data_c8_disc)
    uwerr.(data_08)

    xx = collect(0:length(data_33)-1) .* value(a(ens.beta))
    errorbar(xx, value.(data_33), err.(data_33), fmt="v", color="red", capsize=2, label=L"$G^{(3,3)}$")
    errorbar(xx, value.(data_88), err.(data_88), fmt="^", color="green",  capsize=2, label=L"$10\cdot(G^{(3,3)} - G^{(8,8)})$")
    errorbar(xx, value.(data_cc_conn), err.(data_cc_conn), fmt="^", color="blue",  capsize=2, label=L"$\frac{4}{9}G^{(c,c)}_{\mathrm{conn}}$")
    errorbar(xx, value.(data_cc_disc), err.(data_cc_disc), fmt="^", color="plum",  capsize=2, label=L"$1000\cdot\frac{4}{9}G^{(c,c)}_{\mathrm{disc}}$")
    # errorbar(xx, value.(data_c8_disc), err.(data_c8_disc), fmt="^", color="royalblue",  capsize=2, label=L"$1000\cdot\frac{2}{3\sqrt{3}}G^{(c,8)}_{\mathrm{disc}}$")
    # errorbar(xx, value.(data_08), err.(data_08), fmt="^", color="gold",  capsize=2, label=L"$100\cdot G^{(0,8)}$")
 
    xlabel(L"$t \ [\mathrm{fm}]$")
    ylabel(L"$G^{(d,e)}(t) \cdot K(t, Q^2)$")
    xlim(-0.05,2.0) 
    # ylim(-0.0001, 0.0001)
    legend(ncol=1)
    tight_layout()
    display(gcf())
    # savefig(joinpath(path_plot, "ensembles", ens.id, "$(ens.id)_dalpha_kernel_SD_vs_ILD.pdf"))
    close()
end


#======== SAVE TO BDIO =========#
## SD
[pop!(obs[k],"fvc_pi_SD") for k in eachindex(obs)]
[pop!(obs[k],"fvc_pi_ILD") for k in eachindex(obs)]
[pop!(obs[k],"fvc_pi") for k in eachindex(obs)]
[pop!(obs[k],"fvc_k_SD") for k in eachindex(obs)]
[pop!(obs[k],"fvc_k_ILD") for k in eachindex(obs)]
[pop!(obs[k],"fvc_k") for k in eachindex(obs)]
[pop!(obs[k],"pi_cc_ll_disc") for k in eachindex(obs)]
[pop!(obs[k],"pi_cc_lc_disc") for k in eachindex(obs)]
[pop!(obs[k],"pi_c8_ll_disc") for k in eachindex(obs)]
[pop!(obs[k],"pi_c8_lc_disc") for k in eachindex(obs)]


for (k, o) in enumerate(obs)

    ens = ensinfo[k].id
    p = joinpath(path_bdio, ens,  string(ens, "_PI_hvp_set$(IMPR_SET).bdio"))
    fb = BDIO_open(p, "w", ens)
    uinfo = 0 

    for (key, value) in o
        if isa(value, uwreal)
            write_uwreal(value, fb, uinfo)
        elseif isa(value, Vector{uwreal})
            [write_uwreal(v, fb, uinfo) for v in value]
        end
        uinfo +=1
    end

    BDIO_close!(fb)
end

## ILD
[pop!(obs_ILD[k],"fvc_k") for k in eachindex(obs_SD)]
for (k, o) in enumerate(obs_ILD)

    ens = ensinfo[k].id
    p = joinpath(path_bdio, ens,  string(ens, "_PI_hvp_set$(IMPR_SET)_ILD.bdio"))
    fb = BDIO_open(p, "w", ens)
    uinfo = 0 

    for (key, value) in o
        if isa(value, uwreal)
            write_uwreal(value, fb, uinfo)
        elseif isa(value, Vector{uwreal})
            [write_uwreal(v, fb, uinfo) for v in value]
        end
        uinfo +=1
    end

    BDIO_close!(fb)
end



## test FVC and plotting at different accumulating wrapping

fvc = Vector{Matrix{uwreal}}(undef, length(ensinfo))
fvc_acc = Vector{Vector{Vector{uwreal}}}(undef, length(ensinfo))
for (k, ens) in enumerate(ensinfo)
    fvc[k] = get_fvc(path_fvc, ens.id)
    fvc_acc[k] = Vector{Vector{uwreal}}(undef, 0)

    qmlat = Qmgev * value(t0sqrt_ph^2 / obs_SD[k]["t0"]) /hc^2 *1e6

    qref = 3.0 * value(t0sqrt_ph^2 / obs_SD[k]["t0"]) /hc^2 *1e6
    T = 0
    for (idx,n) in enumerate([1, 2, 6])
        push!(fvc_acc[k],  vcat(sum(fvc[k][1:n,:], dims=1)...))
        T = Float64.(collect(0:length(fvc_acc[k][idx])-1))
        tmr = -1. *fvc_acc[k][idx] .*  KRNL.(T, qref, qmlat)  
        tmp = sum(tmr); 
        # [uwerr(tmp[ll], wpmm) for ll in eachindex(tmp)]
        uwerr(tmp)#, wpmm) 
        tmr .=  tmr ./ value.(t0sqrt_ph) .* sqrt.(obs_SD[k]["t0"]) # convert to phys units
        uwerr.(tmr)
        # plot(T, value.(tmr),  label=L"$\vec{n}^2 = $"* "$(n)")
        fill_between(T, value.(tmr)-err.(tmr), value.(tmr)+err.(tmr), alpha=0.8, lw=50, label=L"$F_\pi(-Q^2),\ \vec{n}^2 \leq$"* " $(n)" )
        println(sum(tmr))
    end
    yy = tmr_integrand(g33_ll[k], qref, qmlat, KRNL, t0ens=obs_SD[k]["t0"], pl=false, data=true)[2] ./ value.(t0sqrt_ph) .* sqrt.(obs_SD[k]["t0"])
    # [uwerr(yy[ll], wpmm) for ll in eachindex(yy)]
    uwerr.(yy)
    fill_between(T, 0, err.(yy), alpha=0.2, color="gray")
    legend()
    xlabel(L"$t/a$")
    ylabel(L"$\delta G^{33}(t)K(t,Q^2) \ [\mathrm{fm}^{-1}]$")
    tight_layout()
    display(gcf())
    t = "$(ens.id)_fvc_Q2_3.pdf"
    # savefig(joinpath(path_plot,t))
    close("all")

end
