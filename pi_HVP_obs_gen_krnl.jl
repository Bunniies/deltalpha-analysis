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
const t0sqrt_ph = uwreal([0.1443, 0.0007], "sqrtt0 [fm]") 

include("./pi_HVP_types.jl")
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
const IMPR_SET  = "2" # either "1" or "2"
const RENORM    = true
const STD_DERIV = false
const KRNL = krnl_dα_qhalf_sub # kernel for isovector component
const KRNL88 = krnl_dα_qhalf   # kernel for isoscalar component
const WindSD = Window("SD")
const WindILD = Window("ILD")

enslist = ["N300", "J303", "E300"]
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

gdelta_iso_ll = Vector{Corr}(undef, length(ensinfo))
gdelta_iso_lc = Vector{Corr}(undef, length(ensinfo))

g08_ll_conn = Vector{Corr}(undef, length(ensinfo))
g08_lc_conn = Vector{Corr}(undef, length(ensinfo))

obs_SD = Vector(undef, length(ensinfo)) # short distance obs
obs_ILD = Vector(undef, length(ensinfo)) # long distance obs



#================ READING DATA ====================#

@info("Reading data...")
@time begin

    for (k, ens) in enumerate(ensinfo)
        
        println("\n - Ensemble: $(ens.id)")

        println("   - G33 ll and lc correlator")
        g33_ll[k], g33_lc[k] = corr33(path_data, ens, sector="light", path_rw=path_rw, impr=IMPR, impr_set=IMPR_SET, cons=true, frw_bcwd=true, std=STD_DERIV)
                
        println("   - G88  ll and lc correlator")
        
        g88_ll_conn[k], g88_lc_conn[k] = corr88_conn(path_data, ens, g33_ll[k], g33_lc=g33_lc[k], path_rw=path_rw, impr=IMPR, impr_set=IMPR_SET, cons=true, frw_bcwd=true, std=STD_DERIV)
        if ens.kappa_l == ens.kappa_s
            println("Symmetric point ensemble") 

        else
            #g88_ll[k], g88_lc[k] = corr88(path_data, ens, g33_ll[k], g33_lc=g33_lc[k], path_rw=path_rw, impr=IMPR, impr_set=IMPR_SET, cons=true, frw_bcwd=true, std=STD_DERIV)
        end

       
        println("   - G08 connected ll and lc correlator")
        g08_ll_conn[k], g08_lc_conn[k] = corr08_conn(g33_ll[k], g88_ll_conn[k], g33_lc=g33_lc[k], g88_lc=g88_lc_conn[k])

        if RENORM

            Z3 = get_Z3(ens, impr_set=IMPR_SET)
            renormalize!(g33_ll[k], Z3^2)
            renormalize!(g33_lc[k], Z3)
            Z8 = get_Z8(ens, impr_set=IMPR_SET)
            renormalize!(g88_ll_conn[k], Z8^2)
            renormalize!(g88_lc_conn[k], Z8)
            Z08 = get_Z08(ens, impr_set=IMPR_SET)
            renormalize!(g08_ll_conn[k], Z8*Z08)
            renormalize!(g08_lc_conn[k], Z08)
        end

        gdelta_iso_ll[k] = Corr(-g88_ll_conn[k].obs + g33_ll[k].obs, g33_ll[k].id, "deltaiso_ll")
        gdelta_iso_lc[k] = Corr(-g88_lc_conn[k].obs + g33_lc[k].obs, g33_lc[k].id, "deltaiso_lc")

        

        println("   - Gradient flow t0")
        obs_SD[k] = OrderedDict()
        obs_ILD[k] = OrderedDict()
        obs_SD[k]["t0"] = get_t0(path_ms, ens, path_rw=path_rw, pl=false)
        obs_ILD[k]["t0"] = obs_SD[k]["t0"]
    end
end
##

#=========== COMPUTE FVC  ==============#
##
for (k, ens) in enumerate(ensinfo)
    
    Qlat = Qgev .* value.(t0sqrt_ph.^2 ./ obs_SD[k]["t0"]) ./hc^2 *1e6
    # Qlat = 1.0 .* 0.1467.^2 ./ obs[k]["t0"] ./hc^2 *1e6

    fvc_raw_pi = get_fvc(joinpath(path_fvc, "JKMPI" ), ens.id)
    fvc_pi = vcat(sum(fvc_raw_pi, dims=1)...)
    
    fvc_raw_k = get_fvc(joinpath(path_fvc, "JKMK" ), ens.id)
    fvc_k = vcat(sum(fvc_raw_k, dims=1)...)
    
    T = Float64.(collect(1:length(fvc_pi)))

    aux_SD_pi = Vector{uwreal}(undef, 0)
    aux_SD_k = Vector{uwreal}(undef, 0)
    aux_ILD_pi = Vector{uwreal}(undef, 0)
    aux_ILD_k = Vector{uwreal}(undef, 0)

    for (j,q) in enumerate(Qlat) 
        qmlat = Qmgev * value.(t0sqrt_ph^2 / obs_SD[k]["t0"]) ./hc^2 *1e6
        krnl = KRNL.(T, q, qmlat)
        # pion FVC
        push!(aux_SD_pi, sum(fvc_pi .* krnl .* WindSD(T.*value(a(ens.beta))))) 
        push!(aux_ILD_pi, sum(fvc_pi .* krnl .* WindILD(T.*value(a(ens.beta)))))
        # kaon FVC
        push!(aux_SD_k, sum(fvc_k .* krnl .* WindSD(T.*value(a(ens.beta))))) 
        push!(aux_ILD_k, sum(fvc_k .* krnl .* WindILD(T.*value(a(ens.beta)))))
    end
    
    # pion FVC
    obs_SD[k]["fvc_pi"]  =  abs.(aux_SD_pi)
    obs_ILD[k]["fvc_pi"] =  abs.(aux_ILD_pi)
    # kaon FVC
    obs_SD[k]["fvc_k"]  =  abs.(aux_SD_k)
    obs_ILD[k]["fvc_k"] =  abs.(aux_ILD_k)
end

#=========== COMPUTE TMR ==============#
##
for (k, ens) in enumerate(ensinfo)
    println("- Ensemble: ", ens.id)

    Qlat = Qgev .* value.(t0sqrt_ph.^2) ./ obs_SD[k]["t0"] ./hc^2 *1e6
    qmlat = Qmgev * value(t0sqrt_ph^2) / obs_SD[k]["t0"] /hc^2 *1e6

    # Qlat = 1.0 .* 0.1467.^2 ./ obs[k]["t0"] ./hc^2 *1e6

    pi_33_ll, pi_33_lc, pi_deltaiso_ll, pi_deltaiso_lc = [Vector{uwreal}(undef, 0) for _ = 1:4]
    pi_33_ll_ILD, pi_33_lc_ILD, pi_deltaiso_ll_ILD, pi_deltaiso_lc_ILD = [Vector{uwreal}(undef, 0) for _ = 1:4]
    

    for (j, q) in enumerate(Qlat)    

        println("   - SD window")
        #isovector
        push!(pi_33_ll, tmr_integrand(g33_ll[k], q, qmlat, KRNL, pl=false, t0ens=obs_SD[k]["t0"], wind=WindSD))
        push!(pi_33_lc, tmr_integrand(g33_lc[k], q, qmlat, KRNL, pl=false, t0ens=obs_SD[k]["t0"], wind=WindSD))
        # isovector - isoscalar
        push!(pi_deltaiso_ll, tmr_integrand(gdelta_iso_ll[k], q, KRNL88, pl=true, t0ens=obs_SD[k]["t0"]))
        push!(pi_deltaiso_lc, tmr_integrand(gdelta_iso_lc[k], q, KRNL88, pl=false, t0ens=obs_SD[k]["t0"]))
        
        println("   - ILD window")
        #isovector
        push!(pi_33_ll_ILD, tmr_integrand(g33_ll[k], q, qmlat, KRNL, pl=false, t0ens=obs_SD[k]["t0"], wind=WindILD))
        push!(pi_33_lc_ILD, tmr_integrand(g33_lc[k], q, qmlat, KRNL, pl=false, t0ens=obs_SD[k]["t0"], wind=WindILD))
        # isovector-isoscalar
        # push!(pi_deltaiso_ll_ILD, tmr_integrand(gdelta_iso_ll[k], q, KRNL88, pl=false, t0ens=obs_SD[k]["t0"]))
        # push!(pi_deltaiso_lc_ILD, tmr_integrand(gdelta_iso_lc[k], q, KRNL88, pl=false, t0ens=obs_SD[k]["t0"]))

    end
    

    # SD window
    obs_SD[k]["pi_33_ll"] = pi_33_ll #.+ obs[k]["fvc"]
    obs_SD[k]["pi_33_lc"] = pi_33_lc #.+ obs[k]["fvc"]

    obs_SD[k]["pi_dltiso_ll"] = ens.kappa_l==ens.kappa_s ? 0.0 : pi_deltaiso_ll #.+ obs[k]["fvc"]
    obs_SD[k]["pi_dltiso_lc"] = ens.kappa_l==ens.kappa_s ? 0.0 : pi_deltaiso_lc #.+ obs[k]["fvc"]

    obs_SD[k]["pi_08_ll_conn"] = ens.kappa_l==ens.kappa_s ? 0.0 : sqrt(3)/2 * (pi_33_ll - pi_deltaiso_ll) # pi_08_ll_conn #.+ (all(iszero.(value.(pi_08_ll_conn))) ? 0.0 : obs[k]["fvc"] ) 
    obs_SD[k]["pi_08_lc_conn"] = ens.kappa_l==ens.kappa_s ? 0.0 : sqrt(3)/2 * (pi_33_lc - pi_deltaiso_lc) #pi_08_lc_conn #.+ (all(iszero.(value.(pi_08_lc_conn))) ? 0.0 : obs[k]["fvc"] )

    # LD window
    obs_ILD[k]["pi_33_ll"] = pi_33_ll_ILD
    obs_ILD[k]["pi_33_lc"] = pi_33_lc_ILD

    obs_ILD[k]["pi_dltiso_ll"] = obs_SD[k]["pi_dltiso_ll"]
    obs_ILD[k]["pi_dltiso_lc"] = obs_SD[k]["pi_dltiso_lc"]

    obs_ILD[k]["pi_08_ll_conn"] = ens.kappa_l==ens.kappa_s ? 0.0 : sqrt(3)/2 * (pi_33_ll_ILD - pi_deltaiso_ll)
    obs_ILD[k]["pi_08_lc_conn"] = ens.kappa_l==ens.kappa_s ? 0.0 : sqrt(3)/2 * (pi_33_lc_ILD - pi_deltaiso_lc)
end
##
#========== ADD FVC ==========#
for k in eachindex(ensinfo)
    if ensinfo[k].kappa_l == ensinfo[k].kappa_s
        println("- Symmetric point ensemble")
        # isovector component
        obs_SD[k]["pi_33_ll"]  .+= 1.5 * obs_SD[k]["fvc_pi"]
        obs_SD[k]["pi_33_lc"]  .+= 1.5 * obs_SD[k]["fvc_pi"]
        obs_ILD[k]["pi_33_ll"] .+= 1.5 * obs_ILD[k]["fvc_pi"]
        obs_ILD[k]["pi_33_lc"] .+= 1.5 * obs_ILD[k]["fvc_pi"]  

    else
        # isovector component
        obs_SD[k]["pi_33_ll"]  .+= obs_SD[k]["fvc_pi"] .+ obs_SD[k]["fvc_k"]
        obs_SD[k]["pi_33_lc"]  .+= obs_SD[k]["fvc_pi"] .+ obs_SD[k]["fvc_k"]
        obs_ILD[k]["pi_33_ll"] .+= obs_ILD[k]["fvc_pi"] .+ obs_ILD[k]["fvc_k"]
        obs_ILD[k]["pi_33_lc"] .+= obs_ILD[k]["fvc_pi"] .+ obs_ILD[k]["fvc_k"]

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
    Qlat = Qgev[3] * value(t0sqrt_ph^2) / obs_SD[k]["t0"] /hc^2 *1e6
    qmlat = Qmgev * value(t0sqrt_ph^2) / obs_SD[k]["t0"] /hc^2 *1e6

    pi33_ild, data_ild = tmr_integrand(g33_ll[k], Qlat, qmlat, KRNL, pl=false, t0ens=obs_SD[k]["t0"], wind=Window("ILD"), data=true)
    pi33_sd, data_sd = tmr_integrand(g33_ll[k], Qlat, qmlat, KRNL, pl=false, t0ens=obs_SD[k]["t0"], wind=Window("SD"), data=true)
    pi_88, data_88 = tmr_integrand(gdelta_iso_ll[k], Qlat, KRNL88, pl=false, t0ens=obs_SD[k]["t0"], data=true)

    data_88 .*=10
    uwerr.(data_ild)
    uwerr.(data_sd)
    uwerr.(data_88)

    xx = collect(1:length(data_ild)) .* value(a(ens.beta))
    errorbar(xx, value.(data_ild), err.(data_ild), fmt="v", color="red", capsize=2, label=L"$G^{(3,3)}_{\mathrm{ILD}}$")
    errorbar(xx, value.(data_sd), err.(data_sd), fmt="^", color="blue",  capsize=2, label=L"$G^{(3,3)}_{\mathrm{SD}}$")
    errorbar(xx, value.(data_88), err.(data_88), fmt="^", color="green",  capsize=2, label=L"$10\cdot(G^{(3,3)} - G^{(8,8)}_{\mathrm{conn}})$")
 
    xlabel(L"$t \ [\mathrm{fm}]$")
    ylabel(L"$G^{(d,e)}(t) \cdot K(t, Q^2)$")
    xlim(-0.05,1.5) 
    # ylim(-0.0001, 0.0001)
    legend()
    tight_layout()
    display(gcf())
    # savefig(joinpath(path_plot, "ensembles", ens.id, "$(ens.id)_dalpha_kernel_SD_vs_ILD.pdf"))
    close()
end


#======== SAVE TO BDIO =========#
## SD
[pop!(obs_SD[k],"fvc_k") for k in eachindex(obs_SD)]
for (k, o) in enumerate(obs_SD)

    ens = ensinfo[k].id
    p = joinpath(path_bdio, ens,  string(ens, "_PI_hvp_set$(IMPR_SET)_SD.bdio"))
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
# ILD
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
