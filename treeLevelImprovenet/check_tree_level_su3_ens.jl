######################################################################################
# This file was created by Alessandro Conigli 
# Here we load the finite-a HVP  stored in BDIO and perform
# the chiral-continuum extrapolations
########################################################################################

using Revise
using HVPobs
using PyPlot, LaTeXStrings
using OrderedCollections
using TimerOutputs
using BDIO, ADerrors
using QuadGK
using ALPHAio
import ADerrors: err

#================ SET UP VARIABLES ===============#

rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
rcParams["text.usetex"] =  true
rcParams["mathtext.fontset"]  = "cm"
rcParams["font.size"] =13
rcParams["axes.labelsize"] =22
rcParams["axes.titlesize"] = 18

plt.rc("text", usetex=true) # set to true if a LaTeX installation is present

#======== INCLUDE, PATHS AND CONSTANTS ==========#


include("../utils/types.jl")
include("../utils/tools.jl")
include("../utils/const.jl")
# include("../utils/data_management.jl")
include("../utils/plot_utils.jl")
include("../utils/IO_BDIO.jl")
include("../PI33/piQ_piQ4/func_comb_PI33.jl")

path_bdio = "/Users/alessandroconigli/MyDrive/postdoc-mainz/projects/deltalpha/data"
path_plot = "/Users/alessandroconigli/MyDrive/postdoc-mainz/projects/deltalpha/plots"
path_store_pi = "/Users/alessandroconigli/MyDrive/postdoc-mainz/projects/deltalpha/PIdata/impr_deriv/high_q_kernel/scale_error_multimom/"

#======= PHYSICAL CONSTANTS ====================#
const MPI_ph = uwreal([134.9768, 0.0005], "mpi phys")
const MK_ph  = uwreal([495.011, 0.01], "mk phys")

const phi2_ph = (sqrt(8)*t0sqrt_ph * MPI_ph / hc)^2
const phi4_ph = (sqrt(8)*t0sqrt_ph)^2 * ((MK_ph/hc)^2 + 0.5*(MPI_ph/hc)^2)


#============== READ BDIO FILES =================#
enslist = sort(["H101", "B450", "N202", "N300", "J500"])

dir_path = filter(x-> basename(x) in enslist, readdir(path_bdio, join=true))
spectrum_path = vcat(filter(!isempty, [filter(x-> occursin("spectrum.bdio", x)  , readdir(dir_path[k], join=true)) for k in eachindex(dir_path)])...)

ensinfo = EnsInfo.(enslist)
NENS = length(ensinfo)

a28t0 = Vector{uwreal}(undef, length(enslist))
for (k, ens) in enumerate(ensinfo)
    println("- Ensemble ", ens.id)
    t0ens = read_BDIO(spectrum_path[k], "spectrum", "t0")[1]
    # a28t0[k] =  1 / (8*t0ens)
    a28t0[k] =  a(ens.beta)^2
end

#============= READ PI 33  FROM BDIO =================#
# set 1 improvement coefficients
pi_33_ll_SD_s1 = Vector{Vector{uwreal}}(undef, 0)
pi_33_lc_SD_s1 = Vector{Vector{uwreal}}(undef, 0)

# set 2 improvement coefficients
pi_33_ll_SD_s2 = Vector{Vector{uwreal}}(undef, 0)
pi_33_lc_SD_s2 = Vector{Vector{uwreal}}(undef, 0)

fb = BDIO_open(joinpath(path_store_pi, "PI_33.bdio"), "r")
res = Dict()
count=0
while ALPHAdobs_next_p(fb)
    count+=1
    d = ALPHAdobs_read_parameters(fb)
    sz = tuple(d["size"]...)
    extra = d["extra"]
    ks = collect(d["keys"])
    res = ALPHAdobs_read_next(fb, size=sz, keys=ks)
    if extra["Ens"] ∉ enslist # != enslist[count]
        @info("Mismatch with EnsID in current BDIO uinfo2 ")
        println(extra["Ens"], " ", enslist)
        continue
    end
    push!(pi_33_ll_SD_s1, res["pi33_ll_SD_s1"])
    push!(pi_33_lc_SD_s1, res["pi33_lc_SD_s1"])
    push!(pi_33_ll_SD_s2, res["pi33_ll_SD_s2"])
    push!(pi_33_lc_SD_s2, res["pi33_lc_SD_s2"])

end
BDIO_close!(fb)
##

# computing the perturbative cl prediction
Qgev = [0.05, 0.1, 0.4, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0] # Q^2
Qmgev = 36.0 # Qm^2
pt_pred = []
for q in Qgev
    f(x) = treelevel_continuum_correlator.(x) .* krnl_dα_qhalf_sub.(x, q*1e6/hc^2, Qmgev*1e6/hc^2) .* Window("SD")(x)
    res, _ = quadgk(f,0,5,rtol=1e-3)
    println(res/2)
    push!(pt_pred, res/2)
end

##  cancelling fluctuations from t0_ph
NOERR = true
if NOERR
    for (k,ens) in enumerate(ensinfo)
        uwerr.(pi_33_ll_SD_s1[k])
        uwerr.(pi_33_lc_SD_s1[k])
        uwerr.(pi_33_ll_SD_s2[k])
        uwerr.(pi_33_lc_SD_s2[k])

        for (j,q) in enumerate(Qgev)
            set_fluc_to_zero!(pi_33_ll_SD_s1[k][j], "sqrtt0 [fm]")
            set_fluc_to_zero!(pi_33_lc_SD_s1[k][j], "sqrtt0 [fm]")
            set_fluc_to_zero!(pi_33_ll_SD_s2[k][j], "sqrtt0 [fm]")
            set_fluc_to_zero!(pi_33_lc_SD_s2[k][j], "sqrtt0 [fm]")

            pi_33_ll_SD_s1[k][j] *= 1.0
            pi_33_lc_SD_s1[k][j] *= 1.0
            pi_33_ll_SD_s2[k][j] *= 1.0
            pi_33_lc_SD_s2[k][j] *= 1.0
        end

        uwerr.(pi_33_ll_SD_s1[k])
        uwerr.(pi_33_lc_SD_s1[k])
        uwerr.(pi_33_ll_SD_s2[k])
        uwerr.(pi_33_lc_SD_s2[k])
    end
end
## reading tree level
fb = BDIO_open(joinpath(path_store_pi,  "treeLevelSet1.bdio"), "r")
res_3l_s1 = Dict()
while ALPHAdobs_next_p(fb)
    d = ALPHAdobs_read_parameters(fb)
    nobs = d["nobs"]
    dims = d["dimensions"]
    sz = tuple(d["size"]...)
    extra = d["extra"]
    if extra["Ens"] ∉ enslist # != enslist[count]
        @info("Mismatch with EnsID in current BDIO uinfo2 ")
        println(extra["Ens"], " ", enslist)
        continue
    end
    ks = collect(d["keys"])
    res_3l_s1[extra["Ens"]] = ALPHAdobs_read_next(fb, size=sz, keys=ks)
end
BDIO_close!(fb)

fb = BDIO_open(joinpath(path_store_pi,  "treeLevelSet2.bdio"), "r")
res_3l_s2 = Dict()
while ALPHAdobs_next_p(fb)
    d = ALPHAdobs_read_parameters(fb)
    nobs = d["nobs"]
    dims = d["dimensions"]
    sz = tuple(d["size"]...)
    extra = d["extra"]
    if extra["Ens"] ∉ enslist # != enslist[count]
        @info("Mismatch with EnsID in current BDIO uinfo2 ")
        println(extra["Ens"], " ", enslist)
        continue
    end
    ks = collect(d["keys"])
    res_3l_s2[extra["Ens"]] = ALPHAdobs_read_next(fb, size=sz, keys=ks)
end
BDIO_close!(fb)

hvp_3l_ll_s1 = Vector{Vector{uwreal}}(undef, 0)
hvp_3l_lc_s1 = Vector{Vector{uwreal}}(undef, 0)
hvp_3l_ll_s2 = Vector{Vector{uwreal}}(undef, 0)
hvp_3l_lc_s2 = Vector{Vector{uwreal}}(undef, 0)

## rearrange 3l impro data
for (k, ens) in enumerate(ensinfo)
    push!(hvp_3l_ll_s1, getindex(res_3l_s1, ens.id)["3l_33_ll"] ./ 2)
    push!(hvp_3l_lc_s1, getindex(res_3l_s1, ens.id)["3l_33_lc"] ./ 2)
    push!(hvp_3l_ll_s2, getindex(res_3l_s2, ens.id)["3l_33_ll"] ./ 2)
    push!(hvp_3l_lc_s2, getindex(res_3l_s2, ens.id)["3l_33_lc"] ./ 2)

end



## CHECK CL


uwerr.(a28t0)
@. cl(x,p) = p[1] + p[2]*x + p[3]*x^(3/2.) #+ p[4]*x^2
idx_beta = findall(x->x>3.34, getfield.(ensinfo, :beta))
# for (k, q) in enumerate(Qgev)
for k in [6,8,12]
    ydata_ll = getindex.(pi_33_ll_SD_s1, k); uwerr.(ydata_ll)
    ydata_lc = getindex.(pi_33_lc_SD_s1, k); uwerr.(ydata_lc)

    #fit_ll = fit_routine(cl, value.(a28t0)[idx_beta], ydata_ll[idx_beta], 3)
    fit_lc = fit_routine(cl, value.(a28t0)[idx_beta], ydata_lc[idx_beta], 3)
    uwerr.(fit_lc.param)    
    #ph_ll = fit_ll.param[1] 
    ph_lc = fit_lc.param[1]
    println(ph_lc)
    
    # errorbar(value.(a28t0), xerr=err.(a28t0), value.(ydata_ll), yerr=err.(ydata_ll), fmt="^", mfc="none", capsize=2, ms=8, label=L"$\mathrm{LL}$", color="orange" )
    # xarr = Float64.(range(0.0, maximum(value.(a28t0)[idx_beta]), length=100))
    # yarr = cl(xarr, fit_ll.param); uwerr.(yarr)
    # fill_between(xarr, value.(yarr)-err.(yarr), value.(yarr)+err.(yarr), alpha=0.2, color="orange")

    errorbar(value.(a28t0), xerr=err.(a28t0), value.(ydata_lc), yerr=err.(ydata_lc), fmt="^", mfc="none", capsize=2, ms=8, label=L"$\mathrm{LC}$", color="red" )
    xarr = Float64.(range(0.0, maximum(value.(a28t0[idx_beta])), length=100))
    yarr = cl(xarr, fit_lc.param); uwerr.(yarr)
    fill_between(xarr, value.(yarr)-err.(yarr), value.(yarr)+err.(yarr), alpha=0.2, color="red")
    
    # tree-level improved
    # additive
    # ydata_ll_3l = ydata_ll .-  value.(getindex.(hvp_3l_33_ll_s1, k)) .+ pt_pred[k]; uwerr.(ydata_ll_3l)
    # ydata_lc_3l = ydata_lc .-  value.(getindex.(hvp_3l_33_lc_s1, k)) .+ pt_pred[k]; uwerr.(ydata_lc_3l)
    # multiplcative
    ydata_ll_3l = ydata_ll .* pt_pred[k] ./ value.(getindex.(hvp_3l_ll_s1, k)) ; uwerr.(ydata_ll_3l)
    ydata_lc_3l = ydata_lc .* pt_pred[k] ./ value.(getindex.(hvp_3l_lc_s1, k)) ; uwerr.(ydata_lc_3l)
    
    #fit_ll_3l = fit_routine(cl, value.(a28t0)[idx_beta], ydata_ll_3l[idx_beta], 3)
    fit_lc_3l = fit_routine(cl, value.(a28t0)[idx_beta], ydata_lc_3l[idx_beta], 3)
    uwerr.(fit_lc_3l.param)
    uwerr.(fit_lc_3l.param)    
    #ph_ll_3l = fit_ll_3l.param[1] 
    ph_lc_3l = fit_lc_3l.param[1]
    println(ph_lc_3l)
    println()

    # errorbar(value.(a28t0), xerr=err.(a28t0), value.(ydata_ll_3l), yerr=err.(ydata_ll_3l), fmt="s", mfc="none", ms=8, capsize=2, label=L"$\mathrm{LL-tl}$", color="red" )
    # xarr = Float64.(range(0.0, maximum(value.(a28t0)[idx_beta]), length=100))
    # yarr = cl(xarr, fit_ll_3l.param); uwerr.(yarr)
    # fill_between(xarr, value.(yarr)-err.(yarr), value.(yarr)+err.(yarr), alpha=0.2, color="red")

    errorbar(value.(a28t0), xerr=err.(a28t0), value.(ydata_lc_3l), yerr=err.(ydata_lc_3l), fmt="s", mfc="none", ms=8, capsize=2, label=L"$\mathrm{LC  \ tree-level \ impr.}$", color="navy" )
    xarr = Float64.(range(0.0, maximum(value.(a28t0)[idx_beta]), length=100))
    yarr = cl(xarr, fit_lc_3l.param); uwerr.(yarr)
    fill_between(xarr, value.(yarr)-err.(yarr), value.(yarr)+err.(yarr), alpha=0.2, color="navy")


    #title(L"$\mathrm{Set \ 1}$")
    axvline(ls="-.", color="black", lw=0.4, alpha=0.7) 
    xlabel(L"$a^2 \ [\mathrm{fm}^2]$")
    ylabel(L"$(\Delta\alpha^{3,3})^{\mathrm{SD}}_{\mathrm{sub}}$")
    # legend(loc="lower left", ncol=2)
    legend(loc="upper left", ncol=1)
    tight_layout()
    display(gcf())
    savefig(joinpath(path_plot,"su3_3level_mult_q2_$(k)_set1.pdf"))
    close("all")
    
end

##

tight_layout()
display(gcf())
close("all")

