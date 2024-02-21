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

#================ SET UP VARIABLES ===============#

rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
rcParams["text.usetex"] =  true
rcParams["mathtext.fontset"]  = "cm"
rcParams["font.size"] =13
rcParams["axes.labelsize"] =22
rcParams["axes.titlesize"] = 18

plt.rc("text", usetex=true) # set to true if a LaTeX installation is present

#======== INCLUDE, PATHS AND CONSTANTS ==========#


include("./pi_HVP_types.jl")
include("./pi_HVP_func_comb.jl")
include("./tools.jl")
include("./data_management.jl")
include("./plot_utils.jl")
include("./IO_BDIO.jl")

path_bdio = "/Users/alessandroconigli/MyDrive/postdoc-mainz/projects/deltalpha/data"
path_plot = "/Users/alessandroconigli/MyDrive/postdoc-mainz/projects/deltalpha/plots"

#======= PHYSICAL CONSTANTS ====================#
const t0sqrt_ph = uwreal([0.1439, 0.0006], "sqrtt0 [fm]") 
const MPI_ph = uwreal([134.9768, 0.0005], "mpi phys")
const MK_ph  = uwreal([495.011, 0.01], "mk phys")

const phi2_ph = (sqrt(8)*t0sqrt_ph * MPI_ph / hc)^2
const phi4_ph = (sqrt(8)*t0sqrt_ph)^2 * ((MK_ph/hc)^2 + 0.5*(MPI_ph/hc)^2)

const muon = uwreal([105.6583755, 0.0000023 ],"muon") * 1e-3
muon_m = uwreal([105.6583755, 0.0000023 ],"muon") * 1e-3
b_sd = 7.422e-10
alpha = 1/137
b_sd / muon^2 / alpha^2 / 4 * 9 * pi^2

(4) / (9*pi^2*25)*log(2)

#============== READ BDIO FILES =================#

dir_path = filter(isdir, readdir(path_bdio, join=true))
hvp_path_s1 = vcat(filter(!isempty, [filter(x-> occursin("hvp_set1.bdio", x)  , readdir(dir_path[k], join=true)) for k in eachindex(dir_path)])...)
hvp_path_s2 = vcat(filter(!isempty, [filter(x-> occursin("hvp_set2.bdio", x)  , readdir(dir_path[k], join=true)) for k in eachindex(dir_path)])...)

spectrum_path = vcat(filter(!isempty, [filter(x-> occursin("spectrum.bdio", x)  , readdir(dir_path[k], join=true)) for k in eachindex(dir_path)])...)

treelevel_path_s1 = vcat(filter(!isempty, [filter(x-> occursin("tree_level_set1.bdio", x)  , readdir(dir_path[k], join=true)) for k in eachindex(dir_path)])...)
treelevel_path_s2 = vcat(filter(!isempty, [filter(x-> occursin("tree_level_set2.bdio", x)  , readdir(dir_path[k], join=true)) for k in eachindex(dir_path)])...)


enslist = [basename(hvp_path_s1[k])[1:4] for k in eachindex(hvp_path_s1)]
idx_su3 = findall(x-> x in ["H101", "B450", "N202", "N300", "J500"], enslist)
enslist = enslist[idx_su3]

hvp_path_s1 = hvp_path_s1[idx_su3]
hvp_path_s2 = hvp_path_s2[idx_su3]
spectrum_path = spectrum_path[idx_su3]
treelevel_path_s1 = treelevel_path_s1[idx_su3]
treelevel_path_s2 = treelevel_path_s2[idx_su3]

ensinfo = EnsInfo.(enslist)
NENS = length(ensinfo)

# set 1 improvement coefficients
pi_33_ll_s1 = Vector{Vector{uwreal}}(undef, NENS)
pi_33_lc_s1 = Vector{Vector{uwreal}}(undef, NENS)

hvp_3l_33_ll_s1 = Vector{Vector{uwreal}}(undef, NENS)
hvp_3l_33_lc_s1 = Vector{Vector{uwreal}}(undef, NENS)

# set 2 improvement coefficients
pi_33_ll_s2 = Vector{Vector{uwreal}}(undef, NENS)
pi_33_lc_s2 = Vector{Vector{uwreal}}(undef, NENS)

hvp_3l_33_ll_s2 = Vector{Vector{uwreal}}(undef, NENS)
hvp_3l_33_lc_s2 = Vector{Vector{uwreal}}(undef, NENS)

a28t0 = Vector{uwreal}(undef, NENS)

for (k, ens) in enumerate(ensinfo)
    println("- Ensemble ", ens.id)

    # set 1 improvement coefficients
    pi_33_ll_s1[k] = read_BDIO(hvp_path_s1[k], "dalpha", "pi_33_ll" )
    pi_33_lc_s1[k] = read_BDIO(hvp_path_s1[k], "dalpha", "pi_33_lc" )

    hvp_3l_33_ll_s1[k] = read_BDIO(treelevel_path_s1[k], "3level", "3l_33_ll")
    hvp_3l_33_lc_s1[k] = read_BDIO(treelevel_path_s1[k], "3level", "3l_33_lc")

    # set 2 improvement coefficients
    pi_33_ll_s2[k] = read_BDIO(hvp_path_s2[k], "dalpha", "pi_33_ll" )
    pi_33_lc_s2[k] = read_BDIO(hvp_path_s2[k], "dalpha", "pi_33_lc" )

    hvp_3l_33_ll_s2[k] = read_BDIO(treelevel_path_s2[k], "3level", "3l_33_ll")
    hvp_3l_33_lc_s2[k] = read_BDIO(treelevel_path_s2[k], "3level", "3l_33_lc")

    t0ens = read_BDIO(spectrum_path[k], "spectrum", "t0")[1]
    a28t0[k] =  1 / (8*t0ens)
end

# reading tree-level prediction in the continuum
pt_pred = []
p = joinpath(path_bdio, string("tree_level_cl_phys_val.bdio"))
fb = BDIO_open(p, "r")
BDIO_seek!(fb)
if BDIO_get_uinfo(fb) == 0
    push!(pt_pred, read_uwreal(fb))
end
while BDIO_seek!(fb, 2)
    if BDIO_get_uinfo(fb) == 0 
        push!(pt_pred, read_uwreal(fb))
    end
end
BDIO_close!(fb)

## CHECK CL
uwerr.(a28t0)
Qgev = [3., 5., 9.]

@. cl(x,p) = p[1] + p[2]*x + p[3]*x^(3/2.) #+ p[4]*x^2
idx_beta = findall(x->x>3.4, getfield.(ensinfo, :beta))
for (k, q) in enumerate(Qgev)
    ydata_ll = getindex.(pi_33_ll_s2, k); uwerr.(ydata_ll)
    ydata_lc = getindex.(pi_33_lc_s2, k); uwerr.(ydata_lc)

    fit_ll = fit_routine(cl, value.(a28t0)[idx_beta], ydata_ll[idx_beta], 3)
    fit_lc = fit_routine(cl, value.(a28t0)[idx_beta], ydata_lc[idx_beta], 3)
    
    ph_ll = fit_ll.param[1] 
    ph_lc = fit_lc.param[1]
    
    errorbar(value.(a28t0), xerr=err.(a28t0), value.(ydata_ll), yerr=err.(ydata_ll), fmt="s", mfc="none", capsize=2, label=L"$\mathrm{LL}$", color="orange" )
    xarr = Float64.(range(0.0, maximum(value.(a28t0)[idx_beta]), length=100))
    yarr = cl(xarr, fit_ll.param); uwerr.(yarr)
    fill_between(xarr, value.(yarr)-err.(yarr), value.(yarr)+err.(yarr), alpha=0.2, color="orange")

    errorbar(value.(a28t0), xerr=err.(a28t0), value.(ydata_lc), yerr=err.(ydata_lc), fmt="d", mfc="none", capsize=2, label=L"$\mathrm{LC}$", color="royalblue" )
    xarr = Float64.(range(0.0, maximum(value.(a28t0[idx_beta])), length=100))
    yarr = cl(xarr, fit_lc.param); uwerr.(yarr)
    fill_between(xarr, value.(yarr)-err.(yarr), value.(yarr)+err.(yarr), alpha=0.2, color="royalblue")
    
    # tree-level improved
    # pt_pred = (1-(sqrt(q)/(2*3))^2 ) * 1 / (4pi^2) * log(2)
    # additive
    # ydata_ll_3l = ydata_ll .-  getindex.(hvp_3l_33_ll_s1, k) .+ value(pt_pred[k]); uwerr.(ydata_ll_3l)
    # ydata_lc_3l = ydata_lc .-  getindex.(hvp_3l_33_lc_s1, k) .+ value(pt_pred[k]); uwerr.(ydata_lc_3l)
    # multiplcative
    ydata_ll_3l = ydata_ll ./  getindex.(hvp_3l_33_ll_s2, k) .* value(pt_pred[k]); uwerr.(ydata_ll_3l)
    ydata_lc_3l = ydata_lc ./  getindex.(hvp_3l_33_lc_s2, k) .* value(pt_pred[k]); uwerr.(ydata_lc_3l)
    
    fit_ll_3l = fit_routine(cl, value.(a28t0)[idx_beta], ydata_ll_3l[idx_beta], 3)
    fit_lc_3l = fit_routine(cl, value.(a28t0)[idx_beta], ydata_lc_3l[idx_beta], 3)

    errorbar(value.(a28t0), xerr=err.(a28t0), value.(ydata_ll_3l), yerr=err.(ydata_ll_3l), fmt="s", capsize=2, label=L"$\mathrm{LL-tl}$", color="orange" )
    xarr = Float64.(range(0.0, maximum(value.(a28t0)[idx_beta]), length=100))
    yarr = cl(xarr, fit_ll_3l.param); uwerr.(yarr)
    fill_between(xarr, value.(yarr)-err.(yarr), value.(yarr)+err.(yarr), alpha=0.2, color="orange")

    errorbar(value.(a28t0), xerr=err.(a28t0), value.(ydata_lc_3l), yerr=err.(ydata_lc_3l), fmt="d", capsize=2, label=L"$\mathrm{LC-tl}$", color="royalblue" )
    xarr = Float64.(range(0.0, maximum(value.(a28t0)[idx_beta]), length=100))
    yarr = cl(xarr, fit_lc_3l.param); uwerr.(yarr)
    fill_between(xarr, value.(yarr)-err.(yarr), value.(yarr)+err.(yarr), alpha=0.2, color="royalblue")


    title(L"$\mathrm{Set \ 2}$")
    axvline(ls="dashed", color="black", lw=0.2, alpha=0.7) 
    xlabel(L"$a^2/8t_0$")
    ylabel(L"$\bar{\Pi}^{33,\mathrm{sub}}(-Q^2)$")
    # legend(loc="lower left", ncol=2)
    legend(loc="upper left", ncol=2)
    tight_layout()
    display(gcf())
    savefig(joinpath(path_plot,"std_vs_3l_imp_su3_cl_q2_$(q)_set2.pdf"))
    close("all")
    
end

##

tight_layout()
display(gcf())
close("all")

