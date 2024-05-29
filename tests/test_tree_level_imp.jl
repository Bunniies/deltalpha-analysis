######################################################################################
# This file was created by Alessandro Conigli 
# Here we compute the finite-a HVP required in the analysis of the running of δα
# The observables are then stored in a BDIO file, following the order 
########################################################################################

using Revise
using HVPobs
using PyPlot, LaTeXStrings
using OrderedCollections
using TimerOutputs
using BDIO, ADerrors
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
# Madrid scale setting
const t0sqrt_ph = uwreal([0.1439, 0.0006], "sqrtt0 [fm]") 

include("./types.jl")
include("tools.jl")
# include("data_management.jl")
include("plot_utils.jl")
include("IO_BDIO.jl")
include("../const.jl")


path_3level = "/Users/alessandroconigli/Lattice/data/HVP/tree_level"
path_bdio = "/Users/alessandroconigli/MyDrive/postdoc-mainz/projects/deltalpha/data"
path_plot = "/Users/alessandroconigli/MyDrive/postdoc-mainz/projects/deltalpha/plots"

dir_path = filter(isdir, readdir(path_bdio, join=true))


IMPR      = true
STD_DERIV = false

spectrum_path = splitpath.(vcat(filter(!isempty, [filter(x-> occursin("spectrum.bdio", x)  , readdir(dir_path[k], join=true)) for k in eachindex(dir_path)])...))
enslist = getindex.(spectrum_path, length(spectrum_path[1])-1 )

idx_su3 = findall(x-> x in ["A653","H101", "B450", "N202", "N300", "J500"], enslist)
enslist = enslist[idx_su3]
spectrum_path = spectrum_path[idx_su3]

ensinfo = EnsInfo.(["A654", "H101", "B450", "N202", "N300", "J500"])
beta_val = getfield.(ensinfo, :beta)
NENS = length(ensinfo)
Qgev = [3., 5., 9.] # Q^2
Qmgev = 9.0 # Qm^2

##
#============ OBSERVABLE ALLOCATIONS ============#
g3l_ll, g3l_lc =  read_tree_level_v33(path_3level, cons=true)
g3l_v3s03_ll, g3l_v3s03_lc = read_tree_level_v3sig03(path_3level, cons=true, massless=true)


vv         = [Corr(uwreal.(g3l_ll)[1:end-1], "tl", "G33ll") for _ in 1:NENS] # no impr
vvc        = [Corr(uwreal.(g3l_lc)[1:end-1], "tl", "G33lc") for _ in 1:NENS] # no impr
vvc_imp_tl = [Corr(uwreal.(g3l_lc)[1:end-1], "tl", "G33lc") for _ in 1:NENS] # tree-level impr cv=0.5
vv_imp_s1  = [Corr(uwreal.(g3l_ll)[1:end-1], "tl", "G33ll") for _ in 1:NENS] # impr with mainz set
vvc_imp_s1 = [Corr(uwreal.(g3l_lc)[1:end-1], "tl", "G33lc") for _ in 1:NENS] # impr with mainz set
vv_imp_s2  = [Corr(uwreal.(g3l_ll)[1:end-1], "tl", "G33ll") for _ in 1:NENS] # impr with SF set
vvc_imp_s2 = [Corr(uwreal.(g3l_lc)[1:end-1], "tl", "G33lc") for _ in 1:NENS] # impr with SF set

v3s03_ll = Corr(g3l_v3s03_ll[1:end-1], "tl", "sigmaV")
v3s03_lc = Corr(g3l_v3s03_lc[1:end-1], "tl", "sigmaVc")

obs = Vector(undef, NENS)

for (k, ens) in enumerate(ensinfo)
    println(" - Ensemble: $(ens.id)")

    obs[k] = OrderedDict()

    # obs[k]["t0"] = read_BDIO(joinpath(spectrum_path[k]), "spectrum", "t0")[1]
    beta = ens.beta

    cv_l_s1 = cv_loc(beta)
    cv_c_s1 = cv_cons(beta)
    cv_l_s2 = -0.346 #cv_loc_set2(beta)
    cv_c_s2 = cv_cons_set2(beta)
    cv_c_tl = 0.5

    # vv improvement
    improve_corr_vkvk!(vv_imp_s1[k], v3s03_ll, 2*cv_l_s1, std=STD_DERIV)
    improve_corr_vkvk!(vv_imp_s2[k], v3s03_ll, 2*cv_l_s2, std=STD_DERIV)
    # vvc improvement
    improve_corr_vkvk_cons!(vvc_imp_s1[k], v3s03_ll, v3s03_lc, cv_l_s1, cv_c_s1, std=STD_DERIV)
    improve_corr_vkvk_cons!(vvc_imp_s2[k], v3s03_ll, v3s03_lc, cv_l_s2, cv_c_s2, std=STD_DERIV)
    improve_corr_vkvk_cons!(vvc_imp_tl[k], v3s03_ll, v3s03_lc, 0.0, cv_c_tl, std=STD_DERIV)
end

##=========== COMPUTE TMR ==============#

for (k, ens) in enumerate(ensinfo)

    KRNL = krnl_dα_qhalf_sub

    Qlat = Qgev .* value.(t0sqrt_ph.^2 ./ obs[k]["t0"]) ./hc^2 *1e6
    qmlat = Qmgev * value(t0sqrt_ph^2 / obs[k]["t0"]) /hc^2 *1e6
    
    
    pi_vv, pi_vvc, pi_vvc_imp_tl, pi_vv_imp_s1, pi_vv_imp_s2, pi_vvc_imp_s1, pi_vvc_imp_s2 = [Vector{uwreal}(undef, 0)  for _ in 1:7]

    for (j,q) in enumerate(Qlat)
        try
            push!(pi3l_33_ll, sum(tmr_integrand(g3l_33_ll[k], q, KRNL, pl=true)))
            push!(pi3l_33_lc, sum(tmr_integrand(g3l_33_lc[k], q, KRNL, pl=false)))
            println("no subtracted kernel")
        catch
            println("- subtracted kernel")
            push!(pi_vv, sum(tmr_integrand_3l(vv[k].obs, q, qmlat, KRNL, pl=false, t0ens=obs[k]["t0"], wind=true, d=0.4)))
            push!(pi_vvc, sum(tmr_integrand_3l(vvc[k].obs, q, qmlat, KRNL, pl=false, t0ens=obs[k]["t0"], wind=true, d=0.4)))
            push!(pi_vvc_imp_tl, sum(tmr_integrand_3l(vvc_imp_tl[k].obs, q, qmlat, KRNL, pl=false, t0ens=obs[k]["t0"], wind=true, d=0.4)))
            push!(pi_vv_imp_s1, sum(tmr_integrand_3l(vv_imp_s1[k].obs, q, qmlat, KRNL, pl=false, t0ens=obs[k]["t0"], wind=true, d=0.4)))
            push!(pi_vvc_imp_s1, sum(tmr_integrand_3l(vvc_imp_s1[k].obs, q, qmlat, KRNL, pl=false, t0ens=obs[k]["t0"], wind=true, d=0.4)))
            push!(pi_vv_imp_s2, sum(tmr_integrand_3l(vv_imp_s2[k].obs, q, qmlat, KRNL, pl=false, t0ens=obs[k]["t0"], wind=true, d=0.4)))
            push!(pi_vvc_imp_s2, sum(tmr_integrand_3l(vvc_imp_s2[k].obs, q, qmlat, KRNL, pl=false, t0ens=obs[k]["t0"], wind=true, d=0.4)))

        end
    end

    obs[k]["pi_vv"] = value.(pi_vv)
    obs[k]["pi_vvc"] = value.(pi_vvc)
    obs[k]["pi_vvc_imp_tl"] = value.(pi_vvc_imp_tl)
    obs[k]["pi_vv_imp_s1"] = value.(pi_vv_imp_s1)
    obs[k]["pi_vvc_imp_s1"] = value.(pi_vvc_imp_s1)
    obs[k]["pi_vv_imp_s2"] = value.(pi_vv_imp_s2)
    obs[k]["pi_vvc_imp_s2"] = value.(pi_vvc_imp_s2)


end
## perturbative prediction
ptpred = []
for q in Qgev
    f(x) = treelevel_continuum_correlator.(x) .* krnl_dα_qhalf_sub.(x, q*1e6/hc^2, Qmgev*1e6/hc^2) .* window_krnl.(x, d=0.4, delta=0.15)
    res, _ = quadgk(f,0,10,rtol=1e-3)
    println(res)
    push!(ptpred, res)
end

##
a28t0 = value.(1. ./ (8 .* getindex.(obs,"t0"))) #; uwerr.(a28t0)
idx_ordered=sortperm(a28t0, rev=true)
a28t0 = a28t0[idx_ordered]

for q in eachindex(Qgev[1:1])

    plot(a28t0, getindex.(getindex.(obs, "pi_vv"), q)[idx_ordered], marker="d", color="black", label="VV")
    plot(a28t0, getindex.(getindex.(obs, "pi_vvc"), q)[idx_ordered], marker="d", color="royalblue", label="VVc")
    plot(a28t0, getindex.(getindex.(obs, "pi_vvc_imp_tl"), q)[idx_ordered], marker="d", color="red", label="VVc, improved")
    plot(a28t0, getindex.(getindex.(obs, "pi_vv_imp_s1"), q)[idx_ordered], marker="d", color="green", label="VV, np improved")
    plot(a28t0, getindex.(getindex.(obs, "pi_vvc_imp_s1"), q)[idx_ordered], marker="d", color="purple", label="VVc, np improved")
    plot(a28t0, getindex.(getindex.(obs, "pi_vv_imp_s2"), q)[idx_ordered], marker="d", color="orange", label="VV, np improved SF")
    plot(a28t0, getindex.(getindex.(obs, "pi_vvc_imp_s2"), q)[idx_ordered], marker="d", color="yellow", label="VVc, np improved SF")

    # plot(0.0, ptpred[q], color="violet", marker="d")
    xlabel(L"$a^2/8t_0$")
    ylabel(L"$\bar{\Pi}^{33, \mathrm{sub, tl}}(-Q^2)$")
    tight_layout()
    legend(ncol=2)
    ylim(0.0075, 0.0165)
    display(gcf())
    savefig(joinpath(path_plot, "su3_tree_level_std_deriv.pdf"))
    close("all")

end