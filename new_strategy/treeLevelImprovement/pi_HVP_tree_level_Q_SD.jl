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

include("../utils/const.jl")
include("../utils/types.jl")
include("../utils/tools.jl")
include("../utils/plot_utils.jl")
include("../utils/IO_BDIO.jl")
# include("../utils/data_management.jl")


path_3level = "/Users/alessandroconigli/Lattice/data/HVP/tree_level"
path_store_3l = "/Users/alessandroconigli/MyDrive/postdoc-mainz/projects/deltalpha/PIdata/impr_deriv/new_strategy/Q_SD"
path_plot = "/Users/alessandroconigli/MyDrive/postdoc-mainz/projects/deltalpha/plots"
path_bdio = "/Users/alessandroconigli/MyDrive/postdoc-mainz/projects/deltalpha/data"

dir_path = filter(isdir, readdir(path_bdio, join=true))


IMPR      = true
IMPR_SET  = "1" # either "1" or "2"
RENORM    = false
STD_DERIV = false


enslist = sort([ "H101", "H102", "N101", "C101", "C102", "D150",
         "B450", "N451", "D450", "D451", "D452", "N452",
         "N202", "N203", "N200", "D251", "D200", "D201", "E250",
          "J307", "J306", "J303", "J304", "E300", "F300",
         "J500", "J501"])

ensinfo = EnsInfo.(enslist)
# ensinfo = EnsInfo.(["H101"])

dir_path = filter(x-> basename(x) in enslist, readdir(path_bdio, join=true))
spectrum_path = vcat(filter(!isempty, [filter(x-> occursin("spectrum.bdio", x)  , readdir(dir_path[k], join=true)) for k in eachindex(dir_path)])...)


beta_val = getfield.(ensinfo, :beta)
NENS = length(ensinfo)
# Qgev = [3., 5., 9.] # Q^2
const Qgev = [4, 5, 6, 7, 8, 9, 12]  # Q^2 / 16 # additional very high values

Qmgev = 9.0 # Qm^2

##
#============ OBSERVABLE ALLOCATIONS ============#
g3l_ll, g3l_lc =  read_tree_level_v33(path_3level, cons=true)
g3l_v3s03_ll, g3l_v3s03_lc = read_tree_level_v3sig03(path_3level, cons=true, massless=true)


g3l_33_ll = [uwreal.(g3l_ll) for _ in 1:NENS]
g3l_33_lc = [uwreal.(g3l_lc) for _ in 1:NENS]

obs = Vector(undef, NENS)

for (k, ens) in enumerate(ensinfo)
    println(" - Ensemble: $(ens.id)")

    obs[k] = OrderedDict()
    obs[k]["t0"] = read_BDIO(joinpath(spectrum_path[k]), "spectrum", "t0")[1]

    if IMPR
        beta = ens.beta
        if IMPR_SET == "1"
            cv_l = cv_loc(beta)
            cv_c = cv_cons(beta)
        elseif IMPR_SET =="2"
            cv_l = cv_loc_set2(beta)
            cv_c = cv_cons_set2(beta)
        end
    end
    improve_corr_vkvk!(g3l_33_ll[k], g3l_v3s03_ll, 2*cv_l, std=STD_DERIV, treelevel=true)
    improve_corr_vkvk_cons!(g3l_33_lc[k], g3l_v3s03_ll, g3l_v3s03_lc, cv_l, cv_c, std=STD_DERIV, treelevel=true)

end


##=========== COMPUTE TMR ==============#

for (k, ens) in enumerate(ensinfo)

    KRNL = krnl_dα_qhalf_sub

    Qlat = Qgev .* value.(t0sqrt_ph.^2 ./ obs[k]["t0"]) ./hc^2 *1e6
    qmlat = Qmgev * value(t0sqrt_ph^2 / obs[k]["t0"]) /hc^2 *1e6
    
    
    pi3l_33_ll, pi3l_33_lc = [Vector{uwreal}(undef, 0)  for _ in 1:2]

    for (j,q) in enumerate(Qlat)
        try
            push!(pi3l_33_ll, sum(tmr_integrand(g3l_33_ll[k], q, KRNL, pl=true)))
            push!(pi3l_33_lc, sum(tmr_integrand(g3l_33_lc[k], q, KRNL, pl=false)))
        catch
            println("- subtracted kernel")
            push!(pi3l_33_ll, sum(tmr_integrand_3l(g3l_33_ll[k], q, qmlat, KRNL, pl=false, t0ens=obs[k]["t0"], wind=true, d=0.4)))
            push!(pi3l_33_lc, sum(tmr_integrand_3l(g3l_33_lc[k], q, qmlat, KRNL, pl=false, t0ens=obs[k]["t0"], wind=true, d=0.4)))
        end
    end

    obs[k]["pi3l_33_ll"] = pi3l_33_ll
    obs[k]["pi3l_33_lc"] = pi3l_33_lc

end

##
#======== SAVE TO BDIO WITH ALPHAIO =========#
@info("Saving tree level improvement TMR")
io = IOBuffer()
write(io, "TMR 3 level improvement set 1")
fb = ALPHAdobs_create(joinpath(path_store_3l, "treeLevelSet1.bdio"), io)

for(k, ens) in enumerate(ensinfo)
    extra = Dict{String, Any}("Ens" => ens.id)
    data = Dict{String, Array{uwreal}}(
        "3l_33_ll" => obs[k]["pi3l_33_ll"],
        "3l_33_lc" => obs[k]["pi3l_33_lc"]
    )
    ALPHAdobs_write(fb, data, extra=extra)
end
ALPHAdobs_close(fb)
println("# Saving complete!")

##
#========= TEST READING =========#
fb = BDIO_open(joinpath(path_store_3l,  "treeLevelSet2.bdio"), "r")
res = Dict()
while ALPHAdobs_next_p(fb)
    d = ALPHAdobs_read_parameters(fb)
    nobs = d["nobs"]
    dims = d["dimensions"]
    sz = tuple(d["size"]...)
    extra = d["extra"]
    ks = collect(d["keys"])
    res[extra["Ens"]] = ALPHAdobs_read_next(fb, size=sz, keys=ks)
end
BDIO_close!(fb)

#======== SAVE TO BDIO OLD FORMAT =========#
##

[pop!(obs[k], "t0") for k in eachindex(obs)]
for (k, o) in enumerate(obs)
    ens = ensinfo[k].id
    p = joinpath(path_bdio, ens,  string(ens, "_tree_level_set$(IMPR_SET).bdio"))
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


###########################################
## testing CL in SU(3) ensembles
##########################################

idx_su3 = findall(x-> x in ["H101", "B450", "N202", "N300", "J500"], enslist)
enslist[idx_su3]
a28t0 = 1. ./ (8 .* getindex.(obs[idx_su3],"t0")); uwerr.(a28t0)

@. cl(x,p) = p[1] + p[2]*x + p[3]*x^(3/2.) #+ p[4]*x^2
ph_val_3l = Vector{uwreal}(undef, 0)

for (k,q) in enumerate(Qgev)
    ydata_ll = [obs[j]["pi3l_33_ll"][k] for j in idx_su3]; uwerr.(ydata_ll)
    ydata_lc = [obs[j]["pi3l_33_lc"][k] for j in idx_su3]; uwerr.(ydata_lc)

    fit = fit_routine(cl, value.(a28t0), ydata_ll, 3)
    ph_val = fit.param[1]; uwerr(ph_val)
    xarr = Float64.(range(0.0, maximum(value.(a28t0)), length=100))
    yarr = cl(xarr, fit.param); uwerr.(yarr)

    errorbar(value.(a28t0), xerr=err.(a28t0), value.(ydata_ll), yerr=err.(ydata_ll), fmt="s", capsize=2, label=L"$\mathrm{LL}$", color="orange", mfc="none" )
    fill_between(xarr, value.(yarr)-err.(yarr), value.(yarr)+err.(yarr), alpha=0.2, color="orange" )

    fit = fit_routine(cl, value.(a28t0), ydata_lc, 3)
    ph_val = fit.param[1]; uwerr(ph_val)
    push!(ph_val_3l, ph_val)
    xarr = Float64.(range(0.0, maximum(value.(a28t0)), length=100))
    yarr = cl(xarr, fit.param); uwerr.(yarr)

    errorbar(value.(a28t0), xerr=err.(a28t0), value.(ydata_lc), yerr=err.(ydata_lc), fmt="s", capsize=2, label=L"$\mathrm{LC}$", color="navy", mfc="none" )
    fill_between(xarr, value.(yarr)-err.(yarr), value.(yarr)+err.(yarr), alpha=0.2, color="navy" )
    
    # PT prediction
    pt_pred = (1-(sqrt(q)/(2*3))^2 ) * 1 / (4pi^2) * log(2)
    #plot(0.0, pt_pred, label="PT", marker="o", color="red")
    axvline(ls="dashed", color="black", lw=0.2, alpha=0.7) 
    xlabel(L"$a^2/8t_0$")
    ylabel(L"$\bar{\Pi}^{33, \mathrm{sub, tl}}(-Q^2)$")
    # ylabel(L"$\bar{\Pi}^{33,\mathrm{sub}}(-Q^2)$")
    legend()
    # ylim(0.0115,0.02)
    # ylim(0.0075, 0.0165)
    title("Set $(IMPR_SET)")
    tight_layout()
    xlim(-0.002, 0.05)
    display(gcf())
    savefig(joinpath(path_plot, "tree-level", "tree_level_su3_cl_q2_$(q)_set_$(IMPR_SET).pdf" ))
    close("all")
end

## save 3l phys val in bdio
p = joinpath(path_bdio, string("tree_level_cl_phys_val.bdio"))

fb = BDIO_open(p, "w", "3l-cl")
uinfo = 0
[write_uwreal(v, fb, uinfo) for v in ph_val_3l]
BDIO_close!(fb)

##
val_read = []
fb = BDIO_open(p, "r")
BDIO_seek!(fb)
if BDIO_get_uinfo(fb) == 0
    push!(val_read, read_uwreal(fb))
end
while BDIO_seek!(fb, 2)
    if BDIO_get_uinfo(fb) == 0 
        push!(val_read, read_uwreal(fb))
    end
end
BDIO_close!(fb)