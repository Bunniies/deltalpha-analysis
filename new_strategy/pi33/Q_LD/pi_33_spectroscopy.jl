# Compute the PI33(-Q^2/16) - PI33(0) contribution with a given kernel and values of the virtualities 
# and store the results in BDIO for each ensemble analysed.
# FVC are included
# this code is speficically suited for subtracted kernel computation for D200 and E250, where spectroscopy is used instead of bounding method 

using Revise, HVPobs
using HDF5
using PyPlot, LaTeXStrings
using BDIO, ADerrors, ALPHAio
using Statistics
import ADerrors: err

#================ SET UP PLOT PARAMETERS ==================#

rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
rcParams["text.usetex"] =  true
rcParams["mathtext.fontset"]  = "cm"
rcParams["font.size"] =20
rcParams["axes.labelsize"] =26
rcParams["axes.titlesize"] = 22
plt.rc("text", usetex=true) # set to true if a LaTeX installation is present


include("../../utils/const.jl")
include("../../utils/types.jl")
include("../../utils/plot_utils.jl")
include("../../utils/IO_BDIO.jl")
include("../../utils/tools.jl")

const path_corr = "/Users/alessandroconigli/MyDrive/postdoc-mainz/projects/deltalpha/corr/impr_deriv/"
const path_bdio_obs = "/Users/alessandroconigli/MyDrive/postdoc-mainz/projects/deltalpha/data"
const path_store_pi = "/Users/alessandroconigli/MyDrive/postdoc-mainz/projects/deltalpha/PIdata/impr_deriv/new_strategy/Q_LD/tmp"
const path_plot = "/Users/alessandroconigli/MyDrive/postdoc-mainz/projects/deltalpha/plots/new_strategy/ensembles"
const path_spec = "/Users/alessandroconigli/Lattice/data/HVP/corr_recostruction"

#======= PHYSICAL CONSTANTS ====================#
const Qgev = [9, 12, 15, 18] ./ 16 # Q^2/16 # 
# const Qgev = [1.0] # Q^2
const Qmgev = 9.0 # Qm^2

const KRNLsub = krnl_dα_sub # subtracted kernel
# const KRNLsub = krnl_dα # non-subtracted kernel

enslist = sort(["D200", "E250"])
ensinfo = EnsInfo.(enslist)

path_ens = vcat([filter(x-> occursin(enslist[k], basename(x)), readdir(path_corr, join=true)) for k in eachindex(enslist)]...)

path_s1 = filter(x-> occursin("set1", x), path_ens)
path_s2 = filter(x-> occursin("set2", x), path_ens)

g33_ll_s1 = Vector{Vector{uwreal}}(undef, length(enslist))
g33_lc_s1 = Vector{Vector{uwreal}}(undef, length(enslist))
g33_ll_s2 = Vector{Vector{uwreal}}(undef, length(enslist))
g33_lc_s2 = Vector{Vector{uwreal}}(undef, length(enslist))

for (k, ens) in enumerate(ensinfo)
    
    fbs1 = BDIO_open(path_s1[k], "r")
    res_s1 = Dict()
    while ALPHAdobs_next_p(fbs1)
        d = ALPHAdobs_read_parameters(fbs1)
        sz = tuple(d["size"]...)
        ks = collect(d["keys"])
        res_s1 = ALPHAdobs_read_next(fbs1, size=sz, keys=ks)
    end
    BDIO_close!(fbs1)
    g33_ll_s1[k] = res_s1["g33_ll"]
    g33_lc_s1[k] = res_s1["g33_lc"]
    
    fbs2 = BDIO_open(path_s2[k], "r")
    res_s2 = Dict()
    while ALPHAdobs_next_p(fbs2)
        d = ALPHAdobs_read_parameters(fbs2)
        sz = tuple(d["size"]...)
        ks = collect(d["keys"])
        res_s2 = ALPHAdobs_read_next(fbs2, size=sz, keys=ks)
    end
    BDIO_close!(fbs2)
    g33_ll_s2[k] = res_s2["g33_ll"]
    g33_lc_s2[k] = res_s2["g33_lc"]
end

#============== READ t0 FROM BDIO FILES =================#
dir_path = filter(x-> basename(x) in enslist, readdir(path_bdio_obs, join=true))
spectrum_path = vcat(filter(!isempty, [filter(x-> occursin("spectrum.bdio", x)  , readdir(dir_path[k], join=true)) for k in eachindex(dir_path)])...)

t0ens = Vector{uwreal}(undef, length(enslist))

for (k, ens) in enumerate(ensinfo)
    t0ens_aux = read_BDIO(spectrum_path[k], "spectrum", "t0")[1]
    t0ens[k] = t0ens_aux
end

## 
#========== RECONSTRUCT CORRELATOR =========#

g33_ll_s1_rec = similar(g33_ll_s1)
g33_lc_s1_rec = similar(g33_lc_s1)
g33_ll_s2_rec = similar(g33_ll_s2)
g33_lc_s2_rec = similar(g33_lc_s2)

for (k, ens) in enumerate(ensinfo)
    E, Z, Z_impr = get_spectr_data(path_spec, ens)

    corr_rec_ll_s1_tmp = reconstr_corr(ens, E, Z, Z_impr, impr_set="1", nmax=Dict("E250"=>4,"D200"=>2)[ens.id], IMPR=true, RENORM=true, total=false, std=true)[end]
    corr_rec_ll_s2_tmp = reconstr_corr(ens, E, Z, Z_impr, impr_set="1", nmax=Dict("E250"=>4,"D200"=>2)[ens.id], IMPR=true, RENORM=true, total=false, std=true)[end]
    
    ratio_lc_ll_s1 = g33_lc_s1[k] ./ g33_ll_s1[k]
    ratio_lc_ll_s2 = g33_lc_s2[k] ./ g33_ll_s2[k]

    corr_rec_lc_s1_tmp = corr_rec_ll_s1_tmp .* ratio_lc_ll_s1 
    corr_rec_lc_s2_tmp = corr_rec_ll_s2_tmp .* ratio_lc_ll_s2 
    
    uwerr.(g33_ll_s1[k]); uwerr.(g33_ll_s2[k]); uwerr.(corr_rec_ll_s1_tmp); uwerr.(corr_rec_ll_s2_tmp)
    trec_s1 = findfirst_uninterrupted(err.(corr_rec_ll_s1_tmp) .< err.(g33_ll_s1[k]))
    trec_s2 = findfirst_uninterrupted(err.(corr_rec_ll_s1_tmp) .< err.(g33_ll_s1[k]))
    println("t_s1 :", trec_s1, " t_s2: ", trec_s2)

    g33_ll_s1_rec[k] = vcat(g33_ll_s1[k][1:trec_s1], corr_rec_ll_s1_tmp[trec_s1+1:end])
    g33_lc_s1_rec[k] = vcat(g33_lc_s1[k][1:trec_s1], corr_rec_lc_s1_tmp[trec_s1+1:end])
    g33_ll_s2_rec[k] = vcat(g33_ll_s2[k][1:trec_s1], corr_rec_ll_s2_tmp[trec_s1+1:end])
    g33_lc_s2_rec[k] = vcat(g33_lc_s2[k][1:trec_s1], corr_rec_lc_s2_tmp[trec_s1+1:end])

end

##
#=========== COMPUTE TMR ==============#
# set 1 improvement coefficients
pi_33_ll_s1 = [Vector{uwreal}(undef, length(Qgev)) for k in eachindex(ensinfo)] 
pi_33_lc_s1 = [Vector{uwreal}(undef, length(Qgev)) for k in eachindex(ensinfo)]

pi_33_ll_s1_rec = [Vector{uwreal}(undef, length(Qgev)) for k in eachindex(ensinfo)] 
pi_33_lc_s1_rec = [Vector{uwreal}(undef, length(Qgev)) for k in eachindex(ensinfo)]
# set 2 improvement coefficients
pi_33_ll_s2  = [Vector{uwreal}(undef, length(Qgev)) for k in eachindex(ensinfo)] 
pi_33_lc_s2  = [Vector{uwreal}(undef, length(Qgev)) for k in eachindex(ensinfo)]

pi_33_ll_s2_rec  = [Vector{uwreal}(undef, length(Qgev)) for k in eachindex(ensinfo)] 
pi_33_lc_s2_rec  = [Vector{uwreal}(undef, length(Qgev)) for k in eachindex(ensinfo)]

for (k, ens) in enumerate(ensinfo)
    println("Ensemble: ", ens.id)

    Qlat = Qgev .* t0sqrt_ph.^2 ./ t0ens[k] ./hc^2 *1e6
    qmlat = Qmgev  * t0sqrt_ph^2 / t0ens[k] /hc^2 *1e6

    for (j,q) in enumerate(Qlat)
        x0 = Float64.(collect(0:length(g33_lc_s1[k])-1))
        kk = KRNLsub.(x0, q, qmlat)

        # pi_33_ll_s1[k][j] = tmr_integrand(g33_ll_s1[k], q, qmlat, KRNLsub, pl=false, t0ens=t0ens[k])
        # pi_33_lc_s1[k][j] = tmr_integrand(g33_lc_s1[k], q, qmlat, KRNLsub, pl=false, t0ens=t0ens[k])
        # pi_33_ll_s2[k][j] = tmr_integrand(g33_ll_s2[k], q, qmlat, KRNLsub, pl=false, t0ens=t0ens[k])
        # pi_33_lc_s2[k][j] = tmr_integrand(g33_lc_s2[k], q, qmlat, KRNLsub, pl=false, t0ens=t0ens[k])

        pi_33_ll_s1_rec[k][j] = boundingMethod(g33_ll_s1_rec[k], ens, kk, "33", pl=true)
        pi_33_lc_s1_rec[k][j] = boundingMethod(g33_lc_s1_rec[k], ens, kk, "33", pl=false)
        pi_33_ll_s2_rec[k][j] = boundingMethod(g33_ll_s2_rec[k], ens, kk, "33", pl=true)
        pi_33_lc_s2_rec[k][j] = boundingMethod(g33_lc_s2_rec[k], ens, kk, "33", pl=false)

    end
end

##
#============== SAVE results to BDIO ===============#
@info("Saving PI 33 results in BDIO")
io = IOBuffer()
write(io, "PI 33 Q_LD. ")
fb = ALPHAdobs_create(joinpath(path_store_pi, "PI_33_D200_E250.bdio"), io)

for (k, ens) in enumerate(ensinfo)
    extra = Dict{String, Any}("Ens" => ens.id)
    data = Dict{String, Array{uwreal}}(
        "pi33_ll_s1" => pi_33_ll_s1_rec[k],
        "pi33_lc_s1" => pi_33_lc_s1_rec[k],
        "pi33_ll_s2" => pi_33_ll_s2_rec[k],
        "pi33_lc_s2" => pi_33_lc_s2_rec[k]
    )

    ALPHAdobs_write(fb, data, extra=extra)
end
ALPHAdobs_close(fb)
println("# Saving complete!")

##
#========= TEST READING =========#
fb = BDIO_open(joinpath(path_store_pi, "PI_33_D200_E250.bdio"), "r")
# fb = BDIO_open(pptest, "r")
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


#########################################
## TEST WITH CORRELATOR RECOSTRUCTION
#########################################

pi_33_ll_s1 = [Vector{uwreal}(undef, length(Qgev)) for k in eachindex(ensinfo)] 
pi_33_ll_s1_rec = [Vector{uwreal}(undef, length(Qgev)) for k in eachindex(ensinfo)] 

# set 2 improvement coefficients
pi_33_ll_s2  = [Vector{uwreal}(undef, length(Qgev)) for k in eachindex(ensinfo)] 
pi_33_ll_s2_rec  = [Vector{uwreal}(undef, length(Qgev)) for k in eachindex(ensinfo)] 

for (k,ens) in enumerate(ensinfo)
    Qlat = Qgev .* t0sqrt_ph.^2 ./ t0ens[k] ./hc^2 *1e6
    qmlat = Qmgev  * t0sqrt_ph^2 / t0ens[k] /hc^2 *1e6

    for (j,q) in enumerate(Qlat)
        x0 = Float64.(collect(0: length(g33_lc_s1[k])-1))
        kk = KRNLsub.(x0, q, qmlat)

        pi_33_ll_s1[k][j], lowBound, upBound, ttime = boundingMethod(g33_ll_s1[k], ens, kk, "33", pl=false, data=true)
        #pi_33_ll_s2[k][j] = boundingMethod(g33_ll_s2[k], ens, kk, "33", pl=false)

        pi_33_ll_s1_rec[k][j], lowBound, upBound, ttime = boundingMethod(g33_ll_s1_rec[k], ens, kk, "33", pl=false, data=true)
        #pi_33_ll_s2_rec[k][j] = boundingMethod(g33_ll_s2_rec, ens, kk, "33", pl=false)
        #pi_33_ll_s1_rec[k][j] = tmr_integrand(g33_ll_s1_rec, q, qmlat, KRNLsub, pl=false, t0ens=t0ens[k])
        #pi_33_ll_s2_rec[k][j] = tmr_integrand(g33_ll_s2_rec, q, qmlat, KRNLsub, pl=false, t0ens=t0ens[k])

        uwerr.(lowBound); uwerr.(upBound)
        uwerr(pi_33_ll_s1_rec[k][j])
        uwerr(pi_33_ll_s1[k][j])
        fig = figure(figsize=(10,7.))
        errorbar(ttime .*value(a(ens.beta)), value.(lowBound), err.(lowBound), fmt="d", mfc="none", color="#009E73", label=L"$M_{\mathrm{eff}}$", capsize=2)
        errorbar(ttime .*value(a(ens.beta)), value.(upBound), err.(upBound), fmt="d", mfc="none", color="#D55E00", label=L"$E_{\pi\pi}$", capsize=2)
        fill_between(ttime .*value(a(ens.beta)), value(pi_33_ll_s1_rec[k][j])- err(pi_33_ll_s1_rec[k][j]), value(pi_33_ll_s1_rec[k][j])+ err(pi_33_ll_s1_rec[k][j]), alpha=0.5, color="#E69F00")
        axhline(xmin=0.5, xmax=1, value(pi_33_ll_s1[k][j]) + err(pi_33_ll_s1[k][j]), ls="dashed", color="#0072B2")
        axhline(xmin=0.5, xmax=1, value(pi_33_ll_s1[k][j]) - err(pi_33_ll_s1[k][j]), ls="dashed", color="#0072B2")
        if ens.id == "D200"
            axvline(32 * value(a(ens.beta)), ls="dashed", color="gray")
        else
            axvline(39 * value(a(ens.beta)), ls="dashed", color="gray")
        end
        ylim(value(pi_33_ll_s1[k][j]) - 5* err(pi_33_ll_s1[k][j]), value(pi_33_ll_s1[k][j]) + 5* err(pi_33_ll_s1[k][j]))
        if ens.id == "D200"
            xlim(0.5,3.5)
        else
            xlim(2.0, 5.0)
        end
        legend()
        xlabel(L"$t \ [\mathrm{fm}]$")
        ylabel(L"$\bar\Pi^{(3,3)}(Q^2/16)$")
        tight_layout()
        display(fig)
        savefig(joinpath(path_plot, ens.id, "bm_vs_rec_q$(Qgev[j]*16).pdf"))
        close("all")
    end
end