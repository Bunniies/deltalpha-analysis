# Compute the PI33(-Q^2/4) - PI33(0) contribution with a given kernel and values of the virtualities 
# and store the results in BDIO for each ensemble analysed.
# FVC are included
# this code is speficically suited for subtracted kernel computation 
using Revise, HVPobs
using HDF5
using PyPlot, LaTeXStrings
using BDIO, ADerrors, ALPHAio
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
# include("./func_comb.jl")

const path_corr = "/Users/alessandroconigli/MyDrive/postdoc-mainz/projects/deltalpha/corr/impr_deriv/"
# const path_corr = "/Users/alessandroconigli/MyDrive/postdoc-mainz/projects/deltalpha/corr/std_deriv/"
const path_bdio_obs = "/Users/alessandroconigli/MyDrive/postdoc-mainz/projects/deltalpha/data"
const path_store_pi = "/Users/alessandroconigli/MyDrive/postdoc-mainz/projects/deltalpha/PIdata/impr_deriv/low_q_kernel/scale_error_artificial/tmp"
const path_fvc  = "/Users/alessandroconigli/Lattice/data/HVP/FSE"
const path_plot = "/Users/alessandroconigli/MyDrive/postdoc-mainz/projects/deltalpha/plots/ensembles"
const path_spec = "/Users/alessandroconigli/Lattice/data/HVP/corr_recostruction"

#======= PHYSICAL CONSTANTS ====================#
const Qgev = [0.05, 0.1, 0.4, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0] ./ 4 # Q^2
# const Qgev = [1.0] # Q^2
const Qmgev = 9.0 # Qm^2

const KRNLsub = krnl_dα_sub # subtracted kernel
# const KRNLsub = krnl_dα # non-subtracted kernel


# enslist = sort([ "H101", "H102", "N101", "C101", "C102", "D150",
        #  "B450", "N451", "D450", "D451", "D452",
        #  "N202", "N203", "N200", "D251", "D200", "D201", "E250",
        #   "J307", "J306", "J303", "J304", "E300", "F300",
        #  "J500", "J501"])
enslist = sort(["D200"])
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

##
#============== READ t0 FROM BDIO FILES =================#
dir_path = filter(x-> basename(x) in enslist, readdir(path_bdio_obs, join=true))
spectrum_path = vcat(filter(!isempty, [filter(x-> occursin("spectrum.bdio", x)  , readdir(dir_path[k], join=true)) for k in eachindex(dir_path)])...)

t0ens = Vector{uwreal}(undef, length(enslist))

for (k, ens) in enumerate(ensinfo)
    t0ens_aux = read_BDIO(spectrum_path[k], "spectrum", "t0")[1]
    t0ens[k] = t0ens_aux
end


##
#=========== COMPUTE FVC  ==============#
@info("Computing Finite-Volume Corrections")
pi_fvc_pi = [Vector{uwreal}(undef, length(Qgev)) for k in eachindex(ensinfo)]
pi_fvc_k  = [Vector{uwreal}(undef, length(Qgev)) for k in eachindex(ensinfo)]

for (k,ens) in enumerate(ensinfo)
    println("Ensemble: $(ens.id)")

    # Qlat = Qgev .* t0sqrt_ph.^2 ./ t0ens[k] ./hc^2 *1e6
    # qmlat = Qmgev  * t0sqrt_ph^2 / t0ens[k] /hc^2 *1e6

    Qlat = Qgev .* t0sqrt_ph.^2 ./ t0ens[k] ./hc^2 *1e6
    qmlat = Qmgev * t0sqrt_ph^2 / t0ens[k] /hc^2 *1e6

    fvc_raw_pi = get_fvc(joinpath(path_fvc, "JKMPI" ), ens.id)
    fvc_pi = vcat(sum(fvc_raw_pi, dims=1)...)
    
    fvc_raw_k = get_fvc(joinpath(path_fvc, "JKMK" ), ens.id)
    fvc_k = vcat(sum(fvc_raw_k, dims=1)...)
    
    T = Float64.(collect(1:length(fvc_pi)))

    for (j,q) in enumerate(Qlat)
        krnl = KRNLsub.(T,q, qmlat)
        pi_fvc_pi[k][j]  = abs(sum(fvc_pi .* krnl))
        pi_fvc_k[k][j]   = abs(sum(fvc_k  .* krnl))
    end
end

##
#=========== COMPUTE TMR ==============#
# set 1 improvement coefficients
pi_33_ll_s1 = [Vector{uwreal}(undef, length(Qgev)) for k in eachindex(ensinfo)] 
pi_33_lc_s1 = [Vector{uwreal}(undef, length(Qgev)) for k in eachindex(ensinfo)]

# set 2 improvement coefficients
pi_33_ll_s2  = [Vector{uwreal}(undef, length(Qgev)) for k in eachindex(ensinfo)] 
pi_33_lc_s2  = [Vector{uwreal}(undef, length(Qgev)) for k in eachindex(ensinfo)]
pi_33_lc_s2_data  = [Vector{Vector{uwreal}}(undef, length(Qgev)) for k in eachindex(ensinfo)]


for (k, ens) in enumerate(ensinfo)
    println("Ensemble: ", ens.id)

    Qlat = Qgev .* t0sqrt_ph.^2 ./ t0ens[k] ./hc^2 *1e6
    qmlat = Qmgev  * t0sqrt_ph^2 / t0ens[k] /hc^2 *1e6

    for (j,q) in enumerate(Qlat)
        x0 = Float64.(collect(0: length(g33_lc_s1[k])-1))
        kk = KRNLsub.(x0, q, qmlat )
        # kk = KRNLsub.(x0, q)
        
        if j == 12
            ss = joinpath(path_plot, ens.id)
            qqq = "9"
        else
            ss = nothing
            qqq = ""
        end
        pi_33_ll_s1[k][j] = boundingMethod(g33_ll_s1[k], ens, kk, "33", path_pl=ss, qval=qqq)
        pi_33_lc_s1[k][j] = boundingMethod(g33_lc_s1[k], ens, kk, "33")
        pi_33_ll_s2[k][j] = boundingMethod(g33_ll_s2[k], ens, kk, "33")
        pi_33_lc_s2[k][j] = boundingMethod(g33_lc_s2[k], ens, kk, "33")

        # pi_33_ll_s1[k][j] = tmr_integrand(g33_ll_s1[k], q, qmlat, KRNLsub, pl=false, t0ens=t0ens[k])
        # pi_33_lc_s1[k][j] = tmr_integrand(g33_lc_s1[k], q, qmlat, KRNLsub, pl=false, t0ens=t0ens[k])
        # pi_33_ll_s2[k][j] = tmr_integrand(g33_ll_s2[k], q, qmlat, KRNLsub, pl=false, t0ens=t0ens[k])
        # pi_33_lc_s2[k][j] = tmr_integrand(g33_lc_s2[k], q, qmlat, KRNLsub, pl=false, t0ens=t0ens[k])
        #_, pi_33_lc_s2_data[k][j] = tmr_integrand(g33_lc_s2[k], q, qmlat, KRNLsub, pl=false, t0ens=t0ens[k], data=true)
    end 
end

##
#============ ADD FVC ============#
@info("Adding Finite-Volume Corrections")
for k in eachindex(ensinfo)

    if ensinfo[k].kappa_l == ensinfo[k].kappa_s
        pi_33_ll_s1[k] .+= 1.5 * pi_fvc_pi[k]
        pi_33_lc_s1[k] .+= 1.5 * pi_fvc_pi[k]
        pi_33_ll_s2[k] .+= 1.5 * pi_fvc_pi[k]
        pi_33_lc_s2[k] .+= 1.5 * pi_fvc_pi[k]   
    else 
        pi_33_ll_s1[k] .+= pi_fvc_pi[k] .+ pi_fvc_k[k]
        pi_33_lc_s1[k] .+= pi_fvc_pi[k] .+ pi_fvc_k[k]
        pi_33_ll_s2[k] .+= pi_fvc_pi[k] .+ pi_fvc_k[k]
        pi_33_lc_s2[k] .+= pi_fvc_pi[k] .+ pi_fvc_k[k]
    end

end

##
#============== SAVE results to BDIO ===============#
@info("Saving PI 33 results in BDIO")
io = IOBuffer()
write(io, "PI delta 33 low q. ")
fb = ALPHAdobs_create(joinpath(path_store_pi, "PI_33_N452_F300.bdio"), io)

for (k, ens) in enumerate(ensinfo)
    extra = Dict{String, Any}("Ens" => ens.id)
    data = Dict{String, Array{uwreal}}(
        "pi33_ll_s1" => pi_33_ll_s1[k],
        "pi33_lc_s1" => pi_33_lc_s1[k],
        "pi33_ll_s2" => pi_33_ll_s2[k],
        "pi33_lc_s2" => pi_33_lc_s2[k]
    )

    ALPHAdobs_write(fb, data, extra=extra)
end
ALPHAdobs_close(fb)
println("# Saving complete!")

##
#========= TEST READING =========#
pptest = "/Users/alessandroconigli/MyDrive/postdoc-mainz/projects/deltalpha/PIdata/impr_deriv/low_q_kernel/scale_error_artificial/PI_33.bdio"
# fb = BDIO_open(joinpath(path_store_pi, "multi_mom", "PI_33.bdio"), "r")
fb = BDIO_open(pptest, "r")
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


##
#============= TEST COMPUTING DERIVATIVES WRT t0_ph =============#
ddtest = derivative(pi_33_ll_SD_s1[1][3], t0sqrt_ph) 

## test with rho meson mass
ensinfo
kk = 2
Thalf = Int64(length(g33_ll_s1[kk]))
@. const_fit(x,p) = p[2] * (exp(-p[1]*(x)) + exp(-p[1]*(Thalf-x)) )# + p[4] * (exp(-p[3]*(x)) + exp(-p[3]*(T-x)) )
#@. const_fit(x,p) = p[2] * (exp(-p[1]*(x)) )# + p[4] * (exp(-p[3]*(x)) + exp(-p[3]*(T-x)) )

xdata = collect(30:Thalf-30)  
fit = fit_routine(const_fit, xdata, g33_ll_s1[kk][30:Thalf-30], 2 )
uwerr(fit.param[1])
println(fit.param[1])
ensinfo[kk]
##
meff(g33_ll_s1[kk], [30,45], pl=true)


## plot comulative sum

fig = figure()
idx=4
cum = cumsum(pi_33_lc_s2_data[1][idx])
xx = collect(1:length(cum)) .* value(a(3.55))
uwerr.(cum)
errorbar(xx, value.(cum), err.(cum), fmt="s", ms=10, mfc="none", color="forestgreen")   
fill_between(xx, value(pi_33_lc_s2[1][idx])-err(pi_33_lc_s2[1][idx]), value(pi_33_lc_s2[1][idx])+err(pi_33_lc_s2[1][idx]), color="red",alpha=0.5)

axvline(55*value(a(3.55)), ls="dashed", color="gray")
xlabel(L"$t \ [\mathrm{fm}]$")
ylabel(L"$\bar{\Pi}^{(3,3)}(Q^2/4)$")
title(L"$Q^2=$" * string(Qgev[idx]*4) * L"$\ \mathrm{GeV}^2$")
display(gcf())
tight_layout()
savefig(joinpath(path_plot,"E250/cum_sum_E250_Qgev$(Qgev[idx]*4).pdf"))
close("all")

#########################################
## TEST WITH CORRELATOR RECOSTRUCTION
#########################################
using HDF5
using Statistics
E, Z, Z_impr = get_spectr_data(path_spec, ensinfo[1])

g33_ll_s1_rec = reconstr_corr(ensinfo[1],E, Z, Z_impr, impr_set="1", nmax=Dict("E250"=>4,"D200"=>2)[ensinfo[1].id], IMPR=true, RENORM=true, total=false, std=true)
uwerr.(g33_ll_s1_rec[end])
uwerr.(g33_ll_s1[1])
trec = findfirst_uninterrupted(err.(g33_ll_s1_rec[end]) .< err.(g33_ll_s1[1][1:96]))
vcat(g33_ll_s1[1][1:trec-1], g33_ll_s1_rec[end][trec:end], g33_ll_s1[1][97:end])


##
#=========== COMPUTE TMR WITH SPECTROSCOPY ==============#
# set 1 improvement coefficients
pi_33_ll_s1 = [Vector{uwreal}(undef, length(Qgev)) for k in eachindex(ensinfo)] 
pi_33_ll_s1_rec = [Vector{uwreal}(undef, length(Qgev)) for k in eachindex(ensinfo)] 

# set 2 improvement coefficients
pi_33_ll_s2  = [Vector{uwreal}(undef, length(Qgev)) for k in eachindex(ensinfo)] 
pi_33_ll_s2_rec  = [Vector{uwreal}(undef, length(Qgev)) for k in eachindex(ensinfo)] 


for (k, ens) in enumerate(ensinfo)
    println("Ensemble: ", ens.id)

    Qlat = Qgev .* t0sqrt_ph.^2 ./ t0ens[k] ./hc^2 *1e6
    qmlat = Qmgev  * t0sqrt_ph^2 / t0ens[k] /hc^2 *1e6

    E, Z, Z_impr = get_spectr_data(path_spec, ens)

    g33_ll_s1_rec_aux = reconstr_corr(ens, E, Z, Z_impr, impr_set="1", nmax=Dict("E250"=>4,"D200"=>2)[ens.id], IMPR=true, RENORM=true, total=false, std=true)[end]
    g33_ll_s2_rec_aux = reconstr_corr(ens, E, Z, Z_impr, impr_set="2", nmax=Dict("E250"=>4,"D200"=>2)[ens.id], IMPR=true, RENORM=true, total=false, std=true)[end]

    uwerr.(g33_ll_s1_rec_aux); uwerr.(g33_ll_s2_rec_aux)
    uwerr.(g33_ll_s1[k]); uwerr.(g33_ll_s2[k])

    trec_s1 = findfirst_uninterrupted(err.(g33_ll_s1_rec_aux) .< err.(g33_ll_s1[k]))
    trec_s2 = findfirst_uninterrupted(err.(g33_ll_s2_rec_aux) .< err.(g33_ll_s2[k]))
    println("set1: ", trec_s1, " set2: ", trec_s2)
    g33_ll_s1_rec = vcat(g33_ll_s1[k][1:trec_s1], g33_ll_s1_rec_aux[trec_s1+1:end])
    g33_ll_s2_rec = vcat(g33_ll_s2[k][1:trec_s2], g33_ll_s2_rec_aux[trec_s2+1:end])

    for (j,q) in enumerate(Qlat)
        x0 = Float64.(collect(0: length(g33_lc_s1[k])-1))
        kk = KRNLsub.(x0, q, qmlat )
        # kk = KRNLsub.(x0, q)
        if j in [2,3,5,7,9,10,11,12]
            #continue
        end
        if j == 12
            ss = joinpath(path_plot, ens.id)
            qqq = "9"
        else
            ss = nothing
            qqq = ""
        end
        pi_33_ll_s1[k][j], lowBound, upBound, ttime = boundingMethod(g33_ll_s1[k], ens, kk, "33", pl=false, data=true)
        pi_33_ll_s2[k][j] = boundingMethod(g33_ll_s2[k], ens, kk, "33", pl=false)

        pi_33_ll_s1_rec[k][j], lowBound, upBound, ttime = boundingMethod(g33_ll_s1_rec, ens, kk, "33", pl=false, data=true)
        pi_33_ll_s2_rec[k][j] = boundingMethod(g33_ll_s2_rec, ens, kk, "33", pl=false)
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
        axvline(trec_s1 * value(a(ens.beta)), ls="dashed", color="gray")

        ylim(value(pi_33_ll_s1[k][j]) - 5* err(pi_33_ll_s1[k][j]), value(pi_33_ll_s1[k][j]) + 5* err(pi_33_ll_s1[k][j]))
        xlim(0.5,3.5)
        legend()
        xlabel(L"$t \ [\mathrm{fm}]$")
        ylabel(L"$\bar\Pi^{(3,3)}(Q^2/4)$")
        tight_layout()
        display(fig)
        savefig(joinpath(path_plot, ens.id, "bm_vs_rec_q$(Qgev[j]*4).pdf"))
        close("all")
    end 

end

## PLOT RES from BM vs RECONSTRUCTION

uwerr.(pi_33_ll_s2[1])
uwerr.(pi_33_ll_s2_rec[1])

fig = figure(figsize=(10,7.))
errorbar(Qgev, value.(pi_33_ll_s2[1]), 20 .*err.(pi_33_ll_s2[1]), fmt="d", ms=10, mfc="none", capsize=2, color="#D55E00",label="BM")
errorbar(Qgev.+0.04, value.(pi_33_ll_s2_rec[1]), 20 .*err.(pi_33_ll_s2_rec[1]), fmt="^", ms=10, mfc="none", capsize=2, color="#0072B2", label="Recostruction")

# xlim(0.9,2.4)
# ylim(0.03, 0.045)
ylabel(L"$\bar\Pi^{(3,3)}(Q^2/4)$")
xlabel(L"$Q^2/4 \ [\mathrm{GeV^2}]$")
title("Set 2")
legend()
tight_layout()
display(fig)
savefig(joinpath(path_plot, ensinfo[1].id, "comparison_bm_rec_set2.pdf"))
close("all")

## RECONSTRUCT THE TMR

# set 1 improvement coefficients
pi_33_ll_s1 = [Vector{uwreal}(undef, length(Qgev)) for k in eachindex(ensinfo)] 
pi_33_ll_s1_rec = [Vector{uwreal}(undef, length(Qgev)) for k in eachindex(ensinfo)] 

for (k, ens) in enumerate(ensinfo)
    println("Ensemble: ", ens.id)

    Qlat = Qgev .* t0sqrt_ph.^2 ./ t0ens[k] ./hc^2 *1e6
    qmlat = Qmgev  * t0sqrt_ph^2 / t0ens[k] /hc^2 *1e6

    E, Z, Z_impr = get_spectr_data(path_spec, ens)

    g33_ll_s1_rec_aux = reconstr_corr(ens, E, Z, Z_impr, impr_set="1", nmax=Dict("E250"=>4,"D200"=>2)[ens.id], IMPR=true, RENORM=true, total=false, std=true)

    uwerr.(g33_ll_s1_rec_aux[end]); uwerr.(g33_ll_s1[k])

    trec_s1 = findfirst_uninterrupted(err.(g33_ll_s1_rec_aux[end]) .< err.(g33_ll_s1[k]))
    println("set1: ", trec_s1)

    for (j,q) in enumerate(Qlat)
        if j != 1
            #continue
        end 
        x0 = Float64.(collect(0: length(g33_lc_s1[k])-1))
        kk = KRNLsub.(x0, q, qmlat )

        
        _, data_orig = tmr_integrand(g33_ll_s1[k], q, qmlat, KRNLsub, pl=false, t0ens=t0ens[k], data=true)
        uwerr.(data_orig)

        fig = figure(figsize=(10,6))
        ttime = collect(0:length(data_orig)-1) .* value(a(ens.beta))
        errorbar(ttime, value.(data_orig), err.(data_orig), fmt="s", color="black", mfc="none", capsize=2, label=L"$\mathrm{LMA}$", lw=2, ms=6)

        colors = ["#E69F00", "#009E73", "#D55E00","#0072B2" ]
        for npi in eachindex(g33_ll_s1_rec_aux)
            _, data_rec = tmr_integrand(g33_ll_s1_rec_aux[npi], q, qmlat, KRNLsub, pl=false, t0ens=t0ens[k], data=true)
            uwerr.(data_rec)
            errorbar(ttime, value.(data_rec), err.(data_rec), fmt="s", mfc="none", capsize=2, color=colors[npi], label=L"$N_{\pi\pi}=\ $"*"$(npi)", lw=2, ms=6)

        end
        axvline(trec_s1 * value(a(ens.beta)), ls="dashed", color="gray")

        legend()
        xlim(0,3.5)
        ylabel(L"$\bar\Pi^{(3,3)}(Q^2/4)$")
        xlabel(L"$t \ [\mathrm{fm}]$")
        tight_layout()
        display(gcf())
        savefig(joinpath(path_plot, ens.id, "tmr_recostruction_mom$(Qgev[j]*4).pdf"))
        close("all")
    end
end