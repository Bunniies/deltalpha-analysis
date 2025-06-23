# Compute the Δ3388(-Q^2/4) - Δ3388(0) contribution with a given kernel and values of the virtualities 
# and store the results in BDIO for each ensemble analysed.
# FVC are included
# this code is speficically suited for subtracted kernel computation 
using Revise, HVPobs
using PyPlot, LaTeXStrings
using BDIO, ADerrors, ALPHAio
import ADerrors: err

#================ SET UP PLOT PARAMETERS ==================#

rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
rcParams["text.usetex"] =  true
rcParams["mathtext.fontset"]  = "cm"
rcParams["font.size"] =13
rcParams["axes.labelsize"] =22
rcParams["axes.titlesize"] = 18
plt.rc("text", usetex=true) # set to true if a LaTeX installation is present


include("../../utils/const.jl")
include("../../utils/types.jl")
include("../../utils/plot_utils.jl")
include("../../utils/IO_BDIO.jl")
include("../../utils/tools.jl")

const path_corr = "/Users/alessandroconigli/MyDrive/postdoc-mainz/projects/deltalpha/corr/impr_deriv/"
const path_bdio_obs = "/Users/alessandroconigli/MyDrive/postdoc-mainz/projects/deltalpha/data"
const path_store_pi = "/Users/alessandroconigli/MyDrive/postdoc-mainz/projects/deltalpha/PIdata/impr_deriv/new_strategy/Q_LD/"
const path_plot = "/Users/alessandroconigli/MyDrive/postdoc-mainz/projects/deltalpha/plots/new_strategy/ensembles"


#======= PHYSICAL CONSTANTS ====================#
const Qgev = [4, 5, 6, 7, 8, 9, 12] ./ 16 # Q^2 / 16 # additional very high values
const Qmgev = 9.0 # Qm^2

const KRNLsub = krnl_dα # non-subtracted kernel Q4 - Q0

#============== READ CORRELATORS FROM BDIO FILES =================#
@info("Reading Correlators")

enslist = sort([#"H101", "H102", "N101", "C101", "C102",  "D150"]) # C102 non computed
        #"B450", "N451", "N452", "D450", "D451", "D452"]) # D451 not computed 
        #"N202", "N203", "N200", "D251", "D200", "D201", "E250"]) # D201 enot computed
        "J307", "J303", "J304", "E300", # J304 not computed
         "J500", "J501"])

ensinfo = EnsInfo.(enslist)

path_ens = vcat([filter(x-> occursin(enslist[k], basename(x)), readdir(path_corr, join=true)) for k in eachindex(enslist)]...)

path_s1 = filter(x-> occursin("set1", x), path_ens)
path_s2 = filter(x-> occursin("set2", x), path_ens)

g3388_dlt_ll_s1 = Vector{Vector{uwreal}}(undef, length(enslist))
g3388_dlt_lc_s1 = Vector{Vector{uwreal}}(undef, length(enslist))
g3388_dlt_ll_s2 = Vector{Vector{uwreal}}(undef, length(enslist))
g3388_dlt_lc_s2 = Vector{Vector{uwreal}}(undef, length(enslist))

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
    g3388_dlt_ll_s1[k] = res_s1["g88_ll"]
    g3388_dlt_lc_s1[k] = res_s1["g88_lc"]
    
    fbs2 = BDIO_open(path_s2[k], "r")
    res_s2 = Dict()
    while ALPHAdobs_next_p(fbs2)
        d = ALPHAdobs_read_parameters(fbs2)
        sz = tuple(d["size"]...)
        ks = collect(d["keys"])
        res_s2 = ALPHAdobs_read_next(fbs2, size=sz, keys=ks)
    end
    BDIO_close!(fbs2)
    g3388_dlt_ll_s2[k] = res_s2["g88_ll"]
    g3388_dlt_lc_s2[k] = res_s2["g88_lc"]

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
#=========== COMPUTE TMR ==============#
@info("Computing TMR 88")
# set 1 improvement coefficients
pi_3388_ll_s1 = [Vector{uwreal}(undef, length(Qgev)) for k in eachindex(ensinfo)] 
pi_3388_lc_s1 = [Vector{uwreal}(undef, length(Qgev)) for k in eachindex(ensinfo)]
# set 2 improvement coefficients
pi_3388_ll_s2 = [Vector{uwreal}(undef, length(Qgev)) for k in eachindex(ensinfo)]
pi_3388_lc_s2 = [Vector{uwreal}(undef, length(Qgev)) for k in eachindex(ensinfo)]

for (k, ens) in enumerate(ensinfo)
    println("Ensemble: ", ens.id)

    Qlat  = Qgev .* t0sqrt_ph.^2 ./ t0ens[k] ./ hc^2 * 1e6
    qmlat = Qmgev * t0sqrt_ph^2 / t0ens[k] / hc^2 * 1e6

    for (j,q) in enumerate(Qlat)

        x0 = Float64.(collect(0: length(g3388_dlt_lc_s1[k])-1))
        kk = KRNLsub.(x0, q)
        if j == 7
            ss = joinpath(path_plot, ens.id)
            qqq = "12"
            pplot = true
        else
            ss = nothing
            qqq = ""
            pplot = false
        end

        pi_3388_ll_s1[k][j] = boundingMethod(g3388_dlt_ll_s1[k], ens, kk, "88", path_pl=ss, qval=qqq)
        pi_3388_lc_s1[k][j] = boundingMethod(g3388_dlt_lc_s1[k], ens, kk, "88", pl=false)
        pi_3388_ll_s2[k][j] = boundingMethod(g3388_dlt_ll_s2[k], ens, kk, "88", pl=false)
        pi_3388_lc_s2[k][j] = boundingMethod(g3388_dlt_lc_s2[k], ens, kk, "88", pl=false) 

        # pi_3388_ll_s1[k][j] =  tmr_integrand(g3388_dlt_ll_s1[k], q, KRNLsub, pl=pplot, t0ens=t0ens[k])
        # pi_3388_lc_s1[k][j] =  tmr_integrand(g3388_dlt_lc_s1[k], q, KRNLsub, pl=false, t0ens=t0ens[k])
        # pi_3388_ll_s2[k][j] =  tmr_integrand(g3388_dlt_ll_s2[k], q, KRNLsub, pl=false, t0ens=t0ens[k])
        # pi_3388_lc_s2[k][j] =  tmr_integrand(g3388_dlt_lc_s2[k], q, KRNLsub, pl=false, t0ens=t0ens[k])
    end
end


##
#============== SAVE results to BDIO ===============#
@info("Saving PI (88) results in BDIO")
io = IOBuffer()
write(io, "PI delta 88. ")
fb = ALPHAdobs_create(joinpath(path_store_pi, "tmp", "PI_88_batch5.bdio"), io)

for (k, ens) in enumerate(ensinfo)
    extra = Dict{String, Any}("Ens" => ens.id)
    data = Dict{String, Array{uwreal}}(
        "pi88_ll_s1" => pi_3388_ll_s1[k],
        "pi88_lc_s1" => pi_3388_lc_s1[k],
        "pi88_ll_s2" => pi_3388_ll_s2[k],
        "pi88_lc_s2" => pi_3388_lc_s2[k]
    )

    ALPHAdobs_write(fb, data, extra=extra)
end
ALPHAdobs_close(fb)
println("# Saving complete!")
##
#========= TEST READING =========#
fb = BDIO_open(joinpath(path_store_pi, "PI_88_beta1.bdio"), "r")
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


## plot comulative sum

fig = figure()
idx=8
cum = cumsum(pi_33_lc_s2_data[1][idx])
xx = collect(1:length(cum)) .* value(a(3.55))
uwerr.(cum)
uwerr.(pi_3388_lc_s2[1][idx])
errorbar(xx, value.(cum), err.(cum), fmt="s", ms=10, mfc="none", color="forestgreen")   
fill_between(xx, value(pi_3388_lc_s2[1][idx])-err(pi_3388_lc_s2[1][idx]), value(pi_3388_lc_s2[1][idx])+err(pi_3388_lc_s2[1][idx]), color="red",alpha=0.5)

axvline(35*value(a(3.55)), ls="dashed", color="gray")
xlabel(L"$t \ [\mathrm{fm}]$")
ylabel(L"$\bar{\Pi}^{(8,8)}_{\mathrm{ILD}}(Q^2/4)$")
title(L"$Q^2=$" * string(Qgev[idx]*4) * L"$\ \mathrm{GeV}^2$")
display(gcf())
tight_layout()
savefig(joinpath(path_plot,"E250/Pi88_cum_sum_E250_Qgev$(Qgev[idx]*4).pdf"))
close("all")