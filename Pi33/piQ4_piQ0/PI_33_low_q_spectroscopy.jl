# Compute the PI33(-Q^2/4) - PI33(0) contribution with a given kernel and values of the virtualities 
# and store the results in BDIO for each ensemble analysed.
# FVC are included
# this code is speficically suited for subtracted kernel computation 
using Revise, HVPobs
using HDF5
using PyPlot, LaTeXStrings
using BDIO, ADerrors, ALPHAio
import ADerrors: err
using Statistics

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
# const KRNLsub = krnl_dα # subtracted kernel


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

##
#============== READ t0 FROM BDIO FILES =================#
dir_path = filter(x-> basename(x) in enslist, readdir(path_bdio_obs, join=true))
spectrum_path = vcat(filter(!isempty, [filter(x-> occursin("spectrum.bdio", x)  , readdir(dir_path[k], join=true)) for k in eachindex(dir_path)])...)

t0ens = Vector{uwreal}(undef, length(enslist))

for (k, ens) in enumerate(ensinfo)
    t0ens_aux = read_BDIO(spectrum_path[k], "spectrum", "t0")[1]
    t0ens[k] = t0ens_aux
end


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

        pi_33_ll_s1[k][j] = tmr_integrand(g33_ll_s1[k], q, qmlat, KRNLsub, pl=false, t0ens=t0ens[k])
        pi_33_lc_s1[k][j] = tmr_integrand(g33_lc_s1[k], q, qmlat, KRNLsub, pl=false, t0ens=t0ens[k])
        pi_33_ll_s2[k][j] = tmr_integrand(g33_ll_s2[k], q, qmlat, KRNLsub, pl=false, t0ens=t0ens[k])
        pi_33_lc_s2[k][j] = tmr_integrand(g33_lc_s2[k], q, qmlat, KRNLsub, pl=false, t0ens=t0ens[k])

        pi_33_ll_s1_rec[k][j] = boundingMethod(g33_ll_s1_rec[k], ens, kk, "33", pl=true)
        pi_33_lc_s1_rec[k][j] = boundingMethod(g33_lc_s1_rec[k], ens, kk, "33", pl=false)
        pi_33_ll_s2_rec[k][j] = boundingMethod(g33_ll_s2_rec[k], ens, kk, "33", pl=true)
        pi_33_lc_s2_rec[k][j] = boundingMethod(g33_lc_s2_rec[k], ens, kk, "33", pl=false)

    end
end
##
#============ ADD FVC ============#
@info("Adding Finite-Volume Corrections")
for k in eachindex(ensinfo)

    if ensinfo[k].kappa_l == ensinfo[k].kappa_s
        pi_33_ll_s1_rec[k] .+= 1.5 * pi_fvc_pi[k]
        pi_33_lc_s1_rec[k] .+= 1.5 * pi_fvc_pi[k]
        pi_33_ll_s2_rec[k] .+= 1.5 * pi_fvc_pi[k]
        pi_33_lc_s2_rec[k] .+= 1.5 * pi_fvc_pi[k]   
    else 
        pi_33_ll_s1_rec[k] .+= pi_fvc_pi[k] .+ pi_fvc_k[k]
        pi_33_lc_s1_rec[k] .+= pi_fvc_pi[k] .+ pi_fvc_k[k]
        pi_33_ll_s2_rec[k] .+= pi_fvc_pi[k] .+ pi_fvc_k[k]
        pi_33_lc_s2_rec[k] .+= pi_fvc_pi[k] .+ pi_fvc_k[k]
    end

end

##
#============== SAVE results to BDIO ===============#
@info("Saving PI 33 results in BDIO")
io = IOBuffer()
write(io, "PI delta 33 low q. ")
fb = ALPHAdobs_create(joinpath(path_store_pi, "PI_33_spectroscopy.bdio"), io)

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

## check results within sigma
uwerr.(pi_33_ll_s1[2]); uwerr.(pi_33_lc_s1[2]); uwerr.(pi_33_ll_s2[2]); uwerr.(pi_33_lc_s2[2])
uwerr.(pi_33_ll_s1_rec[2]); uwerr.(pi_33_lc_s1_rec[2]); uwerr.(pi_33_ll_s2_rec[2]); uwerr.(pi_33_lc_s2_rec[2])

(pi_33_ll_s1_rec[2] .- pi_33_ll_s1[2]) ./ err.(pi_33_ll_s1[2])
(pi_33_ll_s2_rec[2] .- pi_33_ll_s2[2]) ./ err.(pi_33_ll_s2[2])

(pi_33_lc_s1_rec[2] .- pi_33_lc_s1[2]) ./ err.(pi_33_lc_s1[2])
(pi_33_lc_s2_rec[2] .- pi_33_lc_s2[2]) ./ err.(pi_33_lc_s2[2])