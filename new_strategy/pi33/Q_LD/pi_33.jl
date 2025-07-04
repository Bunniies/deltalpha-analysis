# Compute the PI33(-Q^2/16) - PI33(0) contribution with a given kernel and values of the virtualities 
# and store the results in BDIO for each ensemble analysed.
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
const path_store_pi = "/Users/alessandroconigli/MyDrive/postdoc-mainz/projects/deltalpha/PIdata/impr_deriv/new_strategy/Q_LD/tmp"
const path_plot = "/Users/alessandroconigli/MyDrive/postdoc-mainz/projects/deltalpha/plots/new_strategy/ensembles"

#======= PHYSICAL CONSTANTS ====================#
const Qgev = [4, 5, 6, 7, 8, 9, 12] ./ 16 # Q^2/16 # 
# const Qgev = [1.0] # Q^2
const Qmgev = 9.0 # Qm^2

const KRNLsub = krnl_dα_sub # subtracted kernel
# const KRNLsub = krnl_dα # non-subtracted kernel

enslist = sort([ #"H101", "H102", "N101", "C101", "C102", "D150"]
        #   "B450", "N451", "N452", "D450", "D451", "D452"]
        # "N202", "N203", "N200", "D251", "D201"] # D200, E250 removed, use spectroscopy instead
         "J307", "J306", "J303", "J304", "E300","F300"]
        #"J500", "J501"]
)

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
#=========== COMPUTE TMR ==============#
# set 1 improvement coefficients
pi_33_ll_s1 = [Vector{uwreal}(undef, length(Qgev)) for k in eachindex(ensinfo)] 
pi_33_lc_s1 = [Vector{uwreal}(undef, length(Qgev)) for k in eachindex(ensinfo)]

# set 2 improvement coefficients
pi_33_ll_s2  = [Vector{uwreal}(undef, length(Qgev)) for k in eachindex(ensinfo)] 
pi_33_lc_s2  = [Vector{uwreal}(undef, length(Qgev)) for k in eachindex(ensinfo)]


for (k, ens) in enumerate(ensinfo)
    println("Ensemble: ", ens.id)

    Qlat = Qgev .* t0sqrt_ph.^2 ./ t0ens[k] ./hc^2 *1e6
    qmlat = Qmgev  * t0sqrt_ph^2 / t0ens[k] /hc^2 *1e6

    for (j,q) in enumerate(Qlat)
        x0 = Float64.(collect(0: length(g33_lc_s1[k])-1))
        kk = KRNLsub.(x0, q, qmlat)
        # kk = KRNLsub.(x0, q)
        
        if j == 12
            ss = joinpath(path_plot, ens.id)
            qqq = string(q)
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
    end 
end

##
#============== SAVE results to BDIO ===============#
@info("Saving PI 33 results in BDIO")
io = IOBuffer()
write(io, "PI 33 Q_LD ")
fb = ALPHAdobs_create(joinpath(path_store_pi, "PI_33_batch4.bdio"), io)

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
fb = BDIO_open(joinpath(path_store_pi, "PI_33_no_BM.bdio"), "r")
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
