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

const KRNL = krnl_dα # subtracted kernel
# const KRNLsub = krnl_dα # non-subtracted kernel

enslist = sort([ #"H102", "N101", "C101", "C102", "D150"
          #"N451", "N452", "D450", "D451", "D452",
          #"N203", "N200", "D251", "D200", "D201", "E250"
         "J303", "J304", "E300",
          "J501"
])

ensinfo = EnsInfo.(enslist)

path_ens = vcat([filter(x-> occursin(enslist[k], basename(x)), readdir(path_corr, join=true)) for k in eachindex(enslist)]...)

path_s1 = filter(x-> occursin("set1", x), path_ens)
path_s2 = filter(x-> occursin("set2", x), path_ens)

g08_ll_s1 = Vector{Vector{uwreal}}(undef, length(enslist))
g08_lc_s1 = Vector{Vector{uwreal}}(undef, length(enslist))
g08_ll_s2 = Vector{Vector{uwreal}}(undef, length(enslist))
g08_lc_s2 = Vector{Vector{uwreal}}(undef, length(enslist))

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
    # g08_ll_s1[k] = res_s1["g08_ll"]
    g08_lc_s1[k] = res_s1["g08_lc"]
    
    fbs2 = BDIO_open(path_s2[k], "r")
    res_s2 = Dict()
    while ALPHAdobs_next_p(fbs2)
        d = ALPHAdobs_read_parameters(fbs2)
        sz = tuple(d["size"]...)
        ks = collect(d["keys"])
        res_s2 = ALPHAdobs_read_next(fbs2, size=sz, keys=ks)
    end
    BDIO_close!(fbs2)
    # g08_ll_s2[k] = res_s2["g08_ll"]
    g08_lc_s2[k] = res_s2["g08_lc"]
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
@info("Computing TMR (08)")
# set 1 improvement coefficients
pi_08_ll_s1 = [Vector{uwreal}(undef, length(Qgev)) for k in eachindex(ensinfo)] 
pi_08_lc_s1 = [Vector{uwreal}(undef, length(Qgev)) for k in eachindex(ensinfo)]

# set 2 improvement coefficients
pi_08_ll_s2 = [Vector{uwreal}(undef, length(Qgev)) for k in eachindex(ensinfo)]
pi_08_lc_s2 = [Vector{uwreal}(undef, length(Qgev)) for k in eachindex(ensinfo)]

for (k, ens) in enumerate(ensinfo)
    println("Ensemble: ", ens.id)

    # Qlat  = Qgev .* value.(t0sqrt_ph.^2) ./ t0ens[k] ./ hc^2 * 1e6
    # qmlat = Qmgev * value(t0sqrt_ph^2) / t0ens[k] / hc^2 * 1e6

    Qlat  = Qgev .* t0sqrt_ph.^2 ./ t0ens[k] ./ hc^2 * 1e6
    qmlat = Qmgev * t0sqrt_ph^2 / t0ens[k] / hc^2 * 1e6

    for (j,q) in enumerate(Qlat)

        x0 = Float64.(collect(0: length(g08_lc_s1[k])-1))
        kk = KRNL.(x0, q)
        if j == 12
            ss = joinpath(path_plot, ens.id)
            qqq = "9"
        else
            ss = nothing
            qqq = ""
        end

        pi_08_lc_s1[k][j] = boundingMethod(g08_lc_s1[k], ens, kk, "08", path_pl=ss, qval=qqq)
        pi_08_lc_s2[k][j] = boundingMethod(g08_lc_s2[k], ens, kk, "08")

        # pi_08_ll_s1[k][j] =  tmr_integrand(g08_ll_s1[k], q, KRNL, pl=false, t0ens=t0ens[k])
        # pi_08_lc_s1[k][j] =  tmr_integrand(g08_lc_s1[k], q, KRNL, pl=false, t0ens=t0ens[k])
        # pi_08_ll_s2[k][j] =  tmr_integrand(g08_ll_s2[k], q, KRNL, pl=false, t0ens=t0ens[k])
        # pi_08_lc_s2[k][j] =  tmr_integrand(g08_lc_s2[k], q, KRNL, pl=false, t0ens=t0ens[k])
    end
end    

##
#============== SAVE results to BDIO ===============#
@info("Saving PI (08) results in BDIO")
io = IOBuffer()
write(io, "PI 08. ")
fb = ALPHAdobs_create(joinpath(path_store_pi, "PI_08_batch4.bdio"), io)

for (k, ens) in enumerate(ensinfo)
    extra = Dict{String, Any}("Ens" => ens.id)
    data = Dict{String, Array{uwreal}}(
        # "pi08_ll_s1" => pi_08_ll_s1[k],
        "pi08_lc_s1" => pi_08_lc_s1[k],
        # "pi08_ll_s2" => pi_08_ll_s2[k],
        "pi08_lc_s2" => pi_08_lc_s2[k]
    )

    ALPHAdobs_write(fb, data, extra=extra)
end
ALPHAdobs_close(fb)
println("# Saving complete!")