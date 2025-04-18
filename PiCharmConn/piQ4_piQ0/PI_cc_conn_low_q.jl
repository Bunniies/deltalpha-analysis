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


path_corr = "/Users/alessandroconigli/MyDrive/postdoc-mainz/projects/deltalpha/corr/impr_deriv/"
path_bdio_obs = "/Users/alessandroconigli/MyDrive/postdoc-mainz/projects/deltalpha/data"
path_store_pi = "/Users/alessandroconigli/MyDrive/postdoc-mainz/projects/deltalpha/PIdata/impr_deriv/low_q_kernel/scale_error_artificial/tmp/"


#======= PHYSICAL CONSTANTS ====================#
# const Qgev = [3., 5., 9.] # Q^2
const Qgev = [0.05, 0.1, 0.4, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0] ./ 4 # Q^2
const Qmgev = 9.0 # Qm^2

const KRNLsub = krnl_dα_sub # subtracted kernel Q4 - Q0

#============== READ CORRELATORS FROM BDIO FILES =================#

# enslist = sort([ "H101", "H102", "N101", "C101",
        #  "B450", "D450", "D452",
        #  "N202", "N203", "N200", "D200",  "E250",
        #  "N300", "J303", "E300",
        #  "J500"
# ])
enslist = ["E300", "D450"]
ensinfo = EnsInfo.(enslist)

path_ens = vcat([filter(x-> occursin(enslist[k], basename(x)), readdir(path_corr, join=true)) for k in eachindex(enslist)]...)

path_s1 = filter(x-> occursin("set1", x), path_ens)
path_s2 = filter(x-> occursin("set2", x), path_ens)

gcc_ll_s1 = Vector{Vector{uwreal}}(undef, length(enslist))
gcc_lc_s1 = Vector{Vector{uwreal}}(undef, length(enslist))
gcc_ll_s2 = Vector{Vector{uwreal}}(undef, length(enslist))
gcc_lc_s2 = Vector{Vector{uwreal}}(undef, length(enslist))


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
    gcc_ll_s1[k] = res_s1["gcc_ll_conn"]
    gcc_lc_s1[k] = res_s1["gcc_lc_conn"]
    
    fbs2 = BDIO_open(path_s2[k], "r")
    res_s2 = Dict()
    while ALPHAdobs_next_p(fbs2)
        d = ALPHAdobs_read_parameters(fbs2)
        sz = tuple(d["size"]...)
        ks = collect(d["keys"])
        res_s2 = ALPHAdobs_read_next(fbs2, size=sz, keys=ks)
    end
    BDIO_close!(fbs2)
    gcc_ll_s2[k] = res_s2["gcc_ll_conn"]
    gcc_lc_s2[k] = res_s2["gcc_lc_conn"]
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
# set 1 improvement coefficients
pi_cc_ll_s1 = [Vector{uwreal}(undef, length(Qgev)) for k in eachindex(ensinfo)] 
pi_cc_lc_s1 = [Vector{uwreal}(undef, length(Qgev)) for k in eachindex(ensinfo)]

# set 2 improvement coefficients
pi_cc_ll_s2 = [Vector{uwreal}(undef, length(Qgev)) for k in eachindex(ensinfo)]
pi_cc_lc_s2 = [Vector{uwreal}(undef, length(Qgev)) for k in eachindex(ensinfo)]

for (k, ens) in enumerate(ensinfo)
    println("Ensemble: ", ens.id)

    # println("using ensemble by ensemble t0/a^2")
    # Qlat  = Qgev .* t0sqrt_ph.^2 ./ t0ens[k] ./ hc^2 * 1e6
    # qmlat = Qmgev * t0sqrt_ph^2 / t0ens[k] / hc^2 * 1e6
    println("using t0/a^2 at symm point from Regensburg" )
    Qlat  = Qgev .* t0sqrt_ph.^2 ./ t0(ens.beta) ./ hc^2 * 1e6
    qmlat = Qmgev * t0sqrt_ph^2 / t0(ens.beta) / hc^2 * 1e6


    for (j,q) in enumerate(Qlat)
        if j == 12 
            pl_bool=true
        else
            pl_bool=false
        end
        pi_cc_ll_s1[k][j] =  tmr_integrand(gcc_ll_s1[k], q, qmlat, KRNLsub, pl=pl_bool, t0ens=t0ens[k])
        pi_cc_lc_s1[k][j] =  tmr_integrand(gcc_lc_s1[k], q, qmlat, KRNLsub, pl=false, t0ens=t0ens[k])
        pi_cc_ll_s2[k][j] =  tmr_integrand(gcc_ll_s2[k], q, qmlat, KRNLsub, pl=false, t0ens=t0ens[k])
        pi_cc_lc_s2[k][j] =  tmr_integrand(gcc_lc_s2[k], q, qmlat, KRNLsub, pl=false, t0ens=t0ens[k])
    end
end


## 
#============== SAVE results to BDIO ===============#

io = IOBuffer()
write(io, "PI charm-charm connected. ")
fb = ALPHAdobs_create(joinpath(path_store_pi, "PIcc_conn_newEnsStat.bdio"), io)

for (k, ens) in enumerate(ensinfo)
    extra = Dict{String, Any}("Ens" => ens.id)
    data = Dict{String, Array{uwreal}}(
        "picc_ll_s1" => pi_cc_ll_s1[k],
        "picc_lc_s1" => pi_cc_lc_s1[k],
        "picc_ll_s2" => pi_cc_ll_s2[k],
        "picc_lc_s2" => pi_cc_lc_s2[k]
    )

    ALPHAdobs_write(fb, data, extra=extra)
end
ALPHAdobs_close(fb)

##
#========= TEST READING =========#
fb = BDIO_open(joinpath(path_store_pi, "PIcc_conn.bdio"), "r")
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

