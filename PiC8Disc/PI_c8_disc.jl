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


include("../utils/const.jl")
include("../utils/types.jl")
include("../utils/plot_utils.jl")
include("../utils/IO_BDIO.jl")
include("../utils/tools.jl")

path_corr = "/Users/alessandroconigli/MyDrive/postdoc-mainz/projects/deltalpha/corr/impr_deriv/"
path_bdio_obs = "/Users/alessandroconigli/MyDrive/postdoc-mainz/projects/deltalpha/data"
path_store_pi = "/Users/alessandroconigli/MyDrive/postdoc-mainz/projects/deltalpha/PIdata/impr_deriv/"

#======= PHYSICAL CONSTANTS ====================#
const Qgev = [3., 5., 9.] # Q^2
const Qmgev = 9.0 # Qm^2

const KRNL = krnl_dα_qhalf # non-subtracted kernel

#============== READ CORRELATORS FROM BDIO FILES =================#

enslist = sort([  "H102", "N101", "C101", "C102", "D150",
        "N451", "D450", "D451", "D452",
        "N203", "N200", "D200", "D201", "E250",
        "J303", "E300",
        "J501"])


ensinfo = EnsInfo.(enslist)

path_ens = vcat([filter(x-> occursin(enslist[k], basename(x)), readdir(path_corr, join=true)) for k in eachindex(enslist)]...)

path_s1 = filter(x-> occursin("set1", x), path_ens)
path_s2 = filter(x-> occursin("set2", x), path_ens)

gc8_cc_s1 = Vector{Vector{uwreal}}(undef, length(enslist))
gc8_cc_s2 = Vector{Vector{uwreal}}(undef, length(enslist))

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
    gc8_cc_s1[k] = res_s1["gc8_cc_disc"]
    
    fbs2 = BDIO_open(path_s2[k], "r")
    res_s2 = Dict()
    while ALPHAdobs_next_p(fbs2)
        d = ALPHAdobs_read_parameters(fbs2)
        sz = tuple(d["size"]...)
        ks = collect(d["keys"])
        res_s2 = ALPHAdobs_read_next(fbs2, size=sz, keys=ks)
    end
    BDIO_close!(fbs2)
    gc8_cc_s2[k] = res_s2["gc8_cc_disc"]
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
pi_c8_cc_s1 = [Vector{uwreal}(undef, length(Qgev)) for k in eachindex(ensinfo)]

# set 2 improvement coefficients
pi_c8_cc_s2 = [Vector{uwreal}(undef, length(Qgev)) for k in eachindex(ensinfo)]

for (k, ens) in enumerate(ensinfo)
    println("Ensemble: ", ens.id)

    Qlat  = Qgev .* value.(t0sqrt_ph.^2) ./ t0ens[k] ./ hc^2 * 1e6
    qmlat = Qmgev * value(t0sqrt_ph^2) / t0ens[k] / hc^2 * 1e6

    for (j,q) in enumerate(Qlat)
        pi_c8_cc_s1[k][j] =  tmr_integrand(gc8_cc_s1[k], q, KRNL, pl=true, t0ens=t0ens[k], wind=Window("SD"))
        pi_c8_cc_s2[k][j] =  tmr_integrand(gc8_cc_s2[k], q, KRNL, pl=false, t0ens=t0ens[k])
    end
end

##
## 
#============== SAVE results to BDIO ===============#
# popat!(pi_c8_cc_s1, 6 )
# popat!(pi_c8_cc_s2, 6 )
# popat!(ensinfo, 6)
io = IOBuffer()
write(io, "PI charm-8 disconnected. ")
fb = ALPHAdobs_create(joinpath(path_store_pi, "PIc8_disc.bdio"), io)

for (k, ens) in enumerate(ensinfo)
    extra = Dict{String, Any}("Ens" => ens.id)
    data = Dict{String, Array{uwreal}}(
        "pic8_cc_s1" => pi_c8_cc_s1[k],
        "pic8_cc_s2" => pi_c8_cc_s2[k]
    )

    ALPHAdobs_write(fb, data, extra=extra)
end
ALPHAdobs_close(fb)

##
#========= TEST READING =========#
fb = BDIO_open(joinpath(path_store_pi, "PIc8_disc.bdio"), "r")
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

