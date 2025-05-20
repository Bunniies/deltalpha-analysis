# Here we compute the FVE using the Hansen Patella method for the various Q^2 used throughout the analysis
# the FVE are computed with respect to the infinite volume and not a given value of L_{ref}
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

const path_fvc  = "/Users/alessandroconigli/Lattice/data/HVP/FSE"
const path_bdio_obs = "/Users/alessandroconigli/MyDrive/postdoc-mainz/projects/deltalpha/data"
const path_store_pi = "/Users/alessandroconigli/MyDrive/postdoc-mainz/projects/deltalpha/PIdata/impr_deriv/new_strategy/Q_ID"

#======= PHYSICAL CONSTANTS ====================#
const Qgev = [9, 12, 15, 18] ./4  # Q^2 # additional very high values


enslist = sort([ "H101", "H102", "N101", "C101", "C102", "D150",
          "B450", "N451", "N452", "D450", "D451", "D452",
         "N202", "N203", "N200", "D251", "D200", "D201", "E250",
          "J307", "J306", "J303", "J304", "E300", "F300",
         "J500", "J501"]
)

ensinfo = EnsInfo.(enslist)
const Qmgev = 9.0 # Qm^2

const KRNLsub = krnl_dα_qhalf_sub # subtracted kernel
const WindSD = Window("SD")
const WindILD = Window("ILD")

@warn("Always check the kernel that you are using!! \\ It has to match with the one used for the contribution you are considering")
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
pi_fvc_pi_SD = [Vector{uwreal}(undef, length(Qgev)) for k in eachindex(ensinfo)]
pi_fvc_k_SD  = [Vector{uwreal}(undef, length(Qgev)) for k in eachindex(ensinfo)]

pi_fvc_pi_ILD = [Vector{uwreal}(undef, length(Qgev)) for k in eachindex(ensinfo)]
pi_fvc_k_ILD  = [Vector{uwreal}(undef, length(Qgev)) for k in eachindex(ensinfo)]


for (k,ens) in enumerate(ensinfo)
    println("Ensemble: $(ens.id)")

    Qlat = Qgev .* t0sqrt_ph.^2 ./ t0ens[k] ./hc^2 *1e6
    qmlat = Qmgev  * t0sqrt_ph^2 / t0ens[k] /hc^2 *1e6

    # Qlat = Qgev .* value.(t0sqrt_ph.^2 ./ t0ens[k]) ./hc^2 *1e6
    # qmlat = Qmgev * value(t0sqrt_ph^2 / t0ens[k]) /hc^2 *1e6

    fvc_raw_pi = get_fvc(joinpath(path_fvc, "JKMPI" ), ens.id)
    fvc_pi = vcat(sum(fvc_raw_pi, dims=1)...)
    
    fvc_raw_k = get_fvc(joinpath(path_fvc, "JKMK" ), ens.id)
    fvc_k = vcat(sum(fvc_raw_k, dims=1)...)
    
    T = Float64.(collect(1:length(fvc_pi)))

    for (j,q) in enumerate(Qlat)
        krnl = KRNLsub.(T,q, qmlat)
        pi_fvc_pi_SD[k][j]  = abs(sum(fvc_pi .* krnl))
        pi_fvc_k_SD[k][j]   = abs(sum(fvc_k  .* krnl))

        # pi_fvc_pi_SD[k][j]  = abs(sum(fvc_pi .* krnl .* WindSD(T  .* value(a(ens.beta)))))
        # pi_fvc_k_SD[k][j]   = abs(sum(fvc_k  .* krnl .* WindSD(T  .* value(a(ens.beta)))))
        # pi_fvc_pi_ILD[k][j] = abs(sum(fvc_pi .* krnl .* WindILD(T .* value(a(ens.beta)))))
        # pi_fvc_k_ILD[k][j]  = abs(sum(fvc_k  .* krnl .* WindILD(T .* value(a(ens.beta)))))
    end
end

##
#============== SAVE results to BDIO ===============#

@info("Saving HP FVE resulst in BDIO")
io = IOBuffer()
write(io, "FVE HP")

fb = ALPHAdobs_create(joinpath(path_store_pi, "FVE_HP_Q_ID.bdio"), io)

for (k,ens) in enumerate(ensinfo)
    extra = Dict{String, Any}("Ens" => ens.id)
    data = Dict{String, Array{uwreal}}(
        "fvc_pi" => pi_fvc_pi_SD[k],
        # "fvc_pi_ILD" => pi_fvc_pi_ILD[k],
        "fvc_k" => pi_fvc_k_SD[k],
        # "fvc_k_ILD" => pi_fvc_k_ILD[k]
    )
    ALPHAdobs_write(fb, data, extra=extra)
end
ALPHAdobs_close(fb)
println("# Saving complete")

##
#========= TEST READING =========#
fb = BDIO_open(joinpath(path_store_pi, "FVE_HP_Q_LD.bdio"), "r")
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
