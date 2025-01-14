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
include("./func_comb.jl")

path_fvc  = "/Users/alessandroconigli/Lattice/data/HVP/FSE"
path_corr = "/Users/alessandroconigli/MyDrive/postdoc-mainz/projects/deltalpha/corr/impr_deriv/"
path_bdio_obs = "/Users/alessandroconigli/MyDrive/postdoc-mainz/projects/deltalpha/data"
path_store_pi = "/Users/alessandroconigli/MyDrive/postdoc-mainz/projects/deltalpha/PIdata/impr_deriv/high_q_kernel/scale_error_artificial"
##
#======= PHYSICAL CONSTANTS ====================#
# const Qgev = [3., 5., 9.] # Q^2
const Qgev = [0.05, 0.1, 0.4, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0] # Q^2

const Qmgev = 9.0 # Qm^2

const KRNL = krnl_dα_qhalf # non-subtracted kernel

#============== READ CORRELATORS FROM BDIO FILES =================#
@info("Reading Correlators")

enslist = sort([ "H102", "N101", "C101", "C102", "D150",
         "N451", "D450", "D451", "D452",
         "N203", "N200", "D251", "D200", "D201", "E250",
         "J303", "J306", "J304", "E300", "F300",
         "J501"
         ])
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
    g3388_dlt_ll_s1[k] = res_s1["g3388_dlt_ll"]
    g3388_dlt_lc_s1[k] = res_s1["g3388_dlt_lc"]
    
    fbs2 = BDIO_open(path_s2[k], "r")
    res_s2 = Dict()
    while ALPHAdobs_next_p(fbs2)
        d = ALPHAdobs_read_parameters(fbs2)
        sz = tuple(d["size"]...)
        ks = collect(d["keys"])
        res_s2 = ALPHAdobs_read_next(fbs2, size=sz, keys=ks)
    end
    BDIO_close!(fbs2)
    g3388_dlt_ll_s2[k] = res_s2["g3388_dlt_ll"]
    g3388_dlt_lc_s2[k] = res_s2["g3388_dlt_lc"]
end
#
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
    Qlat = Qgev .* t0sqrt_ph.^2 ./ t0ens[k] ./hc^2 *1e6

    fvc_raw_pi = get_fvc(joinpath(path_fvc, "JKMPI" ), ens.id)
    fvc_pi = vcat(sum(fvc_raw_pi, dims=1)...)
    
    fvc_raw_k = get_fvc(joinpath(path_fvc, "JKMK" ), ens.id)
    fvc_k = vcat(sum(fvc_raw_k, dims=1)...)
    
    T = Float64.(collect(1:length(fvc_pi)))

    for (j,q) in enumerate(Qlat)
        krnl = KRNL.(T,q)
        pi_fvc_pi[k][j] = abs(sum(fvc_pi .* krnl))
        pi_fvc_k[k][j]  = abs(sum(fvc_k .* krnl))
    end
end

##
#=========== COMPUTE TMR ==============#
@info("Computing TMR (33-88)")
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

    # Qlat  = Qgev .* value.(t0sqrt_ph.^2) ./ t0ens[k] ./ hc^2 * 1e6
    # qmlat = Qmgev * value(t0sqrt_ph^2) / t0ens[k] / hc^2 * 1e6
    for (j,q) in enumerate(Qlat)
        pi_3388_ll_s1[k][j] =  tmr_integrand(g3388_dlt_ll_s1[k], q, KRNL, pl=false, t0ens=t0ens[k])
        pi_3388_lc_s1[k][j] =  tmr_integrand(g3388_dlt_lc_s1[k], q, KRNL, pl=false, t0ens=t0ens[k])
        pi_3388_ll_s2[k][j] =  tmr_integrand(g3388_dlt_ll_s2[k], q, KRNL, pl=false, t0ens=t0ens[k])
        pi_3388_lc_s2[k][j] =  tmr_integrand(g3388_dlt_lc_s2[k], q, KRNL, pl=false, t0ens=t0ens[k])
    end
end

##
#============ ADD FVC ============#
@info("Adding Finite-Volume Corrections")
for k in eachindex(ensinfo)
    pifvc = pi_fvc_pi[k]
    kfvc = pi_fvc_k[k]

    pi_3388_ll_s1[k] .+= (pifvc + kfvc -2/3*kfvc)
    pi_3388_lc_s1[k] .+= (pifvc + kfvc -2/3*kfvc)
    pi_3388_ll_s2[k] .+= (pifvc + kfvc -2/3*kfvc)
    pi_3388_lc_s2[k] .+= (pifvc + kfvc -2/3*kfvc)
end
##
#============== SAVE results to BDIO ===============#
@info("Saving PI (33-88) results in BDIO")
io = IOBuffer()
write(io, "PI delta 33-88. ")
fb = ALPHAdobs_create(joinpath(path_store_pi, "PI_3388_dlt.bdio"), io)

for (k, ens) in enumerate(ensinfo)
    extra = Dict{String, Any}("Ens" => ens.id)
    data = Dict{String, Array{uwreal}}(
        "pi3388_ll_s1" => pi_3388_ll_s1[k],
        "pi3388_lc_s1" => pi_3388_lc_s1[k],
        "pi3388_ll_s2" => pi_3388_ll_s2[k],
        "pi3388_lc_s2" => pi_3388_lc_s2[k]
    )

    ALPHAdobs_write(fb, data, extra=extra)
end
ALPHAdobs_close(fb)
println("# Saving complete!")
##
#========= TEST READING =========#
fb = BDIO_open(joinpath(path_store_pi, "PI_3388_dlt.bdio"), "r")
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


## test derivatives
uwerr(t0sqrt_ph_noerr)
t0sqrt_ph_noerr.ids
ensembles(t0sqrt_ph_noerr)

uwerr(res["C101"]["pi3388_ll_s1"][1])
res["C101"]["pi3388_ll_s1"][1].ids
ensembles( res["C101"]["pi3388_ll_s1"][1])

derivative( res["C101"]["pi3388_ll_s1"][1], t0sqrt_ph_noerr)



aa = 5. * t0sqrt_ph_noerr
uwerr(aa)
aa.ids
derivative(aa, t0sqrt_ph_noerr)