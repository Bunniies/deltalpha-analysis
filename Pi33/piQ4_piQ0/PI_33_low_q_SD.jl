# Compute the PI33(-Q^2/4) - PI33(0) contribution with a given kernel and values of the virtualities 
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
rcParams["font.size"] = 20
rcParams["axes.labelsize"] = 26
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
const path_store_pi = "/Users/alessandroconigli/MyDrive/postdoc-mainz/projects/deltalpha/PIdata/impr_deriv/low_q_kernel/scale_error_artificial/tmp/"
const path_fvc  = "/Users/alessandroconigli/Lattice/data/HVP/FSE"
const path_plot = "/Users/alessandroconigli/MyDrive/postdoc-mainz/projects/deltalpha/plots/ensembles"

#======= PHYSICAL CONSTANTS ====================#
const Qgev = [0.05, 0.1, 0.4, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0] ./ 4 # Q^2
# const Qgev = [1.0, 3.0, 5.0, 9.0] # Q^2
const Qmgev = 9.0 # Qm^2

const KRNLsub = krnl_dα # non-subtracted kernel
const WindSD = Window("SD")


enslist = sort([ #"H101", "H102", "N101", "C101", "C102", "D150",
          #"B450", "N451", "D450", "D451", "D452"])
          #"N202", "N203", "N200", "D251", "D200", "D201", "E250"])
           "J307", "J306", "J303", "J304", "E300", "F300",
          "J500", "J501"])

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
        krnl = KRNLsub.(T,q)
        pi_fvc_pi[k][j]  = abs(sum(fvc_pi .* krnl .* WindSD(T .* value(a(ens.beta)))))
        pi_fvc_k[k][j]   = abs(sum(fvc_k  .* krnl .* WindSD(T .* value(a(ens.beta)))))
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


for (k, ens) in enumerate(ensinfo)
    println("Ensemble: ", ens.id)

    Qlat = Qgev .* t0sqrt_ph.^2 ./ t0ens[k] ./hc^2 *1e6
    qmlat = Qmgev  * t0sqrt_ph^2 / t0ens[k] /hc^2 *1e6

    # Qlat = Qgev .* 0.1467.^2 ./ t0ens[k] ./hc^2 *1e6
    # qmlat = Qmgev  * 0.1467^2 / t0ens[k] /hc^2 *1e6

    # Qlat  = Qgev .* value.(t0sqrt_ph.^2) ./ t0ens[k] ./ hc^2 * 1e6
    # qmlat = Qmgev * value(t0sqrt_ph^2) / t0ens[k] / hc^2 * 1e6


    for (j,q) in enumerate(Qlat)
        x0 = Float64.(collect(0: length(g33_lc_s1[k])-1))
        kk = KRNLsub.(x0, q)
        
        if j == 12
            ss = joinpath(path_plot, ens.id)
            qqq = "9"
        else
            ss = nothing
            qqq = ""
        end
        # pi_33_ll_s1[k][j] = boundingMethod(g33_ll_s1[k], ens, kk, "33", path_pl=ss, qval=qqq)
        # pi_33_lc_s1[k][j] = boundingMethod(g33_lc_s1[k], ens, kk, "33")
        # pi_33_ll_s2[k][j] = boundingMethod(g33_ll_s2[k], ens, kk, "33")
        # pi_33_lc_s2[k][j] = boundingMethod(g33_lc_s2[k], ens, kk, "33")

        pi_33_ll_s1[k][j] = tmr_integrand(g33_ll_s1[k], q, KRNLsub, pl=true,  t0ens=t0ens[k], wind=WindSD)
        pi_33_lc_s1[k][j] = tmr_integrand(g33_lc_s1[k], q, KRNLsub, pl=false, t0ens=t0ens[k], wind=WindSD)
        pi_33_ll_s2[k][j] = tmr_integrand(g33_ll_s2[k], q, KRNLsub, pl=false, t0ens=t0ens[k], wind=WindSD)
        pi_33_lc_s2[k][j] = tmr_integrand(g33_lc_s2[k], q, KRNLsub, pl=false, t0ens=t0ens[k], wind=WindSD)
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
write(io, "PI delta 33 SD low q. ")
fb = ALPHAdobs_create(joinpath(path_store_pi, "PI_33_SD_beta45.bdio"), io)

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
fb = BDIO_open(joinpath(path_store_pi, "multi_mom", "PI_33.bdio"), "r")
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