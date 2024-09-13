using Revise
using HVPobs
using PyPlot, LaTeXStrings
using OrderedCollections
using TimerOutputs
using BDIO, ADerrors
import ADerrors: err
using ALPHAio

include("../utils/const.jl")
include("../utils/types.jl")
include("../utils/tools.jl")
include("../utils/IO_BDIO.jl")

const KRNL= krnl_dα # non-subtracted kernel


const path_data = "/Users/alessandroconigli/Lattice/data/HVP/2ptdata"
const path_rw   = "/Users/alessandroconigli/Lattice/data/HVP/rwf_deflated"
const path_store_bdio = "/Users/alessandroconigli/MyDrive/postdoc-mainz/projects/deltalpha/corr/impr_deriv"
const path_fvc  = "/Users/alessandroconigli/Lattice/data/HVP/FSE"
const path_bdio_obs = "/Users/alessandroconigli/MyDrive/postdoc-mainz/projects/deltalpha/data"

const IMPR = true
const STD_DERIV = false
const RENORM = true

# enslist = sort([  "H102", "N101", "C101", "C102", "D150",
        #  "N451", "D450", "D451", "D452",
        #  "N203", "N200", "D200", "D201", "E250",
        #  "J303", "E300"])

enslist = sort([ "H102" ])
ensinfo = EnsInfo.(enslist)
Nens = length(ensinfo)

# conn
g88_ll_conn = Vector{Corr}(undef, Nens)
g88_lc_conn = Vector{Corr}(undef, Nens)

g08_ll_conn = Vector{Corr}(undef, Nens)
g08_lc_conn = Vector{Corr}(undef, Nens)

g00_ll_conn = Vector{Corr}(undef, Nens)
g00_lc_conn = Vector{Corr}(undef, Nens)

# disc
g88_ll_disc = Vector{Corr}(undef, Nens)
g88_lc_disc = Vector{Corr}(undef, Nens)

g08_ll_disc = Vector{Corr}(undef, Nens)
g08_lc_disc = Vector{Corr}(undef, Nens)

@time begin
    for (k,ens) in enumerate(ensinfo)
        gll_ll, gll_lc = corrConnected(path_data, ens, "light", path_rw=path_rw, impr=IMPR, impr_set="1", std=STD_DERIV)
        gss_ll, gss_lc = corrConnected(path_data, ens, "strange", path_rw=path_rw, impr=IMPR, impr_set="1", std=STD_DERIV)

        gll_ll_noimp, gll_lc_noimp = corrConnected(path_data, ens, "light", path_rw=path_rw, impr=false, impr_set="1", std=STD_DERIV)
        gss_ll_noimp, gss_lc_noimp = corrConnected(path_data, ens, "strange", path_rw=path_rw, impr=false, impr_set="1", std=STD_DERIV)
        
        g88_ll_conn[k]  = Corr(1/6 .*( gll_ll.obs + 2*gss_ll.obs), ens.id, "G88_ll_conn" )
        g88_lc_conn[k]  = Corr(1/6 .*( gll_lc.obs + 2*gss_lc.obs), ens.id, "G88_lc_conn" )

        g08_ll_conn[k] = Corr(1/(2*sqrt(3)) .* (gll_ll.obs .- gss_ll.obs), ens.id,  "G08_ll_conn")
        g08_lc_conn[k] = Corr(1/(2*sqrt(3)) .* (gll_lc.obs .- gss_lc.obs), ens.id,  "G08_lc_conn")

        g00_lc_conn[k] = Corr(1/4 * (2*gll_lc_noimp.obs .+ gss_lc_noimp.obs), ens.id, "G00_lc_conn")

        if ens.kappa_l != ens.kappa_s
            try
                g88_ll_disc[k], g88_lc_disc[k], _ = corrDisconnected(path_data, ens, "88", path_rw=path_rw, impr=IMPR, impr_set="1", std=STD_DERIV)
                g08_ll_disc[k], g08_lc_disc[k], _ = corrDisconnected(path_data, ens, "80", path_rw=path_rw, impr=IMPR, impr_set="1", std=STD_DERIV)
            catch
                T = HVPobs.Data.get_T(ens.id)
                g08_ll_disc[k] =  g08_lc_disc[k] = g88_lc_disc[k] = g88_ll_disc[k] = Corr(fill(uwreal(0.0), T), ens.id, "Gdisc")
            end
        else
            T = HVPobs.Data.get_T(ens.id)
            g08_ll_disc[k] = g08_lc_disc[k] =  g88_lc_disc[k] = g88_ll_disc[k] = Corr(fill(uwreal(0.0), T), ens.id, "Gdisc")
        end
        
        # if RENORM
            # Z8 = get_Z8(ens, impr_set="1")
            # Z08 = get_Z08(ens, impr_set="1")
            # renormalize!(g88_ll_conn[k], Z8^2)
            # renormalize!(g88_lc_conn[k], Z8)
            # 
            # renormalize!(g08_lc_conn[k], Z8)
            # renormalize!(g08_lc_disc[k], Z8)
            # 
            # renormalize!(g00_lc_conn[k], Z08)
            # 
        # end

        Z8 = get_Z8(ens, impr_set="1")
        Z08 = get_Z08(ens, impr_set="1")

        g88_ll_conn[k].obs[:] = Z8^2 .* g88_ll_conn[k].obs[:] .+ 2 .* Z8 * Z08 .* g08_ll_conn[k].obs[:] 
        g88_lc_conn[k].obs[:] = Z8 .* g88_lc_conn[k].obs[:] .+   Z08 .* g08_lc_conn[k].obs[:] 

        g88_ll_disc[k].obs[:] = Z8^2 .* g88_ll_disc[k].obs[:] .+ 2 .* Z8 * Z08 .* g08_ll_disc[k].obs[:] 
        g88_lc_disc[k].obs[:] = Z8 .* g88_lc_disc[k].obs[:] .+   Z08 .* g08_lc_disc[k].obs[:] 

        g08_lc_conn[k].obs[:] = Z8 .* g08_lc_conn[k].obs[:] .+ Z08 .* g00_lc_conn[k].obs[:]
    end
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

## COMPUTE FVC and TMR

const Qgev = [1.] # Q^2
#=========== COMPUTE FVC  ==============#
@info("Computing Finite-Volume Corrections")
pi_fvc_pi = [Vector{uwreal}(undef, length(Qgev)) for k in eachindex(ensinfo)]
pi_fvc_k  = [Vector{uwreal}(undef, length(Qgev)) for k in eachindex(ensinfo)]

for (k,ens) in enumerate(ensinfo)
    println("Ensemble: $(ens.id)")
    Qlat = Qgev .* value.(t0sqrt_ph.^2 ./ t0ens[k]) ./hc^2 *1e6

    fvc_raw_pi = get_fvc(joinpath(path_fvc, "JKMPI" ), ens.id)
    fvc_pi = vcat(sum(fvc_raw_pi, dims=1)...)
    
    fvc_raw_k = get_fvc(joinpath(path_fvc, "JKMK" ), ens.id)
    fvc_k = vcat(sum(fvc_raw_k, dims=1)...)
    
    T = Float64.(collect(1:length(fvc_pi)))

    for (j,q) in enumerate(Qlat)
        krnl = KRNL.(T,q)
        pi_fvc_pi[k][j]  = abs(sum(fvc_pi .* krnl))
        pi_fvc_k[k][j]   = abs(sum(fvc_k  .* krnl))
    end
end

##
#=========== COMPUTE TMR ==============#
@info("Computing TMR (08)")
# set 1 improvement coefficients
pi_88_ll_conn = [Vector{uwreal}(undef, length(Qgev)) for k in eachindex(ensinfo)]
pi_88_lc_conn = [Vector{uwreal}(undef, length(Qgev)) for k in eachindex(ensinfo)]

pi_88_ll_disc = [Vector{uwreal}(undef, length(Qgev)) for k in eachindex(ensinfo)]
pi_88_lc_disc = [Vector{uwreal}(undef, length(Qgev)) for k in eachindex(ensinfo)]

pi_08_lc_conn = [Vector{uwreal}(undef, length(Qgev)) for k in eachindex(ensinfo)]
pi_08_lc_disc = [Vector{uwreal}(undef, length(Qgev)) for k in eachindex(ensinfo)]

for (k, ens) in enumerate(ensinfo)
    println("Ensemble: ", ens.id)

    Qlat  = Qgev .* value.(t0sqrt_ph.^2) ./ t0ens[k] ./ hc^2 * 1e6

    
    for (j,q) in enumerate(Qlat)
        pi_88_ll_conn[k][j] =  tmr_integrand(g88_ll_conn[k], q, KRNL, pl=true, t0ens=t0ens[k]) *1e5
        pi_88_lc_conn[k][j] =  tmr_integrand(g88_lc_conn[k], q, KRNL, pl=true, t0ens=t0ens[k]) *1e5

        pi_88_ll_disc[k][j] =  tmr_integrand(g88_ll_disc[k], q, KRNL, pl=true, t0ens=t0ens[k], wind=Window("SID")) *1e5
        pi_88_lc_disc[k][j] =  tmr_integrand(g88_lc_disc[k], q, KRNL, pl=true, t0ens=t0ens[k], wind=Window("SID")) *1e5

        pi_08_lc_conn[k][j] =  tmr_integrand(g08_lc_conn[k], q, KRNL, pl=true, t0ens=t0ens[k]) *1e5
        pi_08_lc_disc[k][j] =  tmr_integrand(g08_lc_disc[k], q, KRNL, pl=true, t0ens=t0ens[k], wind=Window("SID")) *1e5
    end
end
## add finite volume
for (k, ens) in enumerate(ensinfo)
    for (j,q) in enumerate(Qgev)
        if ens.kappa_l == ens.kappa_s
            pi_88_ll_conn[k][j] += 1.5 * pi_fvc_pi[k][j] * 1e5 
            pi_88_lc_conn[k][j] += 1.5 * pi_fvc_pi[k][j] * 1e5
        else
            pi_88_ll_conn[k][j] += (1/3 * pi_fvc_pi[k][j] + 2/3 * pi_fvc_k[k][j]) * 1e5 
            pi_88_lc_conn[k][j] += (1/3 * pi_fvc_pi[k][j] + 2/3 * pi_fvc_k[k][j]) * 1e5

            pi_88_ll_disc[k][j] -= (1/3 * pi_fvc_pi[k][j] + 2/3 * pi_fvc_k[k][j]) * 1e5 
            pi_88_lc_disc[k][j] -= (1/3 * pi_fvc_pi[k][j] + 2/3 * pi_fvc_k[k][j]) * 1e5

            pi_08_lc_conn[k][j] += (1/sqrt(3) * pi_fvc_pi[k][j] +  1/sqrt(3) * pi_fvc_k[k][j]) * 1e5

            pi_08_lc_disc[k][j] -= (1/sqrt(3) * pi_fvc_pi[k][j] +  1/sqrt(3) * pi_fvc_k[k][j]) * 1e5

        end
    end
end

##
[uwerr.(pi_88_ll_conn[l]) for l in eachindex(pi_88_ll_conn)]
[uwerr.(pi_88_lc_conn[l]) for l in eachindex(pi_88_lc_conn)]

[uwerr.(pi_88_ll_disc[l]) for l in eachindex(pi_88_ll_disc)]
[uwerr.(pi_88_lc_disc[l]) for l in eachindex(pi_88_lc_disc)]

[uwerr.(pi_08_lc_conn[l]) for l in eachindex(pi_08_lc_conn)]
[uwerr.(pi_08_lc_disc[l]) for l in eachindex(pi_08_lc_disc)]
pi_08_lc_conn
pi_08_lc_disc
sum_cd  = pi_08_lc_conn .+ pi_08_lc_disc
[uwerr.(sum_cd[l]) for l in eachindex(pi_08_lc_disc)]

for (k,ens) in enumerate(ensinfo)
    println("ens: ", ens.id)
    println("88 ll conn: ", pi_88_ll_conn[k])
    println("88 lc conn: ", pi_88_lc_conn[k])

    println("88 ll disc: ", pi_88_ll_disc[k])
    println("88 lc disc: ", pi_88_lc_disc[k])

    println("08 conn:    ", pi_08_lc_conn[k])
    println("08 disc:    ", pi_08_lc_disc[k])
    println("08 tot :    ", sum_cd[k])
    println("")
end

## plot correlator

yyc = g08_lc_conn[1].obs; uwerr.(yyc)
yyd = g08_lc_disc[1].obs; uwerr.(yyd)
xx = collect(1:length(yyc))
# errorbar(xx, value.(yyc), err.(yyc), fmt="s")
errorbar(xx, value.(yyd), err.(yyd), fmt="s")
# xlim(0, Int64(xx[end]/2))
display(gcf())
close("all")