######################################################################################
# This file was created by Alessandro Conigli 
# Computation of mpi and fpi 
########################################################################################

using Revise
using HVPobs
using PyPlot, LaTeXStrings
using OrderedCollections
using TimerOutputs
using BDIO, ADerrors
import ADerrors: err

#================ SET UP VARIABLES ===============#

rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
rcParams["text.usetex"] =  true
rcParams["mathtext.fontset"]  = "cm"
rcParams["font.size"] =13
rcParams["axes.labelsize"] =22
rcParams["axes.titlesize"] = 18

plt.rc("text", usetex=true) # set to true if a LaTeX installation is present

#======== INCLUDE, PATHS AND CONSTANTS ==========#
include("../utils/const.jl")
include("../utils/types.jl")
include("../utils/tools.jl")
# include("data_management.jl")
include("../utils/plot_utils.jl")

path_data = "/Users/alessandroconigli/Lattice/data/HVP/mesons"
path_rw   = "/Users/alessandroconigli/Lattice/data/HVP/rwf_deflated"
path_ms   = "/Users/alessandroconigli/Lattice/data/HVP/ms_t0_dat"
path_fvc  = "/Users/alessandroconigli/Lattice/data/HVP/FVC"

path_plot = "/Users/alessandroconigli/MyDrive/postdoc-mainz/projects/deltalpha/plots/ensembles/"
path_bdio = "/Users/alessandroconigli/MyDrive/postdoc-mainz/projects/deltalpha/data"

IMPR      = true
RENORM    = true
STD_DERIV = false

enslist = ["D200"]
ensinfo = EnsInfo.(enslist)


#============ OBSERVABLE ALLOCATIONS ============#

g_pi = Vector{Corr}(undef, length(ensinfo))
g_k  = Vector{Corr}(undef, length(ensinfo))

g_pi_a0p = Vector{Corr}(undef, length(ensinfo))
g_k_a0p = Vector{Corr}(undef, length(ensinfo))

obs = Vector(undef, length(ensinfo))

#================ READING DATA ====================#

@info("Reading data...")
@time begin
    
    for (k, ens) in enumerate(ensinfo)
        println("\n - Ensemble: $(ens.id)")

        println("   - Pion correlator")
        g_pi[k] = get_mesons_corr(path_data, ens, "ll", "PP", path_rw=path_rw, frw_bcwd=false, L=1 )
        g_pi[k].obs .*= -1
        g_pi_a0p[k] = get_mesons_corr(path_data, ens, "ll", "A0P", path_rw=path_rw, frw_bcwd=false, L=1 )
        
        improve_corr_vkvk!(g_pi_a0p[k], g_pi[k], -ca(ens.beta))
    end
end

##
#======== COMPUTE OBS ============#
#mpi
mpi = meff(g_pi[1], [40,80], pl=true)

#fpi (the function dec_const_ratio is defiend in line 93)

_, y0_ens = findmax(value.(g_pi[1].obs)) # get the source
fpi = dec_const_ratio(g_pi_a0p[1].obs, g_pi[1].obs, mpi, [40,90], y0_ens-1, pl=true) * Za_l_sub(ensinfo[1].beta)
uwerr(fpi)

println("#==== Results ===#")
println("  mpi = ", print_uwreal(mpi))
println("  fpi = ", print_uwreal(fpi))
println("#=== End ===#")

########################
## function definition
#######################
function dec_const_ratio(a0p::Vector{uwreal}, pp::Vector{uwreal}, m::uwreal, plat::Vector{Int64}, y0::Int64; pl::Bool=true)

    a0p = a0p
    T = length(a0p)

    R = (((a0p .* reverse(a0p)) ./ (pp[T-y0]) ).^2).^0.25

    R_av = plat_av(R[2:end-1], plat)

    if pl
        uwerr.(R[2:end-1])
        v = value.(R[2:end-1])
        e = err.(R[2:end-1])
        uwerr(R_av)
        errorbar(collect(1:length(v)), v, e, fmt="s", color="black")
        fill_between(plat, value(R_av)- err(R_av), value(R_av)+err(R_av), alpha=0.5, color="royalblue")
        ylim(value(R_av)*0.9, value(R_av)*1.1)
        display(gcf())
        close("all")
    end

    return sqrt(2/m)  * R_av
end