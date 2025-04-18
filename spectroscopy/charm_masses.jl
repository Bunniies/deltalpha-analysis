######################################################################################
# This file was created by Alessandro Conigli 
# Here we compute the finite-a pion and kaon masses required in the chiral-continuum extrapolations
# The observables are then stored in a BDIO file, following the order 
########################################################################################

using Revise
using HVPobs
using PyPlot, LaTeXStrings
using OrderedCollections
using TimerOutputs
using BDIO, ADerrors
using ALPHAio
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
include("../utils/IO_BDIO.jl")

path_charm_corr = "/Users/alessandroconigli/Lattice/data/HVP/charm_correlators"

path_plot = "/Users/alessandroconigli/MyDrive/postdoc-mainz/projects/deltalpha/plots/ensembles/"
path_bdio = "/Users/alessandroconigli/MyDrive/postdoc-mainz/projects/deltalpha/data"

# enslist = sort([ "H101", "H102", "N101", "C101",
        #  "B450", "D450", "D452",
        #  "N202", "N203", "N200", "D200",  "E250",
        #  "N300", "J303", "E300",
        #  "J500"
# ])
enslist = sort(["D450"])
ensinfo = EnsInfo.(enslist)

path_ens = vcat([filter(x-> occursin(enslist[k], basename(x)), readdir(path_charm_corr, join=true)) for k in eachindex(enslist)]...)
##
# READ CORRELATORS

gcc_k1 = Vector{Vector{uwreal}}(undef, length(enslist))
gcc_k2 = Vector{Vector{uwreal}}(undef, length(enslist))
gcc_k3 = Vector{Vector{uwreal}}(undef, length(enslist))
gcc_k4 = Vector{Vector{uwreal}}(undef, length(enslist))

for (k,ens) in enumerate(ensinfo)

    fb = BDIO_open(path_ens[k], "r")
    res = Dict()
    while ALPHAdobs_next_p(fb)
        d = ALPHAdobs_read_parameters(fb)
        sz = tuple(d["size"]...)
        ks = collect(d["keys"])
        res = ALPHAdobs_read_next(fb, size=sz, keys=ks)
    end
    BDIO_close!(fb)
    gcc_k1[k] = res["sh1"]#[1:end-2]
    gcc_k2[k] = res["sh2"]#[1:end-2]
    gcc_k3[k] = res["sh3"]#[1:end-2]
    gcc_k4[k] = res["sh4"]#[1:end-2]
end

## COMPUTE MASSES
obs = Vector(undef, length(ensinfo))
# mk1 = meff(gcc_k1[1], [10,40], pl=true)

wpmm = Dict{String, Vector{Float64}}()
wpmm["E300"]     = [5.0, -2.0, -1.0, -1.0]
wpmm["J500"]     = [5.0, -2.0, -1.0, -1.0]
wpmm["H400"]     = [5.0, -1.5, -1.0, -1.0]
wpmm["N202"]     = [5.0, -2.0, -1.0, -1.0]
wpmm["N200"]     = [5.0, -2.0, -1.0, -1.0]
wpmm["N203"]     = [5.0, -2.0, -1.0, -1.0]
wpmm["N300"]     = [5.0, -1.5, -1.0, -1.0]
wpmm["J303"]     = [5.0, -2.0, -1.0, -1.0]
wpmm["D450"]     = [2.0, -2.0, -1.0, -1.0]
wpmm["D452"]     = [2.0, -2.0, -1.0, -1.0]
wpmm["F300"]     = [5.0, -2.0, -1.0, -1.0]
wpmm["J306"]     = [5.0, -2.0, -1.0, -1.0]
wpmm["J307"]     = [5.0, -2.0, -1.0, -1.0]
for (k,ens) in enumerate(ensinfo)
    obs[k] = OrderedDict()
    obs[k]["sh1"] = get_meff_BMA([Corr(gcc_k1[k], ens.id, "mDs")], ens, path_plt=path_plot, ll=L"$m_{D_s^{(1)}}$", wpm=wpmm)
    obs[k]["sh2"] = get_meff_BMA([Corr(gcc_k2[k], ens.id, "mDs")], ens, path_plt=path_plot, ll=L"$m_{D_s^{(2)}}$", wpm=wpmm)
    obs[k]["sh3"] = get_meff_BMA([Corr(gcc_k3[k], ens.id, "mDs")], ens, path_plt=path_plot, ll=L"$m_{D_s^{(3)}}$", wpm=wpmm)
    obs[k]["sh4"] = get_meff_BMA([Corr(gcc_k4[k], ens.id, "mDs")], ens, path_plt=path_plot, ll=L"$m_{D_s^{(4)}}$", wpm=wpmm)

end

## 
#============== SAVE results to BDIO ===============#


for (k, ens) in enumerate(ensinfo)
    io = IOBuffer()
    write(io, "mDs masses. ")
    fb = ALPHAdobs_create(joinpath(path_bdio, ens.id, "mDs_masses.bdio"), io)

    extra = Dict{String, Any}("Ens" => ens.id)
    data = Dict{String, Array{uwreal}}(
        "sh1" => obs[k]["sh1"],
        "sh2" => obs[k]["sh2"],
        "sh3" => obs[k]["sh3"],
        "sh4" => obs[k]["sh4"]
    )

    ALPHAdobs_write(fb, data, extra=extra)
    ALPHAdobs_close(fb)
end

## TEST READING
res = Dict()
fb = BDIO_open(joinpath(path_bdio,"H101","mDs_masses.bdio"), "r")
while ALPHAdobs_next_p(fb)
    d = ALPHAdobs_read_parameters(fb)
    sz = tuple(d["size"]...)
    ks = collect(d["keys"])
    res = ALPHAdobs_read_next(fb, size=sz, keys=ks)
end
BDIO_close!(fb)

##
tt = read_phys_res("/Users/alessandroconigli/MyDrive/postdoc-mainz/projects/deltalpha/data/C101","mDs_masses.bdio")