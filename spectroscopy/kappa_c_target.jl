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

include("../utils/const.jl")
include("../utils/types.jl")
include("../utils/plot_utils.jl")
include("../utils/IO_BDIO.jl")
include("../utils/tools.jl")

path_plot = "/Users/alessandroconigli/MyDrive/postdoc-mainz/projects/deltalpha/plots/ensembles/"
path_bdio = "/Users/alessandroconigli/MyDrive/postdoc-mainz/projects/deltalpha/data"

enslist = sort([ "H101", "H102", "N101", "C101",
         "B450", "D450", "D452",
         "N202", "N203", "N200", "D200",  "E250",
         "N300", "J303", "E300",
         "J500"
])

ensinfo = EnsInfo.(enslist)

## READ MASSES & KAPPA VALUES

masses = Dict()
for (k,ens) in enumerate(ensinfo)
    p = joinpath(path_bdio, ens.id, "mDs_masses.bdio")
    fb = BDIO_open(p, "r")
    while ALPHAdobs_next_p(fb)
        d = ALPHAdobs_read_parameters(fb)
        sz = tuple(d["size"]...)
        ks = collect(d["keys"])
        masses[ens.id] = ALPHAdobs_read_next(fb, size=sz, keys=ks)
    end
    BDIO_close!(fb)
end

kappas = Dict()
kappas_aux = get_kappa_values()
for (k,ens) in enumerate(ensinfo)
    kappas[ens.id] = Dict()
    aux = sort(kappas_aux[ens.id][2:end], rev=true)
    kappas[ens.id]["sh1"] = aux[1]
    kappas[ens.id]["sh2"] = aux[2]
    kappas[ens.id]["sh3"] = aux[3]
    kappas[ens.id]["sh4"] = aux[4]
end

## INTERPOLATION 

mDs_ph = uwreal([1983.4, 1e-8], "mDs phys") # = 1968.47 * sqrtt0_bruno / sqrtt0_regensburg
kappa_c_target_tot = Dict{String, Array{uwreal}}()
for (k,ens) in enumerate(ensinfo)

    keystot = collect(keys(masses[ens.id]))
    m_ens = vcat([masses[ens.id][ll] for ll in keystot]...)
    idx = sortperm(value.(m_ens))
    m_ens = m_ens[idx]
    uwerr.(m_ens)
    kappa_ens_inv = sort(1 ./ [kappas[ens.id][ll] for ll in keystot])
    errorbar(kappa_ens_inv, value.(m_ens), err.(m_ens), fmt="s", color="black", mfc="none", capsize=2)

    fitp, chi2exp = lin_fit(kappa_ens_inv, m_ens)
    # mDs_lat = mDs_ph * t0sqrt_ph / sqrt(t0(ens.beta)) / hc
    mDs_lat = mDs_ph * t0sqrt_ph / value(sqrt(t0_bruno(ens.beta))) / hc
    kappa_c_target = x_lin_fit(fitp, mDs_lat)
    uwerr(kappa_c_target)

    kappa_c_target_tot[ens.id] = [1 / kappa_c_target]
    uwerr(kappa_c_target_tot[ens.id][1])

    errorbar(value(kappa_c_target), value(mDs_lat), xerr=err(kappa_c_target), fmt="d", capsize=2, color="royalblue")
    xarr = range(minimum(kappa_ens_inv), maximum(kappa_ens_inv), 100)
    yarr = Vector{uwreal}(undef, 100) #similar(xarr)
    for l in 1:100
        yarr[l] = fitp[1]  + fitp[2] * xarr[l]
    end
    uwerr.(yarr)
    fill_between(xarr, value.(yarr).-err.(yarr), value.(yarr).+err.(yarr), alpha=0.5, color="tomato")

    xlabel(L"$1/\kappa_c$")
    ylabel(L"$m_{D_s}$")
    tight_layout()
    display(gcf())
    # savefig(joinpath(path_plot, ens.id, "kappa_target_interp.pdf"))
    close("all")
end

## SAVING KAPPA_C_TARGET

io = IOBuffer()
write(io, "kappa_charm_target")
fb = ALPHAdobs_create(joinpath("/Users/alessandroconigli/MyDrive/postdoc-mainz/projects/deltalpha/physical_results", "kappa_c_target.bdio"), io)
ALPHAdobs_write(fb, kappa_c_target_tot)
ALPHAdobs_close(fb)

## test reading
res = Dict()
fb = BDIO_open(joinpath("/Users/alessandroconigli/MyDrive/postdoc-mainz/projects/deltalpha/physical_results", "kappa_c_target.bdio"), "r")
while ALPHAdobs_next_p(fb)
    d = ALPHAdobs_read_parameters(fb)
    sz = tuple(d["size"]...)
    ks = collect(d["keys"])
    res = ALPHAdobs_read_next(fb, size=sz, keys=ks)
end
BDIO_close!(fb)