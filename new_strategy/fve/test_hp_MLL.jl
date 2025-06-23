using Revise, HVPobs
using PyPlot, LaTeXStrings
using BDIO, ADerrors, ALPHAio
import ADerrors: err
using Statistics
using DelimitedFiles


#================ SET UP PLOT PARAMETERS ==================#

rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
rcParams["text.usetex"] =  true
rcParams["mathtext.fontset"]  = "cm"
rcParams["font.size"] = 18
rcParams["axes.labelsize"] = 20
rcParams["axes.titlesize"] = 22
plt.rc("text", usetex=true) # set to true if a LaTeX installation is present


include("../../utils/const.jl")
include("../../utils/types.jl")
include("../../utils/plot_utils.jl")
include("../../utils/IO_BDIO.jl")
include("../../utils/tools.jl")


const path_fvc_hp  = "/Users/alessandroconigli/Lattice/data/HVP/FSE"
path_fvc_mll = "/Users/alessandroconigli/MyDrive/postdoc-mainz/projects/deltalpha/llm-fvc/gsFVC"
path_plot = "/Users/alessandroconigli/MyDrive/postdoc-mainz/projects/deltalpha/plots/new_strategy/ensembles/"

KRNL = krnl_dα

ensinfo = EnsInfo.(["D200"])

Qgev = 1.0 
Qmgev = 9.0

##
fvc_hp = Vector{Matrix{uwreal}}(undef, length(ensinfo))
fvc_acc_hp = Vector{Vector{Vector{uwreal}}}(undef, length(ensinfo))

for (k,ens) in enumerate(ensinfo)
    fvc_hp[k] = get_fvc(joinpath(path_fvc_hp, "JKMPI" ), ens.id)
    fvc_acc_hp[k] = Vector{Vector{uwreal}}(undef, 0)

    qlat = Qgev .* t0sqrt_ph.^2 ./ t0(ens.beta) ./hc^2 *1e6
    qmlat = Qmgev  * t0sqrt_ph^2 / t0(ens.beta) /hc^2 *1e6
    T=0

    for (idx,n) in enumerate([1,2,3,6])
        push!(fvc_acc_hp[k],  vcat(sum(fvc_hp[k][1:n,:], dims=1)...))
        T = Float64.(collect(1:length(fvc_acc_hp[k][idx])))
        tmr = -1. *fvc_acc_hp[k][idx] .*  KRNL.(T, qlat)  
        tmp = sum(tmr); 
        # [uwerr(tmp[ll], wpmm) for ll in eachindex(tmp)]
        uwerr.(tmp)
        tmr .=  tmr ./ value.(t0sqrt_ph) .* sqrt.(t0(ens.beta)) # convert to phys units
        uwerr.(tmr)
        fill_between(T, value.(tmr)-err.(tmr), value.(tmr)+err.(tmr), alpha=0.8, lw=50, label=L"$F_\pi(-Q^2),\ \vec{n}^2 \leq$"* " $(n)" )
    end

    ## MLL
    path_fvc_mll_ens = joinpath(path_fvc_mll, "D200_gs_fvc.txt")
    fname = readdlm(path_fvc_mll_ens, comments=true)
    tvals = fname[:,1]
    corr_mll = Vector{uwreal}(undef, length(tvals))
    for t in eachindex(tvals)
        corr_mll[t] = uwreal([fname[t,2], fname[t,3]], "MLL")
    end
    fvc_mll_tmr = -1 * corr_mll .* KRNL.(tvals, qlat)
    fvc_mll_tmr .=  fvc_mll_tmr ./ value.(t0sqrt_ph) .* sqrt.(t0(ens.beta)) # convert to phys units

    uwerr.(fvc_mll_tmr)
    fill_between(tvals, value.(fvc_mll_tmr)-err.(fvc_mll_tmr), value.(fvc_mll_tmr)+err.(fvc_mll_tmr), alpha=0.8, lw=50, label=L"$\mathrm{GS\ } F_\pi(\omega), \ \mathrm{LLM}$" )
    axvline(tvals[1], ls="--", color="black", lw=0.5)

    legend(fontsize=14)
    xlabel(L"$t/a$")
    ylabel(L"$\delta G^{33}(t)K(t,Q^2) \ [\mathrm{fm}^{-1}]$")
    tight_layout()
    xlim(0, Int(T[end]/2))
    display(gcf())
    t = "$(ens.id)_fvc_MLL_HP_1GeV2.pdf"
    savefig(joinpath(path_plot, ens.id ,t))
    close("all")
end