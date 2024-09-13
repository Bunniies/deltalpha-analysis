using HVPobs
using ALPHAio, BDIO
using PyPlot, LaTeXStrings
using ADerrors
import ADerrors:err
using DelimitedFiles
#================ SET UP PLOT PARAMETERS ==================#

rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
rcParams["text.usetex"] =  true
rcParams["mathtext.fontset"]  = "cm"
rcParams["font.size"] =13
rcParams["axes.labelsize"] =22
rcParams["axes.titlesize"] = 18
plt.rc("text", usetex=true) # set to true if a LaTeX installation is present

path_phys_res_highq = "/Users/alessandroconigli/MyDrive/postdoc-mainz/projects/deltalpha/physical_results/scale_error_artificial/piQ_Q4/"
path_phys_res_lowq = "/Users/alessandroconigli/MyDrive/postdoc-mainz/projects/deltalpha/physical_results/scale_error_artificial/piQ4_piQ0/"
path_phys_res_fullq = "/Users/alessandroconigli/MyDrive/postdoc-mainz/projects/deltalpha/physical_results/scale_error_multimom/full_q/"
path_plot = "/Users/alessandroconigli/MyDrive/postdoc-mainz/projects/deltalpha/plots/running"

include("../utils/IO_BDIO.jl")
include("./tools_running.jl")
include("../utils/const.jl")


const Qgev = [0.05, 0.1, 0.4, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0] # Q^2
NMOM = length(Qgev)
Qmgev =  9.0  # GeV^2

######################################
## HIGH Q
#####################################
## read phys res data HIGH Q

# 33 PT
b33qqm_highq = [(q/(4*Qmgev)) * log(2) / (4pi^2) for q in Qgev]

# cc connected
picc_highq_sub = read_phys_res(path_phys_res_highq, "PIcc_conn_physRes.bdio")
bcc_highq = read_phys_res(path_phys_res_highq, "PI_blc_2Qm_physRes.bdio")

picc_highq = picc_highq_sub .+ Qgev ./ (4*Qmgev) .* bcc_highq .+ 2 .* charge_factor["cc"] * b33qqm_highq

######################################
## LOW Q
#####################################
## read phys res data LOW Q
b33qqm_lowq = [(q/(4*Qmgev)) * log(2) / (4pi^2) for q in Qgev ./4]
# cc connected
picc_lowq_sub = read_phys_res(path_phys_res_lowq, "PIcc_conn_physRes.bdio")
picc_lowq = picc_highq_sub .+ 4/3 * Qgev ./4 ./ (4*Qmgev) .* bcc_highq .+ 2 .* charge_factor["cc"] .* b33qqm_lowq

######################
## TOTAL Q
#####################
syst_fullq = read_systematics(path_phys_res_fullq)
picc_tot = read_phys_res(path_phys_res_fullq, "PIcc_conn_physRes.bdio")
picc_tot = picc_tot .+ [uwreal([0.0, syst_fullq["piccconn"][k]], "syst cc fullq") for k in eachindex(Qgev)]
picc_tot ./= charge_factor["cc"]

## comparison of charm plot
picc_tot_sum = picc_highq .+ picc_lowq
picc_tot_sum ./= charge_factor["cc"]

funccc = rational_func(M=2, N=3)
funccc_fullq = rational_func(M=2, N=3)

wpmm = Dict{String, Vector{Float64}}()
wpmm["E300"] = [-1., 2., -1., -1.]
fitcc = fit_routine(funccc, Qgev, picc_tot_sum, 10, wpm=wpmm)
fitcc_fullq = fit_routine(funccc_fullq, Qgev, picc_tot, 10, wpm=wpmm)

ax = gca()
ax.tick_params(right=true)

# cc connected
errorbar(Qgev, value.(picc_tot_sum), err.(picc_tot_sum), fmt="d", color="darkorange", mfc="none", ms=5)
xarr = range(0,9,length=100)
yarrcc = funccc(xarr, fitcc.param); uwerr.(yarrcc)
fill_between(xarr, value.(yarrcc).-err.(yarrcc), value.(yarrcc).+err.(yarrcc), alpha=0.5, label=L"$\mathrm{SD\ strategy}$", color="darkorange")
plot(xarr, value.(yarrcc), color="darkorange")

# cc connected full q
errorbar(Qgev, value.(picc_tot), err.(picc_tot), fmt="d", color="red", mfc="none", ms=5)
xarr = range(0,9,length=100)
yarrcc = funccc(xarr, fitcc_fullq.param); uwerr.(yarrcc)
fill_between(xarr, value.(yarrcc).-err.(yarrcc), value.(yarrcc).+err.(yarrcc), alpha=0.5, label=L"$\mathrm{2022\ strategy}$", color="red")
plot(xarr, value.(yarrcc), color="red")

xlabel(L"$Q^2\ [\mathrm{GeV}^2]$")
ylabel(L"$\mathit{\bar\Pi}^{cc}(-Q^2)$")
legend()
tight_layout()
display(gcf())
savefig(joinpath(path_plot, "charm_comparison_totq.pdf"))
close("all")

## error from t0
uwerr(t0sqrt_ph_err)
mchist(picc_tot_sum[12], "sqrtt0 [fm]") * 1e8 * err(t0sqrt_ph_err)
