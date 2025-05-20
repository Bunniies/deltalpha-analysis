# error comparison between old and new strategy with 2 and 3 windows in the momenta region

using Revise, HVPobs
using PyPlot, LaTeXStrings
using BDIO, ADerrors, ALPHAio
using QuadGK
import ADerrors: err
using Statistics

#================ SET UP PLOT PARAMETERS ==================#

rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
rcParams["text.usetex"] =  true
rcParams["mathtext.fontset"]  = "cm"
rcParams["font.size"] = 20
rcParams["axes.labelsize"] =26
rcParams["axes.titlesize"] = 22
plt.rc("text", usetex=true) # set to true if a LaTeX installation is present


include("../../utils/const.jl")
include("../../utils/IO_BDIO.jl")


path_new_strategy = "/Users/alessandroconigli/MyDrive/postdoc-mainz/projects/deltalpha/physical_results/new_strategy/"
path_old_strategy = "/Users/alessandroconigli/MyDrive/postdoc-mainz/projects/deltalpha/physical_results/scale_error_artificial"

path_plot = "/Users/alessandroconigli/MyDrive/postdoc-mainz/projects/deltalpha/plots/new_strategy/error_comparison"
## reading old strategy
Qgev_old = [0.05, 0.1, 0.4, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0] # Q^2
Qgev_old =vcat(Qgev_old ./4, Qgev_old)
pi33_lowq_old = read_phys_res(joinpath(path_old_strategy,"piQ4_piQ0"), "PI33_physRes.bdio")
pi33_highq_SD_old = read_phys_res(joinpath(path_old_strategy,"piQ_Q4"), "PI33_SD_physRes.bdio")
pi33_highq_ILD_old = read_phys_res(joinpath(path_old_strategy,"piQ_Q4"), "PI33_ILD_physRes.bdio")
pi33_highq_old = pi33_highq_SD_old .+ pi33_highq_ILD_old

pi33_tot_old = vcat([pi33_lowq_old, pi33_lowq_old .+ pi33_highq_old]...)
uwerr.(pi33_tot_old)
uwerr.(pi33_highq_old)

## reading new strategy 
Qgev_new = [9.0, 12.0, 15.0, 18.0] # Q^2
Qgev_new =vcat(Qgev_new ./16, Qgev_new ./4, Qgev_new)

pi33_SD_SD_new = read_phys_res(joinpath(path_new_strategy, "Q_SD"),"PI33_SD_physRes.bdio" )
pi33_SD_ILD_new = read_phys_res(joinpath(path_new_strategy, "Q_SD"),"PI33_ILD_physRes.bdio" )
pi33_SD_new = pi33_SD_SD_new .+ pi33_SD_ILD_new

pi33_ID_new = read_phys_res(joinpath(path_new_strategy, "Q_ID"),"PI33_ID_physRes.bdio" )
pi33_LD_new = read_phys_res(joinpath(path_new_strategy, "Q_LD"),"PI33_LD_physRes.bdio" )

pi33_tot_new = vcat([pi33_LD_new, pi33_LD_new .+ pi33_ID_new, pi33_LD_new .+ pi33_ID_new + pi33_SD_new]...)
uwerr.(pi33_tot_new)

##
fig = figure(figsize=(10,7))

prec_old_lowq = err.(pi33_lowq_old) ./ value.(pi33_lowq_old) .* 100
plot(Qgev_old[1:12], prec_old_lowq, color="red", marker="<", lw=0.0, ms=10, mfc="None")

prec_old_highq = err.(pi33_tot_old) ./ value.(pi33_tot_old) .* 100
plot(Qgev_old[13:end], prec_old_highq[13:end], color="red", marker="<", label="Old strategy", lw=0.0, ms=10)


prec_new = err.(pi33_tot_new) ./ value.(pi33_tot_new) .* 100    
scatter(Qgev_new, prec_new, color="royalblue", marker="^", label="New strategy", 80)

xlabel(L"$Q^2 \ [\mathrm{GeV}^2]$")
ylabel(L"$\delta \Pi^{(3,3)}(Q^2) \ [\%]$")
legend()
tight_layout()
display(fig)
savefig(joinpath(path_plot, "new_vs_old_strategy_err_comparison.pdf"))
close("all")