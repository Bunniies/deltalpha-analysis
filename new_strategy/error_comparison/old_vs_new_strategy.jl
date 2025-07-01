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
include("../../Running/tools_running.jl")



path_new_strategy = "/Users/alessandroconigli/MyDrive/postdoc-mainz/projects/deltalpha/physical_results/new_strategy/"
path_old_strategy = "/Users/alessandroconigli/MyDrive/postdoc-mainz/projects/deltalpha/physical_results/scale_error_artificial"

path_plot = "/Users/alessandroconigli/MyDrive/postdoc-mainz/projects/deltalpha/plots/new_strategy/error_comparison"

############################
## PI 33
###########################
## reading old strategy
@warn("The subtracted piece is not added back here. It almost does not contribute to the error, but when doing the running it should be considered")
Qgev_old = [0.05, 0.1, 0.4, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0] # Q^2
Qgev_old =vcat(Qgev_old ./4, Qgev_old)
pi33_lowq_old = read_phys_res(joinpath(path_old_strategy,"piQ4_piQ0"), "PI33_physRes.bdio")
# add_t0_err!(pi33_lowq_old, t0sqrt_ph_err)

pi33_highq_SD_old = read_phys_res(joinpath(path_old_strategy,"piQ_Q4"), "PI33_SD_physRes.bdio")
pi33_highq_ILD_old = read_phys_res(joinpath(path_old_strategy,"piQ_Q4"), "PI33_ILD_physRes.bdio")
pi33_highq_old = pi33_highq_SD_old .+ pi33_highq_ILD_old
# add_t0_err!(pi33_highq_old, t0sqrt_ph_err)

pi33_tot_old = vcat([pi33_lowq_old, pi33_lowq_old .+ pi33_highq_old]...)
uwerr.(pi33_tot_old)
uwerr.(pi33_highq_old)

## reading new strategy 
Qgev_new = [4, 5, 6, 7, 8, 9, 12] # Q^2
Qgev_new =vcat(Qgev_new ./16, Qgev_new ./4, Qgev_new)

pi33_SD_SD_new = read_phys_res(joinpath(path_new_strategy, "Q_SD"),"PI33_SD_physRes.bdio" )
pi33_SD_ILD_new = read_phys_res(joinpath(path_new_strategy, "Q_SD"),"PI33_ILD_physRes.bdio" )
pi33_SD_new = pi33_SD_SD_new .+ pi33_SD_ILD_new
# add_t0_err!(pi33_SD_new, t0sqrt_ph_err)

pi33_ID_new = read_phys_res(joinpath(path_new_strategy, "Q_ID"),"PI33_ID_physRes.bdio" )
# add_t0_err!(pi33_ID_new, t0sqrt_ph_err)

pi33_LD_new = read_phys_res(joinpath(path_new_strategy, "Q_LD"),"PI33_LD_physRes.bdio" )
# add_t0_err!(pi33_LD_new, t0sqrt_ph_err)

pi33_tot_new = vcat([pi33_LD_new, pi33_LD_new .+ pi33_ID_new, pi33_LD_new .+ pi33_ID_new + pi33_SD_new]...)
uwerr.(pi33_tot_new)

##############
## PLOT PI33
##############
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
ylim(0,2.5)
tight_layout()
display(fig)
# savefig(joinpath(path_plot, "pi33_new_vs_old_strategy_err_comparison.pdf"))
close("all")

#########################
## PI88
#########################
# reading old strategy
pi3388_dlt_lowq_old = read_phys_res(joinpath(path_old_strategy, "piQ4_piQ0"), "PI3388_SD_physRes.bdio")
pi33_lowq_SD_old = read_phys_res(joinpath(path_old_strategy, "piQ4_piQ0"), "PI33_SD_physRes.bdio")
pi88_lowq_ILD_old = read_phys_res(joinpath(path_old_strategy, "piQ4_piQ0"), "PI88_ILD_physRes.bdio")
pi88_lowq_old = pi33_lowq_SD_old * charge_factor["88"] - pi3388_dlt_lowq_old + pi88_lowq_ILD_old
add_t0_err!(pi88_lowq_old, t0sqrt_ph_err)
uwerr.(pi88_lowq_old)

pi3388_dlt_highq_old = read_phys_res(joinpath(path_old_strategy, "piQ_Q4"), "PI3388_physRes.bdio")
pi88_highq_old = pi33_highq_old * charge_factor["88"] - pi3388_dlt_highq_old
add_t0_err!(pi88_highq_old, t0sqrt_ph_err)

pi88_tot_old = vcat([pi88_lowq_old, pi88_lowq_old + pi88_highq_old]...)
uwerr.(pi88_tot_old)


# reading new strategy

pi3388_dlt_QSD_new = read_phys_res(joinpath(path_new_strategy, "Q_SD"), "PI3388_QSD_physRes.bdio")
pi88_QSD_new = pi33_SD_new * charge_factor["88"] - pi3388_dlt_QSD_new
add_t0_err!(pi88_QSD_new, t0sqrt_ph_err)

pi3388_dlt_QID_new = read_phys_res(joinpath(path_new_strategy, "Q_ID"), "PI3388_QID_physRes.bdio")
pi88_QID_new = pi33_ID_new * charge_factor["88"] - pi3388_dlt_QID_new
add_t0_err!(pi88_QID_new, t0sqrt_ph_err)

pi88_QLD_new = read_phys_res(joinpath(path_new_strategy, "Q_LD"), "PI88_LD_physRes.bdio")
add_t0_err!(pi88_QLD_new, t0sqrt_ph_err)
pi88_tot_new = vcat([pi88_QLD_new, pi88_QLD_new .+ pi88_QID_new, pi88_QLD_new .+ pi88_QID_new .+ pi88_QSD_new]...)
uwerr.(pi88_tot_new)

## PLOT PI88
fig = figure(figsize=(10,7))

prec_old_lowq = err.(pi88_lowq_old) ./ value.(pi88_lowq_old) .* 100
plot(Qgev_old[1:12], prec_old_lowq, color="red", marker="<", lw=0.0, ms=10, mfc="None")

prec_old_highq = err.(pi88_tot_old) ./ value.(pi88_tot_old) .* 100
plot(Qgev_old[13:end], prec_old_highq[13:end], color="red", marker="<", label="Old strategy", lw=0.0, ms=10)


prec_new = err.(pi88_tot_new) ./ value.(pi88_tot_new) .* 100    
scatter(Qgev_new, prec_new, color="royalblue", marker="^", label="New strategy", 80)

xlabel(L"$Q^2 \ [\mathrm{GeV}^2]$")
ylabel(L"$\delta \Pi^{(8,8)}(Q^2) \ [\%]$")
legend()
ylim(0,3.0)
tight_layout()
display(fig)
# savefig(joinpath(path_plot, "pi88_new_vs_old_strategy_err_comparison.pdf"))
close("all")

##########
## PI cc
##########
# reading old strategy
picc_old_lowq = read_phys_res(joinpath(path_old_strategy, "PIQ4_piQ0"), "Picc_conn_physRes.bdio")
add_t0_err!(picc_old_lowq, t0sqrt_ph_err); uwerr.(picc_old_lowq)
picc_old_highq = read_phys_res(joinpath(path_old_strategy, "PIQ_Q4"), "Picc_conn_physRes.bdio")
add_t0_err!(picc_old_highq, t0sqrt_ph_err); uwerr.(picc_old_highq)
picc_old_tot = vcat(picc_old_lowq, picc_old_lowq + picc_old_lowq)
uwerr.(picc_old_tot)

# reading new strategy
picc_QSD_new = read_phys_res(joinpath(path_new_strategy, "Q_SD"), "PIcc_conn_QSD_physRes.bdio")
add_t0_err!(picc_QSD_new, t0sqrt_ph_err); 
uwerr.(picc_QSD_new)

picc_QID_new = read_phys_res(joinpath(path_new_strategy, "Q_ID"), "PIcc_conn_QID_physRes.bdio")
add_t0_err!(picc_QID_new, t0sqrt_ph_err); 
uwerr.(picc_QID_new)

picc_QLD_new = read_phys_res(joinpath(path_new_strategy, "Q_LD"), "PIcc_conn_QLD_physRes.bdio")
add_t0_err!(picc_QLD_new, t0sqrt_ph_err); 
uwerr.(picc_QLD_new)

picc_new_tot = vcat(picc_QLD_new, picc_QLD_new + picc_QID_new,picc_QSD_new + picc_QID_new + picc_QLD_new)
uwerr.(picc_new_tot)

## PLOT PICC
fig = figure(figsize=(10,7))

prec_old_lowq = err.(picc_old_lowq) ./ value.(picc_old_lowq) .* 100
plot(Qgev_old[1:12], prec_old_lowq, color="red", marker="<", lw=0.0, ms=10, mfc="None")

prec_old_highq = err.(picc_old_tot) ./ value.(picc_old_tot) .* 100
plot(Qgev_old[13:end], prec_old_highq[13:end], color="red", marker="<", label="Old strategy", lw=0.0, ms=10)


prec_new = err.(picc_new_tot) ./ value.(picc_new_tot) .* 100    
scatter(Qgev_new, prec_new, color="royalblue", marker="^", label="New strategy", 80)

xlabel(L"$Q^2 \ [\mathrm{GeV}^2]$")
ylabel(L"$\delta \Pi^{(c,c)}(Q^2) \ [\%]$")
legend()
ylim(0,3.0)
tight_layout()
display(fig)
savefig(joinpath(path_plot, "picc_new_vs_old_strategy_err_comparison.pdf"))
close("all")


######################
## PI 08
######################

# reading old strategy
pi08_highq_old = read_phys_res(joinpath(path_old_strategy, "piQ_Q4"), "PI08_physRes.bdio")
add_t0_err!(pi08_highq_old, t0sqrt_ph_err); uwerr.(pi08_highq_old)

pi08_lowq_old = read_phys_res(joinpath(path_old_strategy, "piQ4_piQ0"), "PI08_physRes.bdio")
add_t0_err!(pi08_lowq_old, t0sqrt_ph_err); uwerr.(pi08_lowq_old)
pi08_tot_old = vcat(pi08_lowq_old, pi08_lowq_old + pi08_highq_old); uwerr.(pi08_tot_old)

# reading new strategy

pi08_QSD_new = read_phys_res(joinpath(path_new_strategy, "Q_SD"),"PI08_Q_SD_physRes.bdio" )
add_t0_err!(pi08_QSD_new, t0sqrt_ph_err); uwerr.(pi08_QSD_new)

pi08_QID_new = read_phys_res(joinpath(path_new_strategy, "Q_ID"),"PI08_Q_ID_physRes.bdio" )
add_t0_err!(pi08_QID_new, t0sqrt_ph_err); uwerr.(pi08_QID_new)

pi08_QLD_new = read_phys_res(joinpath(path_new_strategy, "Q_LD"),"PI08_Q_LD_physRes.bdio")
add_t0_err!(pi08_QLD_new, t0sqrt_ph_err); uwerr.(pi08_QLD_new)

pi08_tot_new = vcat(pi08_QLD_new, pi08_QLD_new + pi08_QID_new, pi08_QLD_new + pi08_QID_new + pi08_QSD_new)
uwerr.(pi08_tot_new)

## PLOT PI08
fig = figure(figsize=(10,7))

prec_old_lowq = err.(pi08_lowq_old) ./ value.(pi08_lowq_old) .* 100
plot(Qgev_old[1:12], prec_old_lowq, color="red", marker="<", lw=0.0, ms=10, mfc="None")

prec_old_highq = err.(pi08_tot_old) ./ value.(pi08_tot_old) .* 100
plot(Qgev_old[13:end], prec_old_highq[13:end], color="red", marker="<", label="Old strategy", lw=0.0, ms=10)


prec_new = err.(pi08_tot_new) ./ value.(pi08_tot_new) .* 100    
scatter(Qgev_new, prec_new, color="royalblue", marker="^", label="New strategy", 80)

xlabel(L"$Q^2 \ [\mathrm{GeV}^2]$")
ylabel(L"$\delta \Pi^{(0,8)}(Q^2) \ [\%]$")
legend()
ylim(0,3.0)
tight_layout()
display(fig)
savefig(joinpath(path_plot, "pi08_new_vs_old_strategy_err_comparison.pdf"))
close("all")