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

path_phys_res_highq = "/Users/alessandroconigli/MyDrive/postdoc-mainz/projects/deltalpha/physical_results/scale_error_multimom/piQ_Q4/"
path_phys_res_lowq = "/Users/alessandroconigli/MyDrive/postdoc-mainz/projects/deltalpha/physical_results/scale_error_multimom/piQ4_piQ0/"
path_phys_res_fullq = "/Users/alessandroconigli/MyDrive/postdoc-mainz/projects/deltalpha/physical_results/scale_error_multimom/full_q/"
path_plot = "/Users/alessandroconigli/MyDrive/postdoc-mainz/projects/deltalpha/plots/running"

include("../utils/IO_BDIO.jl")
include("./tools_running.jl")

const Qgev = [0.05, 0.1, 0.4, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0] # Q^2
NMOM = length(Qgev)
Qmgev =  36.0  # GeV^2

######################################
## HIGH Q
#####################################
## read phys res data HIGH Q
# reading systematic errors

syst_highq = read_systematics(path_phys_res_highq)

# 33
pi33_SD = read_phys_res(path_phys_res_highq, "PI33_SD_physRes.bdio")
pi33_ILD = read_phys_res(path_phys_res_highq, "PI33_ILD_physRes.bdio")
pi33_SD = pi33_SD .+ [uwreal([0.0, syst_highq["pi33SD"][k]], "syst 33 SD") for k in eachindex(Qgev)]
pi33_ILD = pi33_ILD .+ [uwreal([0.0, syst_highq["pi33ILD"][k]], "syst 33 ILD") for k in eachindex(Qgev)]

b33qqm_highq = [(q/(4*Qmgev)) * log(2) / (4pi^2) for q in Qgev]
pi33_highq = (pi33_ILD + pi33_SD + b33qqm_highq ) 

# 88
pi3388dlt_highq = read_phys_res(path_phys_res_highq, "PI3388_physRes.bdio")
pi3388dlt_highq = pi3388dlt_highq .+ [uwreal([0.0, syst_highq["pi3388"][k]], "syst 3388") for k in eachindex(Qgev)]
pi88_highq = pi33_highq * charge_factor["88"] - pi3388dlt_highq

# cc connected 
picc_highq_sub = read_phys_res(path_phys_res_highq, "PIcc_conn_physRes.bdio")
picc_highq_sub = picc_highq_sub .+ [uwreal([0.0, syst_highq["piccconn"][k]], "syst cc conn") for k in eachindex(Qgev)]

bcc_highq = read_phys_res(path_phys_res_highq, "PI_blc_physRes.bdio")
bcc_highq = bcc_highq .+ [uwreal([0.0, syst_highq["piblc"][k]], "syst blc") for k in eachindex(Qgev)]

picc_highq = picc_highq_sub .+ bcc_highq .+ 2 * charge_factor["cc"] * b33qqm_highq

# 08
pi08_highq = read_phys_res(path_phys_res_highq, "PI08_physRes.bdio")
pi08_highq = pi08_highq .+ [uwreal([0.0, syst_highq["pi80"][k]], "syst 80") for k in eachindex(Qgev)]

# remove charge factor
pi88_highq ./= charge_factor["88"]
picc_highq ./= charge_factor["cc"]
pi08_highq ./= charge_factor["08"]
## fit
func33 = rational_func(M=2, N=3)
func88 = rational_func(M=2, N=3)
funccc = rational_func(M=2, N=3)
func08 = rational_func(M=2, N=3)
fit33 = fit_routine(func33, Qgev, pi33_highq, 10)
fit88 = fit_routine(func88, Qgev, pi88_highq, 10)
fitcc = fit_routine(funccc, Qgev, picc_highq, 10)
fit08 = fit_routine(func08, Qgev, pi08_highq, 10)

## plot 
ax = gca()
ax.tick_params(right=true)

# isovector 
errorbar(Qgev, value.(pi33_highq), err.(pi33_highq), fmt="d", color="dodgerblue", mfc="none", ms=5)
xarr = range(0,9,length=100)
yarr33 = func33(xarr, fit33.param); uwerr.(yarr33)
fill_between(xarr, value.(yarr33).-err.(yarr33), value.(yarr33).+err.(yarr33), alpha=0.5, color="dodgerblue", label=L"$I=1$")
plot(xarr, value.(yarr33), color="dodgerblue")

# isoscalar 
errorbar(Qgev, value.(pi88_highq), err.(pi88_highq), fmt="d", color="forestgreen", mfc="none", ms=5)
xarr = range(0,9,length=100)
yarr88 = func88(xarr, fit88.param); uwerr.(yarr88)
fill_between(xarr, value.(yarr88).-err.(yarr88), value.(yarr88).+err.(yarr88), alpha=0.5, color="forestgreen", label=L"$I=0$")
plot(xarr, value.(yarr88), color="forestgreen")

# cc connected
errorbar(Qgev, value.(picc_highq), err.(picc_highq), fmt="d", color="darkorange", mfc="none", ms=5)
xarr = range(0,9,length=100)
yarrcc = funccc(xarr, fitcc.param); uwerr.(yarrcc)
fill_between(xarr, value.(yarrcc).-err.(yarrcc), value.(yarrcc).+err.(yarrcc), alpha=0.5, label=L"$\mathrm{charm}$", color="darkorange")
plot(xarr, value.(yarrcc), color="darkorange")

# 08
errorbar(Qgev, value.(pi08_highq), err.(pi08_highq), fmt="d", color="plum", mfc="none", ms=5)
xarr = range(0,9,length=100)
yarr08 = func08(xarr, fit08.param); uwerr.(yarr08)
fill_between(xarr, value.(yarr08).-err.(yarr08), value.(yarr08).+err.(yarr08), alpha=0.5, label=L"$08$", color="plum")
plot(xarr, value.(yarr08), color="plum")

xlabel(L"$Q^2\ [\mathrm{GeV}^2]$")
ylabel(L"$\Pi(-Q^2) -\Pi(-Q^2/4) $")
legend()
tight_layout()
display(gcf())
savefig(joinpath(path_plot, "hvp_running_highq.pdf"))
close("all")


######################################
## LOW Q
#####################################
## read phys res data LOW Q
# reading systematic errors

syst_lowq = read_systematics(path_phys_res_lowq)

# 33
pi33_lowq_sub = read_phys_res(path_phys_res_lowq, "PI33_physRes.bdio") 
pi33_lowq_sub = pi33_lowq_sub .+ [uwreal([0.0, syst_lowq["pi33"][k]], "syst 33 lowq") for k in eachindex(Qgev)]
b33qqm_lowq = [(q/(4*Qmgev)) * log(2) / (4pi^2) for q in Qgev ./4]
pi33_lowq = pi33_lowq_sub .+ b33qqm_lowq

# 88
pi88_lowq_sub = read_phys_res(path_phys_res_lowq, "PI88_physRes.bdio")
pi88_lowq_sub = pi88_lowq_sub .+ [uwreal([0.0, syst_lowq["pi88"][k]], "syst 88 lowq") for k in eachindex(Qgev)]
pi88_lowq = pi88_lowq_sub .+ b33qqm_lowq

# cc connected
picc_lowq_sub = read_phys_res(path_phys_res_lowq, "PIcc_conn_physRes.bdio")
picc_lowq_sub = picc_lowq_sub .+ [uwreal([0.0, syst_lowq["piccconn"][k]], "syst cc lowq") for k in eachindex(Qgev)]

bcc_lowq = read_phys_res(path_phys_res_lowq, "PI_blc_physRes.bdio")
bcc_lowq = bcc_lowq .+ [uwreal([0.0, syst_lowq["piblc"][k]], "syst blc lowq") for k in eachindex(Qgev)]

picc_lowq = picc_lowq_sub  .+ 2 .* charge_factor["cc"] .* b33qqm_lowq #.+ bcc_lowq

# 08
pi08_lowq = read_phys_res(path_phys_res_lowq, "PI08_physRes.bdio")
pi08_lowq = pi08_lowq .+ [uwreal([0.0, syst_lowq["pi80"][k]], "syst 80 lowq") for k in eachindex(Qgev)]

# remove charge factor
pi88_lowq ./= charge_factor["88"]
picc_lowq ./= charge_factor["cc"]
pi08_lowq ./= charge_factor["08"]
## fit
func33 = rational_func(M=2, N=3)
func88 = rational_func(M=2, N=3)
funccc = rational_func(M=2, N=3)
func08 = rational_func(M=2, N=3)
fit33 = fit_routine(func33, Qgev ./ 4, pi33_lowq, 10)
fit88 = fit_routine(func88, Qgev ./ 4, pi88_lowq, 10)
fitcc = fit_routine(funccc, Qgev ./ 4, picc_lowq, 10)
fit08 = fit_routine(func08, Qgev ./ 4, pi08_lowq, 10)

## plot
ax = gca()
ax.tick_params(right=true)
# isovector 
errorbar(Qgev ./ 4, value.(pi33_lowq), err.(pi33_lowq), fmt="d", color="dodgerblue", mfc="none", ms=5)
xarr = range(0,9/4,length=100)
yarr33 = func33(xarr, fit33.param); uwerr.(yarr33)
fill_between(xarr, value.(yarr33).-err.(yarr33), value.(yarr33).+err.(yarr33), alpha=0.5, color="dodgerblue", label=L"$I=1$")
plot(xarr, value.(yarr33), color="dodgerblue")

# isoscalar 
errorbar(Qgev ./ 4, value.(pi88_lowq), err.(pi88_lowq), fmt="d", color="forestgreen", mfc="none", ms=5)
xarr = range(0,9/4,length=100)
yarr88 = func88(xarr, fit88.param); uwerr.(yarr88)
fill_between(xarr, value.(yarr88).-err.(yarr88), value.(yarr88).+err.(yarr88), alpha=0.5, color="forestgreen", label=L"$I=0$")
plot(xarr, value.(yarr88), color="forestgreen")

# cc connected
errorbar(Qgev ./ 4, value.(picc_lowq), err.(picc_lowq), fmt="d", color="darkorange", mfc="none", ms=5)
xarr = range(0,9/4,length=100)
yarrcc = funccc(xarr, fitcc.param); uwerr.(yarrcc)
fill_between(xarr, value.(yarrcc).-err.(yarrcc), value.(yarrcc).+err.(yarrcc), alpha=0.5, label=L"$\mathrm{charm}$", color="darkorange")
plot(xarr, value.(yarrcc), color="darkorange")

# 08
errorbar(Qgev ./ 4 , value.(pi08_lowq), err.(pi08_lowq), fmt="d", color="plum", mfc="none", ms=5)
xarr = range(0,9/4,length=100)
yarr08 = func08(xarr, fit08.param); uwerr.(yarr08)
fill_between(xarr, value.(yarr08).-err.(yarr08), value.(yarr08).+err.(yarr08), alpha=0.5, label=L"$08$", color="plum")
plot(xarr, value.(yarr08), color="plum")

xlabel(L"$Q^2\ [\mathrm{GeV}^2]$")
ylabel(L"$\Pi(-Q^2/4) -\Pi(0) $")
legend()
tight_layout()
display(gcf())
savefig(joinpath(path_plot, "hvp_running_lowq.pdf"))
close("all")


######################
## TOTAL Q
#####################
syst_fullq = read_systematics(path_phys_res_fullq)

pi33_tot = pi33_highq .+ pi33_lowq
pi88_tot = pi88_highq .+ pi88_lowq
pi08_tot = pi08_highq .+ pi08_lowq

# picc_tot = picc_highq .+ picc_lowq
picc_tot = read_phys_res(path_phys_res_fullq, "PIcc_conn_physRes.bdio")
picc_tot = picc_tot .+ [uwreal([0.0, syst_fullq["piccconn"][k]], "syst cc fullq") for k in eachindex(Qgev)]
picc_tot ./= charge_factor["cc"] # remove charge factor 

wpmm = Dict{String, Vector{Float64}}()
wpmm["E300"] = [-1., 2., -1., -1.]
func33 = rational_func(M=2, N=3)
func88 = rational_func(M=2, N=3)
funccc = rational_func(M=2, N=3)
func08 = rational_func(M=2, N=3)
fit33 = fit_routine(func33, Qgev, pi33_tot, 10)
fit88 = fit_routine(func88, Qgev, pi88_tot, 10)
fitcc = fit_routine(funccc, Qgev, picc_tot, 10, wpm=wpmm)
fit08 = fit_routine(func08, Qgev, pi08_tot, 10)

##

## figure
# isovector 
ax = gca()
ax.tick_params(right=true)

errorbar(Qgev, value.(pi33_tot), err.(pi33_tot), fmt="d", color="dodgerblue", mfc="none", ms=5)
xarr = range(0,9,length=100)
yarr33 = func33(xarr, fit33.param); uwerr.(yarr33)
fill_between(xarr, value.(yarr33).-err.(yarr33), value.(yarr33).+err.(yarr33), alpha=0.5, color="dodgerblue", label=L"$I=1$")
plot(xarr, value.(yarr33), color="dodgerblue")

# isoscalar 
errorbar(Qgev, value.(pi88_tot), err.(pi88_tot), fmt="d", color="forestgreen", mfc="none", ms=5)
xarr = range(0,9,length=100)
yarr88 = func88(xarr, fit88.param); uwerr.(yarr88)
fill_between(xarr, value.(yarr88).-err.(yarr88), value.(yarr88).+err.(yarr88), alpha=0.5, color="forestgreen", label=L"$I=0$")
plot(xarr, value.(yarr88), color="forestgreen")

# cc connected
errorbar(Qgev, value.(picc_tot), err.(picc_tot), fmt="d", color="darkorange", mfc="none", ms=5)
xarr = range(0,9,length=100)
yarrcc = funccc(xarr, fitcc.param); uwerr.(yarrcc)
fill_between(xarr, value.(yarrcc).-err.(yarrcc), value.(yarrcc).+err.(yarrcc), alpha=0.5, label=L"$\mathrm{charm}$", color="darkorange")
plot(xarr, value.(yarrcc), color="darkorange")

# 08
errorbar(Qgev, value.(pi08_tot), err.(pi08_tot), fmt="d", color="plum", mfc="none", ms=5)
xarr = range(0,9,length=100)
yarr08 = func08(xarr, fit08.param); uwerr.(yarr08)
fill_between(xarr, value.(yarr08).-err.(yarr08), value.(yarr08).+err.(yarr08), alpha=0.5, label=L"$08$", color="plum")
plot(xarr, value.(yarr08), color="plum")

xlabel(L"$Q^2\ [\mathrm{GeV}^2]$")
ylabel(L"$\mathit{\bar\Pi}(-Q^2)$")
legend()
tight_layout()
display(gcf())
savefig(joinpath(path_plot, "hvp_running_totq.pdf"))
close("all")

## comparison of charm plot

picc_tot = picc_highq .+ picc_lowq
picc_tot_fullq = read_phys_res(path_phys_res_fullq, "PIcc_conn_physRes.bdio")
picc_tot_fullq = picc_tot_fullq .+ [uwreal([0.0, syst_fullq["piccconn"][k]], "syst cc fullq") for k in eachindex(Qgev)]
picc_tot_fullq ./= charge_factor["cc"] # remove charge factor 

funccc = rational_func(M=2, N=3)
funccc_fullq = rational_func(M=2, N=3)

fitcc = fit_routine(funccc, Qgev, picc_tot, 10, wpm=wpmm)
fitcc_fullq = fit_routine(funccc_fullq, Qgev, picc_tot_fullq, 10, wpm=wpmm)

ax = gca()
ax.tick_params(right=true)

# cc connected
errorbar(Qgev, value.(picc_tot), err.(picc_tot), fmt="d", color="darkorange", mfc="none", ms=5)
xarr = range(0,9,length=100)
yarrcc = funccc(xarr, fitcc.param); uwerr.(yarrcc)
fill_between(xarr, value.(yarrcc).-err.(yarrcc), value.(yarrcc).+err.(yarrcc), alpha=0.5, label=L"$\mathrm{SD\ strategy}$", color="darkorange")
plot(xarr, value.(yarrcc), color="darkorange")

# cc connected full q
errorbar(Qgev, value.(picc_tot_fullq), err.(picc_tot_fullq), fmt="d", color="red", mfc="none", ms=5)
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
