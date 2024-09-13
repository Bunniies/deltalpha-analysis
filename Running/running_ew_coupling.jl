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

picc_highq = picc_highq_sub  .+ 2 * charge_factor["cc"] * b33qqm_highq #.+ bcc_highq

# 08
pi08_highq = read_phys_res(path_phys_res_highq, "PI08_physRes.bdio")
pi08_highq = pi08_highq .+ [uwreal([0.0, syst_highq["pi80"][k]], "syst 80") for k in eachindex(Qgev)]


## runing of alpha HIGH Q

alpha_lo = 4*pi / 137.035999084
alpha_33_highq = alpha_lo .* pi33_highq
alpha_88_highq = alpha_lo .* pi88_highq
alpha_cc_highq = alpha_lo .* picc_highq

func33 = rational_func(M=2, N=3)
func88 = rational_func(M=2, N=3)
funccc = rational_func(M=2, N=3)

wpmm = Dict{String, Vector{Float64}}()
wpmm["N451"] = [1., -1., -1., -1.]
fit33 = fit_routine(func33, Qgev, alpha_33_highq, 10)#, wpm=wpmm)
fit88 = fit_routine(func88, Qgev, alpha_88_highq, 10)#, wpm=wpmm)
fitcc = fit_routine(funccc, Qgev, alpha_cc_highq, 10)#, wpm=wpmm)


## plot
# isovector 
ax = gca()
ax.tick_params(right=true)
errorbar(Qgev, value.(alpha_33_highq), err.(alpha_33_highq), fmt="d", color="dodgerblue", mfc="none", ms=5)
xarr = range(0,9,length=100)
yarr33 = func33(xarr, fit33.param); uwerr.(yarr33)
fill_between(xarr, value.(yarr33).-err.(yarr33), value.(yarr33).+err.(yarr33), alpha=0.5, color="dodgerblue", label=L"$I=1$")
plot(xarr, value.(yarr33), color="dodgerblue")

# isoscalar 
errorbar(Qgev, value.(alpha_88_highq), err.(alpha_88_highq), fmt="d", color="forestgreen", mfc="none", ms=5)
xarr = range(0,9,length=100)
yarr88 = func88(xarr, fit88.param); uwerr.(yarr88)
fill_between(xarr, value.(yarr88).-err.(yarr88), value.(yarr88).+err.(yarr88), alpha=0.5, color="forestgreen", label=L"$I=0$")
plot(xarr, value.(yarr88), color="forestgreen")

# cc connected
errorbar(Qgev , value.(alpha_cc_highq), err.(alpha_cc_highq), fmt="d", color="darkorange", mfc="none", ms=5)
xarr = range(0,9,length=100)
yarrcc = funccc(xarr, fitcc.param); uwerr.(yarrcc)
fill_between(xarr, value.(yarrcc).-err.(yarrcc), value.(yarrcc).+err.(yarrcc), alpha=0.5, label=L"$\mathrm{charm}$", color="darkorange")
plot(xarr, value.(yarrcc), color="darkorange")

# total
yarr = yarr33 .+ yarrcc .+ yarr88 ; uwerr.(yarr)
fill_between(xarr, value.(yarr).-err.(yarr), value.(yarr).+err.(yarr), alpha=0.5, label=L"$\mathrm{total}$", color="black")
plot(xarr, value.(yarr), color="black")

xlabel(L"$Q^2\ [\mathrm{GeV}^2]$")
ylabel(L"$\Delta \alpha(-Q^2) - \Delta \alpha(-Q^2/4)$")
#ylabel(L"$\Delta \alpha_{\mathrm{had}}(-Q^2)$")
legend()
tight_layout()
display(gcf())
savefig(joinpath(path_plot, "alpha_runing_highq.pdf"))
close("all")

## SIN2 OMEGA running  HIGH Q

# sinW_tompLim = (0.5 - uwreal([0.22337, 0.0001], "sin2W"))
sinW_tompLim = (0.5 - uwreal([0.23857, 0.00005], "sin2W"))

norm_sinW = -alpha_lo / sinW_tompLim

sinW_33_highq = norm_sinW * (0.5 - sinW_tompLim) * pi33_highq
sinW_88_highq = norm_sinW * (0.5 - sinW_tompLim) * pi88_highq
sinW_cc_highq = norm_sinW * ((0.5 - sinW_tompLim) * picc_highq - 1/18*picc_highq)
sinW_08_highq = - norm_sinW * pi08_highq

func33 = rational_func(M=2, N=3)
func88 = rational_func(M=2, N=3)
funccc = rational_func(M=2, N=3)
func08 = rational_func(M=2, N=3)
wpmm["E300"] = [-1., 1., -1., -1.]

fit33 = fit_routine(func33, Qgev, sinW_33_highq, 10, wpm=wpmm)
fit88 = fit_routine(func88, Qgev, sinW_88_highq, 10, wpm=wpmm)
fitcc = fit_routine(funccc, Qgev, sinW_cc_highq, 10, wpm=wpmm)
fit08 = fit_routine(func08, Qgev, sinW_08_highq, 10, wpm=wpmm)

## plot sinW high q
# isovector 
ax = gca()
ax.tick_params(right=true)
errorbar(Qgev, value.(sinW_33_highq), err.(sinW_33_highq), fmt="d", color="dodgerblue", mfc="none", ms=5)
xarr = range(0,9,length=100)
yarr33 = func33(xarr, fit33.param); uwerr.(yarr33)
fill_between(xarr, value.(yarr33).-err.(yarr33), value.(yarr33).+err.(yarr33), alpha=0.5, color="dodgerblue", label=L"$I=1$")
plot(xarr, value.(yarr33), color="dodgerblue")

# isoscalar 
errorbar(Qgev, value.(sinW_88_highq), err.(sinW_88_highq), fmt="d", color="forestgreen", mfc="none", ms=5)
xarr = range(0,9,length=100)
yarr88 = func88(xarr, fit88.param); uwerr.(yarr88)
fill_between(xarr, value.(yarr88).-err.(yarr88), value.(yarr88).+err.(yarr88), alpha=0.5, color="forestgreen", label=L"$I=0$")
plot(xarr, value.(yarr88), color="forestgreen")

# cc connected
errorbar(Qgev , value.(sinW_cc_highq), err.(sinW_cc_highq), fmt="d", color="darkorange", mfc="none", ms=5)
xarr = range(0,9,length=100)
yarrcc = funccc(xarr, fitcc.param); uwerr.(yarrcc)
fill_between(xarr, value.(yarrcc).-err.(yarrcc), value.(yarrcc).+err.(yarrcc), alpha=0.5, label=L"$\mathrm{charm}$", color="darkorange")
plot(xarr, value.(yarrcc), color="darkorange")

# 08
errorbar(Qgev, value.(sinW_08_highq), err.(sinW_08_highq), fmt="d", color="plum", mfc="none", ms=5)
xarr = range(0,9,length=100)
yarr08 = func08(xarr, fit08.param); uwerr.(yarr08)
fill_between(xarr, value.(yarr08).-err.(yarr08), value.(yarr08).+err.(yarr08), alpha=0.5, label=L"$08$", color="plum")
plot(xarr, value.(yarr08), color="plum")

# total
yarr = yarr33 .+ yarrcc .+ yarr88 .+ yarr08 ; uwerr.(yarr)
fill_between(xarr, value.(yarr).-err.(yarr), value.(yarr).+err.(yarr), alpha=0.5, label=L"$\mathrm{total}$", color="black")
plot(xarr, value.(yarr), color="black")

axhline(ls="-.", color="black", lw=0.8)

xlabel(L"$Q^2\ [\mathrm{GeV}^2]$")
ylabel(L"$\Delta\sin^2\theta_W(-Q^2) - \Delta\sin^2\theta_W(-Q^2/4)$", fontsize=16)
# ylabel(L"$\Delta_{\mathrm{had}}\sin^2\theta_W(-Q^2)$")
legend()
tight_layout()
display(gcf())
savefig(joinpath(path_plot, "sinW_runing_highq.pdf"))
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


## running of ALPHA LOW Q

alpha_33_lowq = alpha_lo .* pi33_lowq
alpha_88_lowq = alpha_lo .* pi88_lowq
alpha_cc_lowq = alpha_lo .* picc_lowq

fit33 = fit_routine(func33, Qgev ./4, alpha_33_lowq, 10)#, wpm=wpmm)
fit88 = fit_routine(func88, Qgev ./4, alpha_88_lowq, 10)#, wpm=wpmm)
fitcc = fit_routine(funccc, Qgev ./4, alpha_cc_lowq, 10)#, wpm=wpmm)

## plot
# isovector 
ax = gca()
ax.tick_params(right=true)
errorbar(Qgev ./4 , value.(alpha_33_lowq), err.(alpha_33_lowq), fmt="d", color="dodgerblue", mfc="none", ms=5)
xarr = range(0,9/4,length=100)
yarr33 = func33(xarr, fit33.param); uwerr.(yarr33)
fill_between(xarr, value.(yarr33).-err.(yarr33), value.(yarr33).+err.(yarr33), alpha=0.5, color="dodgerblue", label=L"$I=1$")
plot(xarr, value.(yarr33), color="dodgerblue")

# isoscalar 
errorbar(Qgev ./4, value.(alpha_88_lowq), err.(alpha_88_lowq), fmt="d", color="forestgreen", mfc="none", ms=5)
xarr = range(0,9/4,length=100)
yarr88 = func88(xarr, fit88.param); uwerr.(yarr88)
fill_between(xarr, value.(yarr88).-err.(yarr88), value.(yarr88).+err.(yarr88), alpha=0.5, color="forestgreen", label=L"$I=0$")
plot(xarr, value.(yarr88), color="forestgreen")

# cc connected
errorbar(Qgev ./4 , value.(alpha_cc_lowq), err.(alpha_cc_lowq), fmt="d", color="darkorange", mfc="none", ms=5)
xarr = range(0,9/4,length=100)
yarrcc = funccc(xarr, fitcc.param); uwerr.(yarrcc)
fill_between(xarr, value.(yarrcc).-err.(yarrcc), value.(yarrcc).+err.(yarrcc), alpha=0.5, label=L"$\mathrm{charm}$", color="darkorange")
plot(xarr, value.(yarrcc), color="darkorange")

# total
yarr = yarr33 .+ yarrcc .+ yarr88 ; uwerr.(yarr)
fill_between(xarr, value.(yarr).-err.(yarr), value.(yarr).+err.(yarr), alpha=0.5, label=L"$\mathrm{total}$", color="black")
plot(xarr, value.(yarr), color="black")

xlabel(L"$Q^2/4\ [\mathrm{GeV}^2]$")
ylabel(L"$\Delta \alpha_{\mathrm{had}}(-Q^2/4) - \Delta \alpha_{\mathrm{had}}(0)$")
legend()
tight_layout()
display(gcf())
savefig(joinpath(path_plot, "alpha_runing_lowq.pdf"))
close("all")

## SIN2 OMEGA running  LOW Q

sinW_33_lowq = norm_sinW * (0.5 - sinW_tompLim) * pi33_lowq
sinW_88_lowq = norm_sinW * (0.5 - sinW_tompLim) * pi88_lowq
sinW_cc_lowq = norm_sinW * ((0.5 - sinW_tompLim) * picc_lowq - 1/18*picc_lowq)
sinW_08_lowq = - norm_sinW * pi08_lowq

func33 = rational_func(M=2, N=3)
func88 = rational_func(M=2, N=3)
funccc = rational_func(M=2, N=3)
func08 = rational_func(M=2, N=3)
fit33 = fit_routine(func33, Qgev ./4, sinW_33_lowq, 10, wpm=wpmm)
fit88 = fit_routine(func88, Qgev ./4, sinW_88_lowq, 10, wpm=wpmm)
fitcc = fit_routine(funccc, Qgev ./4, sinW_cc_lowq, 10, wpm=wpmm)
fit08 = fit_routine(func08, Qgev ./4, sinW_08_lowq, 10, wpm=wpmm)

# plot sinW low q
# isovector 
ax = gca()
ax.tick_params(right=true)
errorbar(Qgev ./4, value.(sinW_33_lowq), err.(sinW_33_lowq), fmt="d", color="dodgerblue", mfc="none", ms=5)
xarr = range(0,9/4,length=100)
yarr33 = func33(xarr, fit33.param); uwerr.(yarr33)
fill_between(xarr, value.(yarr33).-err.(yarr33), value.(yarr33).+err.(yarr33), alpha=0.5, color="dodgerblue", label=L"$I=1$")
plot(xarr, value.(yarr33), color="dodgerblue")

# isoscalar 
errorbar(Qgev ./4, value.(sinW_88_lowq), err.(sinW_88_lowq), fmt="d", color="forestgreen", mfc="none", ms=5)
xarr = range(0,9/4,length=100)
yarr88 = func88(xarr, fit88.param); uwerr.(yarr88)
fill_between(xarr, value.(yarr88).-err.(yarr88), value.(yarr88).+err.(yarr88), alpha=0.5, color="forestgreen", label=L"$I=0$")
plot(xarr, value.(yarr88), color="forestgreen")

# cc connected
errorbar(Qgev ./4 , value.(sinW_cc_lowq), err.(sinW_cc_lowq), fmt="d", color="darkorange", mfc="none", ms=5)
xarr = range(0,9/4,length=100)
yarrcc = funccc(xarr, fitcc.param); uwerr.(yarrcc)
fill_between(xarr, value.(yarrcc).-err.(yarrcc), value.(yarrcc).+err.(yarrcc), alpha=0.5, label=L"$\mathrm{charm}$", color="darkorange")
plot(xarr, value.(yarrcc), color="darkorange")

# 08
errorbar(Qgev ./4, value.(sinW_08_lowq), err.(sinW_08_lowq), fmt="d", color="plum", mfc="none", ms=5)
xarr = range(0,9/4,length=100)
yarr08 = func08(xarr, fit08.param); uwerr.(yarr08)
fill_between(xarr, value.(yarr08).-err.(yarr08), value.(yarr08).+err.(yarr08), alpha=0.5, label=L"$08$", color="plum")
plot(xarr, value.(yarr08), color="plum")

# total
yarr = yarr33 .+ yarrcc .+ yarr88 .+ yarr08 ; uwerr.(yarr)
fill_between(xarr, value.(yarr).-err.(yarr), value.(yarr).+err.(yarr), alpha=0.5, label=L"$\mathrm{total}$", color="black")
plot(xarr, value.(yarr), color="black")

axhline(ls="-.", color="black", lw=0.8)

xlabel(L"$Q^2/4\ [\mathrm{GeV}^2]$")
ylabel(L"$\Delta_{\mathrm{had}}\sin^2\theta_W(-Q^2)$")
legend()
tight_layout()
display(gcf())
savefig(joinpath(path_plot, "sinW_runing_lowq.pdf"))
close("all")

#########################
## RUNNING ALPHA FULL Q
#########################

alpha_33 = alpha_33_highq .+ alpha_33_lowq
alpha_88 = alpha_88_highq .+ alpha_88_lowq
alpha_cc = alpha_cc_highq .+ alpha_cc_lowq

fit33 = fit_routine(func33, Qgev, alpha_33, 10)#, wpm=wpmm)
fit88 = fit_routine(func88, Qgev, alpha_88, 10)#, wpm=wpmm)
fitcc = fit_routine(funccc, Qgev, alpha_cc, 10)#, wpm=wpmm)

ax = gca()
ax.tick_params(right=true)
errorbar(Qgev, value.(alpha_33), err.(alpha_33), fmt="d", color="dodgerblue", mfc="none", ms=5)
xarr = range(0,9,length=100)
yarr33 = func33(xarr, fit33.param); uwerr.(yarr33)
fill_between(xarr, value.(yarr33).-err.(yarr33), value.(yarr33).+err.(yarr33), alpha=0.5, color="dodgerblue", label=L"$I=1$")
plot(xarr, value.(yarr33), color="dodgerblue")

# isoscalar 
errorbar(Qgev, value.(alpha_88), err.(alpha_88), fmt="d", color="forestgreen", mfc="none", ms=5)
xarr = range(0,9,length=100)
yarr88 = func88(xarr, fit88.param); uwerr.(yarr88)
fill_between(xarr, value.(yarr88).-err.(yarr88), value.(yarr88).+err.(yarr88), alpha=0.5, color="forestgreen", label=L"$I=0$")
plot(xarr, value.(yarr88), color="forestgreen")

# cc connected
errorbar(Qgev , value.(alpha_cc), err.(alpha_cc), fmt="d", color="darkorange", mfc="none", ms=5)
xarr = range(0,9,length=100)
yarrcc = funccc(xarr, fitcc.param); uwerr.(yarrcc)
fill_between(xarr, value.(yarrcc).-err.(yarrcc), value.(yarrcc).+err.(yarrcc), alpha=0.5, label=L"$\mathrm{charm}$", color="darkorange")
plot(xarr, value.(yarrcc), color="darkorange")

# total
yarr = yarr33 .+ yarrcc .+ yarr88 ; uwerr.(yarr)
fill_between(xarr, value.(yarr).-err.(yarr), value.(yarr).+err.(yarr), alpha=0.5, label=L"$\mathrm{total}$", color="black")
plot(xarr, value.(yarr), color="black")

xlabel(L"$Q^2\ [\mathrm{GeV}^2]$")
ylabel(L"$\Delta \alpha_{\mathrm{had}}(-Q^2)$")
legend()
tight_layout()
display(gcf())
savefig(joinpath(path_plot, "alpha_runing_fullq.pdf"))
close("all")


#########################
## RUNNING SIN W FULL Q
#########################

sinW_33 = sinW_33_highq .+ sinW_33_lowq
sinW_88 = sinW_88_highq .+ sinW_88_lowq
sinW_cc = sinW_cc_highq .+ sinW_cc_lowq
sinW_08 = sinW_08_highq .+ sinW_08_lowq

func33 = rational_func(M=2, N=3)
func88 = rational_func(M=2, N=3)
funccc = rational_func(M=2, N=3)
func08 = rational_func(M=2, N=3)
fit33 = fit_routine(func33, Qgev, sinW_33, 10, wpm=wpmm)
fit88 = fit_routine(func88, Qgev, sinW_88, 10, wpm=wpmm)
fitcc = fit_routine(funccc, Qgev, sinW_cc, 10, wpm=wpmm)
fit08 = fit_routine(func08, Qgev, sinW_08, 10, wpm=wpmm)

# plot sinW full q
# isovector 
ax = gca()
ax.tick_params(right=true)
errorbar(Qgev, value.(sinW_33), err.(sinW_33), fmt="d", color="dodgerblue", mfc="none", ms=5)
xarr = range(0,9,length=100)
yarr33 = func33(xarr, fit33.param); uwerr.(yarr33)
fill_between(xarr, value.(yarr33).-err.(yarr33), value.(yarr33).+err.(yarr33), alpha=0.5, color="dodgerblue", label=L"$I=1$")
plot(xarr, value.(yarr33), color="dodgerblue")

# isoscalar 
errorbar(Qgev, value.(sinW_88), err.(sinW_88), fmt="d", color="forestgreen", mfc="none", ms=5)
xarr = range(0,9,length=100)
yarr88 = func88(xarr, fit88.param); uwerr.(yarr88)
fill_between(xarr, value.(yarr88).-err.(yarr88), value.(yarr88).+err.(yarr88), alpha=0.5, color="forestgreen", label=L"$I=0$")
plot(xarr, value.(yarr88), color="forestgreen")

# cc connected
errorbar(Qgev , value.(sinW_cc), err.(sinW_cc), fmt="d", color="darkorange", mfc="none", ms=5)
xarr = range(0,9,length=100)
yarrcc = funccc(xarr, fitcc.param); uwerr.(yarrcc)
fill_between(xarr, value.(yarrcc).-err.(yarrcc), value.(yarrcc).+err.(yarrcc), alpha=0.5, label=L"$\mathrm{charm}$", color="darkorange")
plot(xarr, value.(yarrcc), color="darkorange")

# 08
errorbar(Qgev, value.(sinW_08), err.(sinW_08), fmt="d", color="plum", mfc="none", ms=5)
xarr = range(0,9,length=100)
yarr08 = func08(xarr, fit08.param); uwerr.(yarr08)
fill_between(xarr, value.(yarr08).-err.(yarr08), value.(yarr08).+err.(yarr08), alpha=0.5, label=L"$08$", color="plum")
plot(xarr, value.(yarr08), color="plum")

# total
yarr = yarr33 .+ yarrcc .+ yarr88 .+ yarr08 ; uwerr.(yarr)
fill_between(xarr, value.(yarr).-err.(yarr), value.(yarr).+err.(yarr), alpha=0.5, label=L"$\mathrm{total}$", color="black")
plot(xarr, value.(yarr), color="black")

axhline(ls="-.", color="black", lw=0.8)

xlabel(L"$Q^2\ [\mathrm{GeV}^2]$")
ylabel(L"$\Delta_{\mathrm{had}}\sin^2\theta_W(-Q^2)$")
legend()
tight_layout()
display(gcf())
savefig(joinpath(path_plot, "sinW_runing_fullq.pdf"))
close("all")

################################
## check error and contributions
################################
uwerr.(pi33_highq)
uwerr.(pi88_highq)
uwerr.(picc_highq)
# pi08_highq .*= 1e5
uwerr.(pi08_highq)
Qgev # 4, 8, 12

Qgev[12]
details(pi08_highq[12])