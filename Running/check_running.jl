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

# Qgev = [1.0, 5.0, 9.0 ] # GeV^2
const Qgev = [0.05, 0.1, 0.4, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0] # Q^2
NMOM = length(Qgev)
Qmgev =  36.0  # GeV^2


#####################
## HIGH Q
#####################
## read phys res data HIGH Q
# reading systematic errors
f = readdlm(joinpath(path_phys_res_highq, "systematics.txt"))
delim = findall(x-> x=="#", f)
dict_syst = Dict()
for k in delim
    idx = k.I[1]
    name = join(f[idx,2:end])
    dict_syst[name] = Float64.(f[idx+1:idx+12, 2])
end

# 33
pi33_SD = read_phys_res(path_phys_res_highq, "PI33_SD_physRes.bdio")
pi33_ILD = read_phys_res(path_phys_res_highq, "PI33_ILD_physRes.bdio")

b33qqm_highq = [(q/(4*Qmgev)) * log(2) / (4pi^2) for q in Qgev]
pi33_highq = pi33_ILD + pi33_SD + b33qqm_highq

# 88
pi3388dlt_highq = Vector{uwreal}(undef, 0)
[push!(pi3388dlt_highq, read_phys_res(path_phys_res_highq, "PI3388_physRes.bdio")[k]) for k in 1:NMOM]
pi88_highq = pi33_highq/3 - pi3388dlt_highq 

# cc connected
picc_highq_sub = Vector{uwreal}(undef, 0)
[push!(picc_highq_sub, read_phys_res(path_phys_res_highq, "PIcc_conn_physRes.bdio")[k]) for k in 1:NMOM]
bcc_highq = Vector{uwreal}(undef, 0)
[push!(bcc_highq, read_phys_res(path_phys_res_highq, "PI_blc_physRes.bdio")[k]) for k in 1:NMOM]
picc_highq = picc_highq_sub .+ bcc_highq  .+ 2 * 4/9 * b33qqm_highq 

# 08
pi08_highq = Vector{uwreal}(undef, 0)
[push!(pi08_highq, read_phys_res(path_phys_res_highq, "PI08_physRes.bdio")[k]) for k in 1:NMOM]


function rational_func(; M=3, N=3)
    num = (x,p) -> sum(p[k] .* x.^k for k in 1:M)
    den = (x,p) -> 1. .+ sum(p[k+M] .* x.^k for k in 1:N)

    num1 = (x,p) -> sum(p[k+M+N] .* (x ./4).^k for k in 1:M)
    den1 = (x,p) -> 1. .+ sum(p[k+M+M+N] .* (x ./4).^k for k in 1:N)

    func(x,p) = num(x,p) ./ den(x,p) .- num1(x,p)./ den1(x,p)
    return func
end


## test fit

func33 = rational_func(M=2, N=3)
func88 = rational_func(M=2, N=3)
funccc = rational_func(M=2, N=3)
func08 = rational_func(M=2, N=3)
wpmm = Dict{String, Vector{Float64}}()
wpmm["N451"] = [1., -1., -1., -1.]
fit33 = fit_routine(func33, Qgev, pi33_highq, 10 , wpm=wpmm)
fit88 = fit_routine(func88, Qgev, pi88_highq, 10, wpm=wpmm)
fitcc = fit_routine(funccc, Qgev, picc_highq, 10, wpm=wpmm)
fit08 = fit_routine(func08, Qgev, pi08_highq, 10, wpm=wpmm)

## figure
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


#total
# yarr = yarr33 .+ yarrcc .+ yarr88 .+ yarr08; uwerr.(yarr)
# fill_between(xarr, value.(yarr).-err.(yarr), value.(yarr).+err.(yarr), alpha=0.5, label=L"$\mathrm{total}$", color="black")
# plot(xarr, value.(yarr), color="black")


xlabel(L"$Q^2\ [\mathrm{GeV}^2]$")
ylabel(L"$\Pi(-Q^2) -\Pi(-Q^2/4) $")
legend()
tight_layout()
display(gcf())
savefig(joinpath(path_plot, "hvp_running_highq.pdf"))
close("all")

#########################
## LOW Q
#########################
# 33
pi33_lowq_sub = read_phys_res(path_phys_res_lowq, "PI33_physRes.bdio")

b33qqm_lowq = [(q/(4*Qmgev)) * log(2) / (4pi^2) for q in Qgev ./4]
pi33_lowq = pi33_lowq_sub  + b33qqm_lowq

# 88
pi88_aux = Vector{uwreal}(undef, 0)
[push!(pi88_aux, read_phys_res(path_phys_res_lowq, "PI88_physRes.bdio")[k]) for k in 1:NMOM]
pi88_lowq = pi88_aux .+ b33qqm_lowq

# cc connected
picc_lowq_sub = Vector{uwreal}(undef, 0)
[push!(picc_lowq_sub, read_phys_res(path_phys_res_lowq, "PIcc_conn_physRes.bdio")[k]) for k in 1:NMOM]
bcc_lowq = Vector{uwreal}(undef, 0)
# [push!(bcc_lowq, read_phys_res(path_phys_res_lowq, "PI_blc_physRes.bdio")[k]) for k in 1:NMOM]
picc_lowq = picc_lowq_sub .+ 2*4/9*b33qqm_lowq

# 08
pi08_lowq = Vector{uwreal}(undef, 0)
[push!(pi08_lowq, read_phys_res(path_phys_res_lowq, "PI08_physRes.bdio")[k]) for k in 1:NMOM]


## test fit

func33 = rational_func(M=2, N=3)
func88 = rational_func(M=2, N=3)
funccc = rational_func(M=2, N=3)
func08 = rational_func(M=2, N=3)
wpmm = Dict{String, Vector{Float64}}()
wpmm["N451"] = [1., -1., -1., -1.]
fit33 = fit_routine(func33, Qgev./4, pi33_lowq, 10 , wpm=wpmm)
fit88 = fit_routine(func88, Qgev./4, pi88_lowq, 10, wpm=wpmm)
fitcc = fit_routine(funccc, Qgev./4, picc_lowq_sub, 10, wpm=wpmm)
fit08 = fit_routine(func08, Qgev./4, pi08_lowq, 10, wpm=wpmm)

## figure
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
errorbar(Qgev ./ 4, value.(picc_lowq_sub), err.(picc_lowq_sub), fmt="d", color="darkorange", mfc="none", ms=5)
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


#total
# yarr = yarr33 .+ yarrcc .+ yarr88 .+ yarr08; uwerr.(yarr)
# fill_between(xarr, value.(yarr).-err.(yarr), value.(yarr).+err.(yarr), alpha=0.5, label=L"$\mathrm{total}$", color="black")
# plot(xarr, value.(yarr), color="black")


xlabel(L"$Q^2\ [\mathrm{GeV}^2]$")
ylabel(L"$\Pi(-Q^2/4) -\Pi(0) $")
legend()
tight_layout()
display(gcf())
savefig(joinpath(path_plot, "hvp_running_lowq.pdf"))
close("all")


###############
## TOTAL Q
###############
pi33_tot = pi33_highq .+ pi33_lowq
pi88_tot = (pi88_highq .+ pi88_lowq) * 3
# picc_tot = (picc_highq .+ picc_lowq) *9/4
picc_tot = Vector{uwreal}(undef, 0)
[push!(picc_tot, read_phys_res(path_phys_res_fullq, "PIcc_conn_physRes.bdio")[k]) for k in 1:NMOM]
picc_tot .*= (9 / 4)
pi08_tot = (pi08_highq .+ pi08_lowq) * (6*sqrt(3))


func33 = rational_func(M=2, N=3)
func88 = rational_func(M=2, N=3)
funccc = rational_func(M=2, N=3)
func08 = rational_func(M=2, N=3)
wpmm = Dict{String, Vector{Float64}}()
wpmm["N451"] = [1., -1., -1., -1.]
fit33 = fit_routine(func33, Qgev, pi33_tot, 10 , wpm=wpmm)
fit88 = fit_routine(func88, Qgev, pi88_tot, 10, wpm=wpmm)
fitcc = fit_routine(funccc, Qgev, picc_tot, 10, wpm=wpmm)
fit08 = fit_routine(func08, Qgev, pi08_tot, 10, wpm=wpmm)


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


#total
# yarr = yarr33 .+ yarrcc .+ yarr88 .+ yarr08; uwerr.(yarr)
# fill_between(xarr, value.(yarr).-err.(yarr), value.(yarr).+err.(yarr), alpha=0.5, label=L"$\mathrm{total}$", color="black")
# plot(xarr, value.(yarr), color="black")

xlabel(L"$Q^2\ [\mathrm{GeV}^2]$")
ylabel(L"$\mathit{\bar\Pi}(-Q^2)$")
legend()
tight_layout()
display(gcf())
savefig(joinpath(path_plot, "hvp_running_totq.pdf"))
close("all")



#################################
## RUNNING OF ALPHA
################################
alpha_lo = 1 / 137.035999084
alpha_33 = 4*pi* alpha_lo * pi33_tot; uwerr.(alpha_33)
alpha_88 = 4*pi* alpha_lo * pi88_tot / 3; uwerr.(alpha_88)
alpha_cc = 4*pi* alpha_lo * picc_tot * 4 / 9 ; uwerr.(alpha_cc)

alpha_tot = alpha_33 .+ alpha_88 .+ alpha_cc

##

func33 = rational_func(M=2, N=3)
func88 = rational_func(M=2, N=3)
funccc = rational_func(M=2, N=3)
# func08 = rational_func(M=2, N=3)
wpmm = Dict{String, Vector{Float64}}()
wpmm["N451"] = [1., -1., -1., -1.]
fit33 = fit_routine(func33, Qgev, alpha_33, 10 , wpm=wpmm)
fit88 = fit_routine(func88, Qgev, alpha_88, 10, wpm=wpmm)
fitcc = fit_routine(funccc, Qgev, alpha_cc, 10, wpm=wpmm)
# fit08 = fit_routine(func08, Qgev, pi08, 10, wpm=wpmm)

## plot
# isovector 
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
savefig(joinpath(path_plot, "alpha_runing.pdf"))
close("all")



##

ax = gca()
xx = collect(1:0.01:2pi)
plot(sin.(xx))
ax.tick_params(right=true)
# loc, labels= yticks()
# yticks(loc, labels)
# ax1.set_yticks(loc, labels)
display(gcf())
close("all")