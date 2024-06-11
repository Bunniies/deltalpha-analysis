using HVPobs
using ALPHAio, BDIO
using PyPlot, LaTeXStrings
using ADerrors
import ADerrors:err
#================ SET UP PLOT PARAMETERS ==================#

rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
rcParams["text.usetex"] =  true
rcParams["mathtext.fontset"]  = "cm"
rcParams["font.size"] =13
rcParams["axes.labelsize"] =22
rcParams["axes.titlesize"] = 18
plt.rc("text", usetex=true) # set to true if a LaTeX installation is present

path_phys_res = "/Users/alessandroconigli/MyDrive/postdoc-mainz/projects/deltalpha/physical_results/multi_mom/"
path_plot = "/Users/alessandroconigli/MyDrive/postdoc-mainz/projects/deltalpha/plots/running"
include("../utils/IO_BDIO.jl")

# Qgev = [1.0, 5.0, 9.0 ] # GeV^2
const Qgev = [0.05, 0.1, 0.4, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0] # Q^2
NMOM = length(Qgev)
Qmgev =  9.0  # GeV^2

## read phys res data
# 33
pi33_SD = read_phys_res(path_phys_res, "PI33_SD_physRes.bdio")
pi33_ILD = read_phys_res(path_phys_res, "PI33_ILD_physRes.bdio")

bqqm = [(q/(4*Qmgev)) * log(2) / (4pi^2) for q in Qgev]
pi33 = pi33_ILD + pi33_SD + bqqm

# 88
pi3388dlt = Vector{uwreal}(undef, 0)
[push!(pi3388dlt, read_phys_res(path_phys_res, "PI3388_physRes.bdio")[k]) for k in 1:NMOM]
pi88 = pi33/3 - pi3388dlt 

# cc connected
picc = Vector{uwreal}(undef, 0)
[push!(picc, read_phys_res(path_phys_res, "PIcc_conn_physRes.bdio")[k]) for k in 1:NMOM]

# 08
pi08 = Vector{uwreal}(undef, 0)
[push!(pi08, read_phys_res(path_phys_res, "PI08_physRes.bdio")[k]) for k in 1:NMOM]




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
fit33 = fit_routine(func33, Qgev, pi33, 10 , wpm=wpmm)
fit88 = fit_routine(func88, Qgev, pi88, 10, wpm=wpmm)
fitcc = fit_routine(funccc, Qgev, picc, 10, wpm=wpmm)
fit08 = fit_routine(func08, Qgev, pi08, 10, wpm=wpmm)

## figure
# isovector 
errorbar(Qgev, value.(pi33), err.(pi33), fmt="d", color="dodgerblue", mfc="none", ms=5)
xarr = range(0,9,length=100)
yarr33 = func33(xarr, fit33.param); uwerr.(yarr33)
fill_between(xarr, value.(yarr33).-err.(yarr33), value.(yarr33).+err.(yarr33), alpha=0.5, color="dodgerblue", label=L"$I=1$")
plot(xarr, value.(yarr33), color="dodgerblue")

# isoscalar 
errorbar(Qgev, value.(pi88), err.(pi88), fmt="d", color="forestgreen", mfc="none", ms=5)
xarr = range(0,9,length=100)
yarr88 = func88(xarr, fit88.param); uwerr.(yarr88)
fill_between(xarr, value.(yarr88).-err.(yarr88), value.(yarr88).+err.(yarr88), alpha=0.5, color="forestgreen", label=L"$I=0$")
plot(xarr, value.(yarr88), color="forestgreen")

# cc connected
errorbar(Qgev, value.(picc), err.(picc), fmt="d", color="darkorange", mfc="none", ms=5)
xarr = range(0,9,length=100)
yarrcc = funccc(xarr, fitcc.param); uwerr.(yarrcc)
fill_between(xarr, value.(yarrcc).-err.(yarrcc), value.(yarrcc).+err.(yarrcc), alpha=0.5, label=L"$\mathrm{charm}$", color="darkorange")
plot(xarr, value.(yarrcc), color="darkorange")

# 08
errorbar(Qgev, value.(pi08), err.(pi08), fmt="d", color="plum", mfc="none", ms=5)
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
savefig(joinpath(path_plot, "hvp_running.pdf"))
close("all")