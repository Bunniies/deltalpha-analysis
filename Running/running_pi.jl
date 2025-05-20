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
rcParams["font.size"] =20
rcParams["axes.labelsize"] =26
rcParams["axes.titlesize"] = 22
plt.rc("text", usetex=true) # set to true if a LaTeX installation is present

path_phys_res_highq = "/Users/alessandroconigli/MyDrive/postdoc-mainz/projects/deltalpha/physical_results/scale_error_artificial/piQ_Q4/"
path_phys_res_lowq = "/Users/alessandroconigli/MyDrive/postdoc-mainz/projects/deltalpha/physical_results/scale_error_artificial/piQ4_piQ0/"
path_phys_res_fullq = "/Users/alessandroconigli/MyDrive/postdoc-mainz/projects/deltalpha/physical_results/scale_error_artificial/full_q/"
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
# reading systematic errors

# syst_highq = read_systematics(path_phys_res_highq)

# 33
pi33_SD = read_phys_res(path_phys_res_highq, "PI33_SD_physRes.bdio")
pi33_ILD = read_phys_res(path_phys_res_highq, "PI33_ILD_physRes.bdio")

b33qm_pt_5loops = uwreal([1896.11*1e-5, 2.67*1e-5], "b33_qm_pt") # from wolfram code
b33qqm_highq = [(q/(4*Qmgev)) * b33qm_pt_5loops for q in Qgev] # using PT result
# b33qqm_highq = [(q/(4*Qmgev)) * log(2) / (4pi^2) for q in Qgev] # at LO
pi33_highq = (pi33_ILD + pi33_SD + b33qqm_highq ) 
#add_t0_err!(pi33_highq, t0sqrt_ph_err)

# 88
pi3388dlt_highq = read_phys_res(path_phys_res_highq, "PI3388_physRes.bdio")
# add_t0_err!(pi3388dlt_highq, t0sqrt_ph_err)
pi88_highq = pi33_highq * charge_factor["88"] - pi3388dlt_highq

# cc connected 
picc_highq_sub = read_phys_res(path_phys_res_highq, "PIcc_conn_physRes.bdio")
bcc_highq = read_phys_res(path_phys_res_highq, "PI_blc_2Qm_physRes.bdio")

picc_highq = picc_highq_sub .+ Qgev ./ (4*Qmgev) .* bcc_highq .+ 2 * charge_factor["cc"] * b33qqm_highq
#add_t0_err!(picc_highq, t0sqrt_ph_err)

# 08
pi08_highq = Vector{uwreal}(read_phys_res(path_phys_res_highq, "PI08_physRes.bdio"))

# add_t0_err!(pi08_highq, t0sqrt_ph_err)

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
fig = figure(figsize=(10,7.))
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
ylabel(L"$\widehat{\Pi}(Q^2)$")
legend()
tight_layout()
display(fig)
savefig(joinpath(path_plot, "hvp_running_highq.pdf"))
close("all")


######################################
## LOW Q
#####################################
## read phys res data LOW Q
# reading systematic errors

# syst_lowq = read_systematics(path_phys_res_lowq)

# 33
pi33_lowq_sub = read_phys_res(path_phys_res_lowq, "PI33_physRes.bdio") 
b33qqm_lowq = [(q/(4*Qmgev)) * b33qm_pt_5loops for q in Qgev ./4]
# b33qqm_lowq = [(q/(4*Qmgev^2)) * log(2) / (4pi^2) for q in Qgev ./4] # at LO
pi33_lowq = pi33_lowq_sub .+ b33qqm_lowq
# add_t0_err!(pi33_lowq, t0sqrt_ph_err)

# 88
pi33_lowq_SD = read_phys_res(path_phys_res_lowq, "PI33_SD_physRes.bdio")
# add_t0_err!(pi33_lowq_SD, t0sqrt_ph_err)
pi3388dlt_lowq = read_phys_res(path_phys_res_lowq, "PI3388_SD_physRes.bdio")
# add_t0_err!(pi3388dlt_lowq, t0sqrt_ph_err)
pi88_lowq_ILD = read_phys_res(path_phys_res_lowq, "PI88_ILD_physRes.bdio")
# add_t0_err!(pi88_lowq_ILD, t0sqrt_ph_err)

pi88_lowq = pi33_lowq_SD * charge_factor["88"] - pi3388dlt_lowq + pi88_lowq_ILD
# add_t0_err!(pi88_lowq, t0sqrt_ph_err) # use either this or the previous three lines

# cc connected
picc_lowq_sub = read_phys_res(path_phys_res_lowq, "PIcc_conn_physRes.bdio")
bcc_lowq = read_phys_res(path_phys_res_highq, "PI_blc_2Qm_physRes.bdio")
picc_lowq = picc_lowq_sub .+ 4/3 * Qgev ./ 4 ./ (4*Qmgev) .* bcc_lowq .+ 2 .* charge_factor["cc"] .* b33qqm_lowq #.+ bcc_lowq
# add_t0_err!(picc_lowq, t0sqrt_ph_err)

# 08
pi08_lowq = Vector{uwreal}(read_phys_res(path_phys_res_lowq, "PI08_physRes.bdio"))
# add_t0_err!(pi08_lowq, t0sqrt_ph_err)

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
fig = figure(figsize=(10,7.))
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
ylabel(L"$\bar\Pi(Q^2/4)$")
legend()
tight_layout()
display(fig)
savefig(joinpath(path_plot, "hvp_running_lowq.pdf"))
close("all")


######################
## TOTAL Q
#####################

# add_t0_err!(pi33_tot, t0sqrt_ph_err)
## test with t0err
# uwerr.(pi33_tot)
# pi33_tot
# add_t0_phi2_phi4_err!(pi33_tot, t0sqrt_ph_err)
# uwerr.(pi33_tot)
# pi33_tot
# details(pi33_tot[12])
##

pi33_tot = pi33_highq .+ pi33_lowq
add_t0_err!(pi33_tot, t0sqrt_ph_err)

pi88_tot = pi88_highq .+ pi88_lowq
add_t0_err!(pi88_tot, t0sqrt_ph_err)

pi08_tot = pi08_highq .+ pi08_lowq
add_t0_err!(pi08_tot, t0sqrt_ph_err)

picc_tot = picc_highq .+ picc_lowq
add_t0_err!(picc_tot, t0sqrt_ph_err)

# picc_tot = read_phys_res(path_phys_res_fullq, "PIcc_conn_physRes.bdio")
# picc_tot = picc_tot .+ [uwreal([0.0, syst_fullq["piccconn"][k]], "syst cc fullq") for k in eachindex(Qgev)]
# picc_tot ./= charge_factor["cc"] # remove charge factor 

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
fig = figure(figsize=(10,7.))
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
ylabel(L"$\bar\Pi(Q^2)$")
legend()
tight_layout()
display(fig)
savefig(joinpath(path_plot, "hvp_running_totq_scale_err.pdf"))
# savefig(joinpath(path_plot, "hvp_running_totq.pdf"))
close("all")


##
pi_tot = pi33_tot .+ pi88_tot .+ picc_tot
uwerr.(pi_tot)
mchist(pi_tot[10], "F300")

uwerr.(pi33_SD)
mchist(pi33_SD[10], "F300")
## SAVING FINAL RESULTS IN BDIO

@info("Saving final PI results")
using ALPHAio
io = IOBuffer()
write(io, "PI phys final results")
path_store_pi_fin = "/Users/alessandroconigli/MyDrive/postdoc-mainz/projects/deltalpha/final_pi_phys_results"
fb = ALPHAdobs_create(joinpath(path_store_pi_fin, "PI_all_phys_res.bdio"), io)

data = Dict{String, Array{uwreal}}(
    "pi_33_highq" => pi33_highq,
    "pi_33_lowq" => pi33_lowq,
    "pi_88_highq" => pi88_highq,
    "pi_88_lowq" => pi88_lowq,
    "pi_08_highq" => pi08_highq,
    "pi_08_lowq" => pi08_lowq,
    "pi_cc_highq" => picc_highq,
    "pi_cc_lowq" => picc_lowq
)
ALPHAdobs_write(fb, data)
ALPHAdobs_close(fb)
println("# Saving complete!")

## test reading
fb = BDIO_open(joinpath(path_store_pi_fin, "PI_all_phys_res.bdio"), "r")
res = Dict{String, Array{uwreal}}()
while ALPHAdobs_next_p(fb)
    d = ALPHAdobs_read_parameters(fb)
    nobs = d["nobs"]
    dims = d["dimensions"]
    sz = tuple(d["size"]...)
    ks = collect(d["keys"])
    # push!(res, ALPHAdobs_read_next(fb, size=sz, keys=ks))
    res =  ALPHAdobs_read_next(fb, size=sz, keys=ks)
end
BDIO_close!(fb)

## test with delta alpha
alpha_lo = 4*pi / 137.035999084
deltaalpha_running = alpha_lo .* ( pi33_tot + 4/9*picc_tot + 1/3*pi88_tot); uwerr.(deltaalpha_running)

pt_running_fix = uwreal([0.000045,0.000002], "run_Mz_minus_Mz")
pt_running = [
    uwreal([0.023945,0.000223], "q0.05"),
    uwreal([0.023945,0.000223], "q0.1"),
    uwreal([0.023945,0.000223], "q0.4"),
    uwreal([0.023945,0.000223], "q1"),
    uwreal([0.022565,0.000147], "q2"),
    uwreal([0.021729,0.000085], "q3"),
    uwreal([0.021117,0.000062], "q4"),
    uwreal([0.020628,0.000051], "q5"),
    uwreal([0.020220,0.000044], "q6"),
    uwreal([0.019868,0.000040], "q7"),
    uwreal([0.019558,0.000036], "q8"),
    uwreal([0.019280,0.000034], "q9")
]


dalpha_fin = deltaalpha_running .+ pt_running .+ pt_running_fix; uwerr.(dalpha_fin)
dalpha_fin

##
fill_between(Qgev[4:end], value.(dalpha_fin[4:end]) .-err.(dalpha_fin[4:end]), value.(dalpha_fin[4:end]) .+ err.(dalpha_fin[4:end]), alpha=0.4)
errorbar(5, 0.02773, 0.00015, fmt="s")
ylim(0.0270, 0.0282)
display(gcf())
close("all")


## comparison of charm plot (OLD)
picc_tot = picc_highq .+ picc_lowq
change_id(from="D450", to="D450_3")
change_id(from="E300", to="E300_3")

picc_tot_fullq = read_phys_res("/Users/alessandroconigli/MyDrive/postdoc-mainz/projects/deltalpha/physical_results/scale_error_multimom/full_q", "PIcc_conn_physRes.bdio")
picc_tot_fullq = picc_tot_fullq .+ [uwreal([0.0, 0.0], "syst cc fullq") for k in eachindex(Qgev)]
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
# savefig(joinpath(path_plot, "charm_comparison_totq.pdf"))
close("all")
