using HVPobs
using ALPHAio
using BDIO
using ADerrors
import ADerrors: err
using PyPlot, LaTeXStrings

path_store_pi_fin = "/Users/alessandroconigli/MyDrive/postdoc-mainz/projects/deltalpha/final_pi_phys_results/"
path_plot = "/Users/alessandroconigli/MyDrive/postdoc-mainz/projects/deltalpha/plots/error_comparison"
#================ SET UP PLOT PARAMETERS ==================#

rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
rcParams["text.usetex"] =  true
rcParams["mathtext.fontset"]  = "cm"
rcParams["font.size"] =13
rcParams["axes.labelsize"] =22
rcParams["axes.titlesize"] = 18
plt.rc("text", usetex=true) # set to true if a LaTeX installation is present

const Qgev = [0.05, 0.1, 0.4, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0] # Q^2


include("../utils/IO_BDIO.jl")
include("./tools_running.jl")
include("../utils/const.jl")

fb = BDIO_open(joinpath(path_store_pi_fin, "PI_all_phys_res.bdio"), "r")
res = Dict{String, Array{uwreal}}()
while ALPHAdobs_next_p(fb)
    d = ALPHAdobs_read_parameters(fb)
    nobs = d["nobs"]
    dims = d["dimensions"]
    sz = tuple(d["size"]...)
    ks = collect(d["keys"])
    res =  ALPHAdobs_read_next(fb, size=sz, keys=ks)
end
BDIO_close!(fb)

pi33_lowq = res["pi_33_lowq"]
pi33_highq = res["pi_33_highq"]
add_t0_err!(pi33_lowq, t0sqrt_ph_err)
add_t0_err!(pi33_highq, t0sqrt_ph_err)
pi33 = pi33_highq .+ pi33_lowq

pi88_lowq = res["pi_88_lowq"]
pi88_highq = res["pi_88_highq"]
add_t0_err!(pi88_lowq, t0sqrt_ph_err)
add_t0_err!(pi88_highq, t0sqrt_ph_err)
pi88 = pi88_highq .+ pi88_lowq

picc_lowq = res["pi_cc_lowq"]
picc_highq = res["pi_cc_highq"]
add_t0_err!(picc_highq, t0sqrt_ph_err)
add_t0_err!(picc_lowq, t0sqrt_ph_err)
picc = picc_highq .+ picc_lowq

pi08_lowq = res["pi_08_lowq"]
pi08_highq = res["pi_08_highq"]
add_t0_err!(pi08_highq, t0sqrt_ph_err)
add_t0_err!(pi08_lowq, t0sqrt_ph_err)
pi08 = pi08_highq .+ pi08_lowq

pi_tot_lowq = pi33_lowq + charge_factor["88"] * pi88_lowq + charge_factor["cc"]*picc_lowq
pi_tot_highq = pi33_highq + charge_factor["88"] * pi88_highq + charge_factor["cc"]*picc_highq
pi_tot = pi_tot_lowq + pi_tot_highq


## variance over Q for the total piece

var_highq_perc = []
var_lowq_perc = []
for k in eachindex(pi_tot)
    uwerr(pi_tot[k])
    uwerr(pi_tot_highq[k])
    uwerr(pi_tot_lowq[k])

    var_lowq = err(pi_tot_lowq[k])^2
    var_highq = err(pi_tot_highq[k])^2
    var_tot = err(pi_tot[k])^2
    push!(var_highq_perc, var_highq /  (var_highq+var_lowq))
    push!(var_lowq_perc, var_lowq / (var_highq+var_lowq))
end
##
fill_between(Qgev, zeros(length(Qgev)), var_lowq_perc, label=L"$\sigma^2(\bar\Pi(Q^2/4))$")
fill_between(Qgev, var_lowq_perc, var_highq_perc+var_lowq_perc, label=L"$\sigma^2(\widehat\Pi(Q^2))$")

tight_layout()
legend()
xlabel(L"$Q^2 \ \mathrm{GeV}^2$")
ylabel(L"$\sigma^2(\bar\Pi^{(\gamma,\gamma)})$")
tight_layout()
display(gcf())
savefig(joinpath(path_plot, "variance_pigammagamma_tot.pdf"))
close("all")

## plotting the precision in % for lowq highq tot
uwerr.(pi33_lowq); uwerr.(pi33_highq)
uwerr.(pi88_lowq); uwerr.(pi88_highq)
uwerr.(picc_lowq); uwerr.(picc_highq)
uwerr.(pi33)
uwerr.(pi88)
uwerr.(picc)
uwerr.(pi08)
fig = figure(figsize=(10, 8))
# subplots_adjust(wspace=0.0)

subplot(311) # lowq
rel_prec_tot = err.(pi_tot_lowq) ./ value.(pi_tot_lowq) .*100
rel_prec_33 = err.(pi33_lowq) ./ value.(pi33_lowq) .*100
rel_prec_88 = err.(pi88_lowq) ./ value.(pi88_lowq) .*100
rel_prec_cc = err.(picc_lowq) ./ value.(picc_lowq) .*100

plot(Qgev./4, rel_prec_tot, marker="o", lw=3, ms=4, label="tot")
plot(Qgev./4, rel_prec_33, ms=4, marker="o", lw=2, ls="dashed", label="I=1")
plot(Qgev./4, rel_prec_88, ms=4, marker="o", lw=2, ls="dashed", label="I=0")
plot(Qgev./4, rel_prec_cc, ms=4, marker="o", lw=2, ls="dashed", label="charm")

ylabel(L"$\delta \bar{\Pi}({Q^2/4}) \ [\%] $")

legend()

subplot(312) # highq
rel_prec_tot = err.(pi_tot_highq) ./ value.(pi_tot_highq) .*100
rel_prec_33 = err.(pi33_highq) ./ value.(pi33_highq) .*100
rel_prec_88 = err.(pi88_highq) ./ value.(pi88_highq) .*100
rel_prec_cc = err.(picc_highq) ./ value.(picc_highq) .*100

plot(Qgev, rel_prec_tot, label="tot", lw=3, ms=4, marker="o")
plot(Qgev, rel_prec_33, ms=4, marker="o", lw=2, ls="dashed", label="I=1")
plot(Qgev, rel_prec_88, ms=4, marker="o", lw=2, ls="dashed", label="I=0")
plot(Qgev, rel_prec_cc, ms=4, marker="o", lw=2, ls="dashed", label="charm")
xlim(0.5,9.5)
ylim(0,2)
ylabel(L"$\delta \widehat{\Pi}({Q^2}) \ [\%] $")

legend(ncol=2)

subplot(313) # tot
rel_prec_tot = err.(pi_tot) ./ value.(pi_tot) .*100
rel_prec_33 = err.(pi33) ./ value.(pi33) .*100
rel_prec_88 = err.(pi88) ./ value.(pi88) .*100
rel_prec_cc = err.(picc) ./ value.(picc) .*100

plot(Qgev, rel_prec_tot, label="tot", lw=3, ms=4, marker="o")
plot(Qgev, rel_prec_33,  ms=4, marker="o", lw=2, ls="dashed", label="I=1")
plot(Qgev, rel_prec_88,  ms=4, marker="o", lw=2, ls="dashed", label="I=0")
plot(Qgev, rel_prec_cc,  ms=4, marker="o", lw=2, ls="dashed", label="charm")
xlim(0.5,9.5)
ylim(0,2)
legend(ncol=2)
ylabel(L"$\delta \bar{\Pi}({Q^2}) \ [\%] $")
xlabel(L"Q^2 \ [\mathrm{GeV}]^2")


tight_layout()
display(fig)
savefig(joinpath(path_plot, "err_vs_Q_per_component.pdf"))
close("all")


##
for kk in collect(keys(res))
    uwerr.(res[kk])
end
##

uwerr.(pi33)
uwerr.(pi88)
uwerr.(picc)
uwerr.(pi08)
fig = figure(figsize=(10, 8))
subplots_adjust(wspace=0.0)

subplot(131) # 3 GeV^2
title(L"$Q^2 = 3 \ \mathrm{GeV}^2$")
# pi33
rel_err = ( err(pi33[6]) ) / value(pi33[6]) * 100
rel_err_no_scale = 0.00024227952970388953 / value(pi33[6]) * 100
rel_err_paper = 0.0006 / 0.0488 * 100
rel_err_paper_no_scale=0.0005 / 0.0488 * 100
axhline(5, 0.0, rel_err / 2.5, color="#0072B2", lw=3,alpha=0.5, label=L"$\mathrm{This \ work}$")
axhline(5, 0.0, rel_err_no_scale / 2.5, color="#0072B2", lw=5)

axhline(4.9, 0.0, rel_err_paper / 2.5, color="#009E73", lw=3, alpha=0.5, label=L"$\mathrm{2022\ paper}$")
axhline(4.9, 0.0, rel_err_paper_no_scale / 2.5, color="#009E73", lw=5)

# pi88
rel_err = ( err(pi88[6]) ) / value(pi88[6]) * 100
rel_err_no_scale = 0.00022784386239015593 / value(pi88[6]) * 100
rel_err_paper = 0.0006 / 0.0374 * 100
rel_err_paper_noscale = 0.0005 / 0.0374 * 100

axhline(4, 0.0, rel_err / 2.5, color="#0072B2", alpha=0.5, lw=3)
axhline(4, 0.0, rel_err_no_scale / 2.5, color="#0072B2", lw=5)

axhline(3.9, 0.0, rel_err_paper / 2.5, alpha=0.5, color="#009E73", lw=3)
axhline(3.9, 0.0, rel_err_paper_no_scale / 2.5, color="#009E73", lw=5)

# pi 08
rel_err = ( err(pi08[6]) ) / value(pi08[6]) * 100
rel_err_no_scale = 7.868738982574486e-5 / value(pi08[6]) * 100
rel_err_paper = 0.00017 / 0.00690 * 100
rel_err_paper_no_scale = 0.00008 / 0.00698 * 100
axhline(3, 0.0, rel_err / 2.5, color="#0072B2", lw=3, alpha=0.5)
axhline(3, 0.0, rel_err_no_scale / 2.5, color="#0072B2", lw=5, zorder=1)

axhline(2.9, 0.0, rel_err_paper / 2.5, color="#009E73", lw=3, alpha=0.5)
axhline(2.9, 0.0, rel_err_paper_no_scale / 2.5, color="#009E73", lw=5, zorder=1)

# pi cc
# println("YOU ARE SUBTRACTING 1 FROM PICC ERR AT 3 GEV2")
rel_err = ( err(picc[6]) ) / value(picc[6]) * 100 
rel_err_no_scale = 0.0001003605583073771 / value(picc[6]) * 100 

rel_err_paper = 0.00019 / 0.01064 * 100
rel_err_paper_no_scale = 0.00005 / 0.01064 * 100

axhline(2, 0.0, rel_err / 2.5, color="#0072B2", lw=3, alpha=0.5)
axhline(2, 0.0, rel_err_no_scale / 2.5, color="#0072B2", lw=5, zorder=1)

axhline(1.9, 0.0, rel_err_paper / 2.5, color="#009E73", lw=3, alpha=0.5)
axhline(1.9, 0.0, rel_err_paper_no_scale / 2.5, color="#009E73", lw=5, zorder=1)

xlim(0, 2.5)
ylim(1.5,5.5)
yticks([2,3,4,5], [L"$\bar\mathit{\Pi}^{cc}$", L"$\bar\mathit{\Pi}^{08}$", L"$\bar\mathit{\Pi}^{88}$", L"$\bar\mathit{\Pi}^{33}$"], size=22)
xlabel(L"$\mathrm{err} \ [\%]$")
legend()


subplot(132) # 5 GeV^2
ax1 = gca()
title(L"$Q^2 = 5 \ \mathrm{GeV}^2$")
setp(ax1.get_yticklines(),visible=false) # Disable x tick lines
setp(ax1.get_yticklabels(),visible=false) # Disable x tick labels

# pi33
rel_err = ( err(pi33[8]) ) / value(pi33[8]) * 100
rel_err_no_scale = 0.00031762426868129327 / value(pi33[8]) * 100

rel_err_paper = 0.0008 / 0.0560 * 100
rel_err_paper_no_scale = 0.00061 / 0.0560 * 100

axhline(5, 0.0, rel_err / 2.5, color="#0072B2", lw=3, alpha=0.5)
axhline(5, 0.0, rel_err_no_scale / 2.5, color="#0072B2", lw=5)
axhline(4.9, 0.0, rel_err_paper / 2.5, color="#009E73", lw=3, alpha=0.5)
axhline(4.9, 0.0, rel_err_paper_no_scale / 2.5, color="#009E73", lw=3)

# pi88
rel_err = ( err(pi88[8]) ) / value(pi88[8]) * 100
rel_err_no_scale = 0.0002517374410905558 / value(pi88[8]) * 100

rel_err_paper = 0.0008 / 0.0445 * 100
rel_err_paper_no_scale = 0.00062 / 0.0445 * 100
axhline(4, 0.0, rel_err / 2.5, color="#0072B2", lw=3, alpha=0.5)
axhline(4, 0.0, rel_err_no_scale / 2.5, color="#0072B2", lw=5)

axhline(3.9, 0.0, rel_err_paper / 2.5, color="#009E73", lw=3, alpha=0.5)
axhline(3.9, 0.0, rel_err_paper_no_scale / 2.5, color="#009E73", lw=5)

# pi 08
rel_err = ( err(pi08[8]) ) / value(pi08[8]) * 100
rel_err_no_scale = 8.475841226332736e-5 / value(pi08[8]) * 100

rel_err_paper = 0.00017 / 0.00701 * 100
rel_err_paper_no_scale = 0.00009 / 0.00701 * 100

axhline(3, 0.0, rel_err / 2.5, color="#0072B2", lw=3, alpha=0.5)
axhline(3, 0.0, rel_err_no_scale / 2.5, color="#0072B2", lw=5, zorder=1)

axhline(2.9, 0.0, rel_err_paper / 2.5, color="#009E73", lw=3, alpha=0.5)
axhline(2.9, 0.0, rel_err_paper_no_scale / 2.5, color="#009E73", lw=5, zorder=1)

# pi cc
# println("YOU ARE SUBTRACTING 1 FROM PICC ERR AT 5 GEV2")
rel_err = ( err(picc[8]) ) / value(picc[8]) * 100 
rel_err_no_scale =  0.00015662126616213247 / value(picc[8]) * 100 

rel_err_paper = 0.00027 / 0.01608 * 100
rel_err_paper_no_scale = 0.00007 / 0.01608 * 100

axhline(2, 0.0, rel_err / 2.5, color="#0072B2", lw=3, alpha=0.5)
axhline(2, 0.0, rel_err_no_scale / 2.5, color="#0072B2", lw=5)

axhline(1.9, 0.0, rel_err_paper / 2.5, color="#009E73", lw=3, alpha=0.5)
axhline(1.9, 0.0, rel_err_paper_no_scale / 2.5, color="#009E73", lw=5)

xlim(0, 2.5)
ylim(1.5,5.5)
xlabel(L"$\mathrm{err} \ [\%]$")


subplot(133)
title(L"$Q^2 = 7 \ \mathrm{GeV}^2$")
ax1 = gca()
setp(ax1.get_yticklines(),visible=false) # Disable x tick lines
setp(ax1.get_yticklabels(),visible=false) # Disable x tick labels

# pi33
rel_err = ( err(pi33[10]) ) / value(pi33[10]) * 100
rel_err_no_scale = 0.00034206154465770176 / value(pi33[10]) * 100

rel_err_paper = 0.0008 / 0.0608 * 100
rel_err_paper_no_scale = 0.00063 / 0.0608 * 100
axhline(5, 0.0, rel_err / 2.5, color="#0072B2", lw=3, alpha=0.5)
axhline(5, 0.0, rel_err_no_scale / 2.5, color="#0072B2", lw=5)

axhline(4.9, 0.0, rel_err_paper / 2.5, color="#009E73", lw=3, alpha=0.5)
axhline(4.9, 0.0, rel_err_paper_no_scale / 2.5, color="#009E73", lw=5)

# pi88
rel_err = ( err(pi88[10]) ) / value(pi88[10]) * 100
rel_err_no_scale = 0.00026019080979986913 / value(pi88[10]) * 100

rel_err_paper = 0.0008 / 0.0493 * 100
rel_err_paper_no_scale = 0.00064 / 0.0493 * 100

axhline(4, 0.0, rel_err / 2.5, color="#0072B2", lw=3, alpha=0.5)
axhline(4, 0.0, rel_err_no_scale / 2.5, color="#0072B2", lw=5)

axhline(3.9, 0.0, rel_err_paper / 2.5, color="#009E73", lw=3, alpha=0.5)
axhline(3.9, 0.0, rel_err_paper_no_scale / 2.5, color="#009E73", lw=5)

# pi 08
rel_err = ( err(pi08[10]) ) / value(pi08[10]) * 100
rel_err_no_scale = 9.312547526919385e-5 / value(pi08[10]) * 100

rel_err_paper = 0.00017 / 0.00704 * 100
rel_err_paper_no_scale = 0.00008 / 0.00704 * 100

axhline(3, 0.0, rel_err / 2.5, color="#0072B2", lw=3, alpha=0.5)
axhline(3, 0.0, rel_err_no_scale / 2.5, color="#0072B2", lw=5)

axhline(2.9, 0.0, rel_err_paper / 2.5, color="#009E73", lw=3, alpha=0.5)
axhline(2.9, 0.0, rel_err_paper_no_scale / 2.5, color="#009E73", lw=5)

# pi cc
# println("YOU ARE SUBTRACTING 1 FROM PICC ERR AT 7 GEV2")
rel_err = ( err(picc[10]) ) / value(picc[10]) * 100 
rel_err_no_scale = 0.00019218180866509472 / value(picc[10]) * 100 

rel_err_paper = 0.00032 / 0.02066 * 100
rel_err_paper_no_scale = 0.00008 / 0.02066 * 100

axhline(2, 0.0, rel_err / 2.5, color="#0072B2", lw=3, alpha=0.5)
axhline(2, 0.0, rel_err_no_scale / 2.5, color="#0072B2", lw=5)

axhline(1.9, 0.0, rel_err_paper / 2.5, color="#009E73", lw=3, alpha=0.5)
axhline(1.9, 0.0, rel_err_paper_no_scale / 2.5, color="#009E73", lw=5)

xlim(0, 2.5)
ylim(1.5,5.5)
xlabel(L"$\mathrm{err} \ [\%]$")



tight_layout()
display(fig)
savefig(joinpath(path_plot, "error_comparison.pdf"))
close("all")



## PIE CHART PI33

fig = figure(figsize=(10, 4))
subplots_adjust(wspace=0.0)

subplot(131) # 3 GeV^2
title(L"$ \mathit{\bar{\Pi}}^{(3,3)}(Q^2 = 3 \ \mathrm{GeV}^2)$", fontsize=12)

labels = ("Systematics",  string("Stat. +",L"$\chi$", "-cont. limit"), "Scale setting", "Others")
sizes = [60.11, 22.63, 15.7]
push!(sizes, 100-sum(sizes))
explode = fill(0.03,4)
colors = ["#E69F00", "#56B4E9", "#009E73", "#CC79A7"]
pie(sizes, explode=explode, colors=colors)
legend(labels, fontsize=8, bbox_to_anchor=(0.5,0.9), loc="upper right")

subplot(132) # 5 GeV^2
title(L"$ \mathit{\bar{\Pi}}^{(3,3)}(Q^2 = 5 \ \mathrm{GeV}^2)$", fontsize=12)
labels = ("Systematics",  string("Stat. +",L"$\chi$", "-cont. limit"), "Scale setting", "Others")
sizes = [64.7, 24.5, 8.92]
push!(sizes, 100-sum(sizes))
explode = fill(0.03,4)
colors = ["#E69F00", "#56B4E9", "#009E73", "#CC79A7"]
pie(sizes, explode=explode, colors=colors)
legend(labels, fontsize=8, bbox_to_anchor=(0.5,0.9), loc="upper right")

subplot(133) # 9 GeV^2
title(L"$ \mathit{\bar{\Pi}}^{(3,3)}(Q^2 = 9 \ \mathrm{GeV}^2)$", fontsize=12)
labels = ("Systematics",  string("Stat. +",L"$\chi$", "-cont. limit"), "Scale setting", "Others")
sizes = [64, 26, 7.8]
push!(sizes, 100-sum(sizes))
explode = fill(0.03,4)
colors = ["#E69F00", "#56B4E9", "#009E73", "#CC79A7"]
pie(sizes, explode=explode, colors=colors)
legend(labels, fontsize=8, bbox_to_anchor=(0.5,0.9), loc="upper right")

tight_layout()
display(fig)
savefig(joinpath(path_plot, "piechart_pi33.pdf"))
close("all")


## PIE CHART PI88

fig = figure(figsize=(10, 4))
subplots_adjust(wspace=0.0)

subplot(131) # 3 GeV^2
title(L"$ \mathit{\bar{\Pi}}^{(8,8)}(Q^2 = 3 \ \mathrm{GeV}^2)$", fontsize=12)

labels = ("Systematics",  string("Stat. +",L"$\chi$", "-cont. limit"), "Scale setting", "Others")
sizes = [37.5, 27.2 , 34.12]
push!(sizes, 100-sum(sizes))
explode = fill(0.03,4)
colors = ["#E69F00", "#56B4E9", "#009E73", "#CC79A7"]
pie(sizes, explode=explode, colors=colors)
legend(labels, fontsize=8, bbox_to_anchor=(0.5,0.9), loc="upper right")

subplot(132) # 5 GeV^2
title(L"$ \mathit{\bar{\Pi}}^{(8,8)}(Q^2 = 5 \ \mathrm{GeV}^2)$", fontsize=12)
labels = ("Systematics",  string("Stat. +",L"$\chi$", "-cont. limit"), "Scale setting", "Others")
sizes = [36.4, 26.0, 36.9]
push!(sizes, 100-sum(sizes))
explode = fill(0.03,4)
colors = ["#E69F00", "#56B4E9", "#009E73", "#CC79A7"]
pie(sizes, explode=explode, colors=colors)
legend(labels, fontsize=8, bbox_to_anchor=(0.5,0.9), loc="upper right")

subplot(133) # 9 GeV^2
title(L"$ \mathit{\bar{\Pi}}^{(8,8)}(Q^2 = 9 \ \mathrm{GeV}^2)$", fontsize=12)
labels = ("Systematics",  string("Stat. +",L"$\chi$", "-cont. limit"), "Scale setting", "Others")
sizes = [31.5, 25.6, 41.8]
push!(sizes, 100-sum(sizes))
explode = fill(0.03,4)
colors = ["#E69F00", "#56B4E9", "#009E73", "#CC79A7"]
pie(sizes, explode=explode, colors=colors)
legend(labels, fontsize=8, bbox_to_anchor=(0.5,0.9), loc="upper right")

tight_layout()
display(fig)
savefig(joinpath(path_plot, "piechart_pi88.pdf"))
close("all")


## PIE CHART PI08

fig = figure(figsize=(10, 4))
subplots_adjust(wspace=0.0)

subplot(131) # 3 GeV^2
title(L"$ \mathit{\bar{\Pi}}^{(0,8)}(Q^2 = 3 \ \mathrm{GeV}^2)$", fontsize=12)

labels = ("Systematics",  string("Stat. +",L"$\chi$", "-cont. limit"), "Scale setting", "Others")
sizes = [3, 75.2 , 20.69]
push!(sizes, 100-sum(sizes))
explode = fill(0.03,4)
colors = ["#E69F00", "#56B4E9", "#009E73", "#CC79A7"]
pie(sizes, explode=explode, colors=colors)
legend(labels, fontsize=8, bbox_to_anchor=(0.5,0.9), loc="upper right")

subplot(132) # 5 GeV^2
title(L"$ \mathit{\bar{\Pi}}^{(0,8)}(Q^2 = 5 \ \mathrm{GeV}^2)$", fontsize=12)
labels = ("Systematics",  string("Stat. +",L"$\chi$", "-cont. limit"), "Scale setting", "Others")
sizes = [2.8, 75.2, 20.77]
push!(sizes, 100-sum(sizes))
explode = fill(0.03,4)
colors = ["#E69F00", "#56B4E9", "#009E73", "#CC79A7"]
pie(sizes, explode=explode, colors=colors)
legend(labels, fontsize=8, bbox_to_anchor=(0.5,0.9), loc="upper right")

subplot(133) # 9 GeV^2
title(L"$ \mathit{\bar{\Pi}}^{(0,8)}(Q^2 = 9 \ \mathrm{GeV}^2)$", fontsize=12)
labels = ("Systematics",  string("Stat. +",L"$\chi$", "-cont. limit"), "Scale setting", "Others")
sizes = [5, 75.1, 18.7]
push!(sizes, 100-sum(sizes))
explode = fill(0.03,4)
colors = ["#E69F00", "#56B4E9", "#009E73", "#CC79A7"]
pie(sizes, explode=explode, colors=colors)
legend(labels, fontsize=8, bbox_to_anchor=(0.5,0.9), loc="upper right")

tight_layout()
display(fig)
savefig(joinpath(path_plot, "piechart_pi08.pdf"))
close("all")


## PIE CHART PICC

fig = figure(figsize=(10, 4))
subplots_adjust(wspace=0.0)

subplot(131) # 3 GeV^2
title(L"$ \mathit{\bar{\Pi}}^{(c,c)}(Q^2 = 3 \ \mathrm{GeV}^2)$", fontsize=12)

labels = ("Systematics",  string("Stat. +",L"$\chi$", "-cont. limit"), "Scale setting", "Others")
sizes = [16, 57.6 , 25.56]
push!(sizes, 100-sum(sizes))
explode = fill(0.03,4)
colors = ["#E69F00", "#56B4E9", "#009E73", "#CC79A7"]
pie(sizes, explode=explode, colors=colors)
legend(labels, fontsize=8, bbox_to_anchor=(0.5,0.9), loc="upper right")

subplot(132) # 5 GeV^2
title(L"$ \mathit{\bar{\Pi}}^{(c,c)}(Q^2 = 5 \ \mathrm{GeV}^2)$", fontsize=12)
labels = ("Systematics",  string("Stat. +",L"$\chi$", "-cont. limit"), "Scale setting", "Others")
sizes = [11, 63., 25.2]
push!(sizes, 100-sum(sizes))
explode = fill(0.03,4)
colors = ["#E69F00", "#56B4E9", "#009E73", "#CC79A7"]
pie(sizes, explode=explode, colors=colors)
legend(labels, fontsize=8, bbox_to_anchor=(0.5,0.9), loc="upper right")

subplot(133) # 9 GeV^2
title(L"$ \mathit{\bar{\Pi}}^{(c,c)}(Q^2 = 9 \ \mathrm{GeV}^2)$", fontsize=12)
labels = ("Systematics",  string("Stat. +",L"$\chi$", "-cont. limit"), "Scale setting", "Others")
sizes = [6, 60.7, 32.3]
push!(sizes, 100-sum(sizes))
explode = fill(0.03,4)
colors = ["#E69F00", "#56B4E9", "#009E73", "#CC79A7"]
pie(sizes, explode=explode, colors=colors)
legend(labels, fontsize=8, bbox_to_anchor=(0.5,0.9), loc="upper right")

tight_layout()
display(fig)
savefig(joinpath(path_plot, "piechart_picc.pdf"))
close("all")





## PIE CHART PIgamma gamma at low energies

fig = figure(figsize=(10, 4))
subplots_adjust(wspace=0.0)

subplot(131) # 0.0125 GeV^2 (100 MeV)
title(L"$ \mathit{\bar{\Pi}}^{(\gamma,\gamma)}(Q^2 = 0.0125 \ \mathrm{GeV}^2) \ [1.9 \ \%]$", fontsize=12)

labels = ("Systematics",  string("Stat. +",L"$\chi$", "-cont. limit"), "Scale setting", "Others")
sizes = [60.4, 26.63, 11.7]
push!(sizes, 100-sum(sizes))
explode = fill(0.03,4)
colors = ["#E69F00", "#56B4E9", "#009E73", "#CC79A7"]
pie(sizes, explode=explode, colors=colors)
legend(labels, fontsize=8, bbox_to_anchor=(0.5,0.9), loc="upper right")

subplot(132) # 0.1 GeV^2
title(L"$ \mathit{\bar{\Pi}}^{(\gamma,\gamma)}(Q^2 = 0.1 \ \mathrm{GeV}^2) \ [1.3\ \%]$", fontsize=12)
labels = ("Systematics",  string("Stat. +",L"$\chi$", "-cont. limit"), "Scale setting", "Others")
sizes = [52, 28.5, 18]
push!(sizes, 100-sum(sizes))
explode = fill(0.03,4)
colors = ["#E69F00", "#56B4E9", "#009E73", "#CC79A7"]
pie(sizes, explode=explode, colors=colors)
legend(labels, fontsize=8, bbox_to_anchor=(0.5,0.9), loc="upper right")

subplot(133) # 0.5 GeV^2
title(L"$ \mathit{\bar{\Pi}}^{(\gamma,\gamma)}(Q^2 = 0.5 \ \mathrm{GeV}^2) \ [0.9\ \%]$", fontsize=12)
labels = ("Systematics",  string("Stat. +",L"$\chi$", "-cont. limit"), "Scale setting", "Others")
sizes = [48, 31, 19]
push!(sizes, 100-sum(sizes))
explode = fill(0.03,4)
colors = ["#E69F00", "#56B4E9", "#009E73", "#CC79A7"]
pie(sizes, explode=explode, colors=colors)
legend(labels, fontsize=8, bbox_to_anchor=(0.5,0.9), loc="upper right")

tight_layout()
display(fig)
savefig(joinpath(path_plot, "piechart_pigammagamma_lowq.pdf"))
close("all")