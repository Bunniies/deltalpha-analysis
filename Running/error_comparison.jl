using HVPobs
using ALPHAio
using BDIO
using ADerrors
import ADerrors: err
using PyPlot, LaTeXStrings

path_store_pi_fin = "/Users/alessandroconigli/MyDrive/postdoc-mainz/projects/deltalpha/final_pi_phys_results"
path_plot = "/Users/alessandroconigli/MyDrive/postdoc-mainz/projects/deltalpha/plots/error_comparison"
#================ SET UP PLOT PARAMETERS ==================#

rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
rcParams["text.usetex"] =  true
rcParams["mathtext.fontset"]  = "cm"
rcParams["font.size"] =13
rcParams["axes.labelsize"] =22
rcParams["axes.titlesize"] = 18
plt.rc("text", usetex=true) # set to true if a LaTeX installation is present



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


##
for kk in collect(keys(res))
    uwerr.(res[kk])
end
##


fig = figure(figsize=(10, 8))
subplots_adjust(wspace=0.0)

subplot(131) # 3 GeV^2
title(L"$Q^2 = 3 \ \mathrm{GeV}^2$")
# pi33
rel_err = ( err(res["pi_33"][6]) ) / value(res["pi_33"][6]) * 100
rel_err_no_scale = 0.00024227952970388953 / value(res["pi_33"][6]) * 100
rel_err_paper = 0.0006 / 0.0488 * 100
rel_err_paper_no_scale=0.0005 / 0.0488 * 100
axhline(5, 0.0, rel_err / 2.5, color="#0072B2", lw=3,alpha=0.5, label=L"$\mathrm{This \ work}$")
axhline(5, 0.0, rel_err_no_scale / 2.5, color="#0072B2", lw=5, label=L"$\mathrm{This \ work}$")

axhline(4.9, 0.0, rel_err_paper / 2.5, color="#009E73", lw=3, alpha=0.5, label=L"$\mathrm{2022\ paper}$")
axhline(4.9, 0.0, rel_err_paper_no_scale / 2.5, color="#009E73", lw=5, label=L"$\mathrm{2022\ paper}$")

# pi88
rel_err = ( err(res["pi_88"][6]) ) / value(res["pi_88"][6]) * 100
rel_err_no_scale = 0.00022784386239015593 / value(res["pi_88"][6]) * 100
rel_err_paper = 0.0006 / 0.0374 * 100
rel_err_paper_noscale = 0.0005 / 0.0374 * 100

axhline(4, 0.0, rel_err / 2.5, color="#0072B2", alpha=0.5, lw=3)
axhline(4, 0.0, rel_err_no_scale / 2.5, color="#0072B2", lw=5)

axhline(3.9, 0.0, rel_err_paper / 2.5, alpha=0.5, color="#009E73", lw=3)
axhline(3.9, 0.0, rel_err_paper_no_scale / 2.5, color="#009E73", lw=5)

# pi 08
rel_err = ( err(res["pi_08"][6]) ) / value(res["pi_08"][6]) * 100
rel_err_no_scale = 7.868738982574486e-5 / value(res["pi_08"][6]) * 100
rel_err_paper = 0.00017 / 0.00690 * 100
rel_err_paper_no_scale = 0.00008 / 0.00698 * 100
axhline(3, 0.0, rel_err / 2.5, color="#0072B2", lw=3, alpha=0.5)
axhline(3, 0.0, rel_err_no_scale / 2.5, color="#0072B2", lw=5, zorder=1)

axhline(2.9, 0.0, rel_err_paper / 2.5, color="#009E73", lw=3, alpha=0.5)
axhline(2.9, 0.0, rel_err_paper_no_scale / 2.5, color="#009E73", lw=5, zorder=1)

# pi cc
# println("YOU ARE SUBTRACTING 1 FROM PICC ERR AT 3 GEV2")
rel_err = ( err(res["pi_cc"][6]) ) / value(res["pi_cc"][6]) * 100 
rel_err_no_scale = 0.0001003605583073771 / value(res["pi_cc"][6]) * 100 

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
rel_err = ( err(res["pi_33"][8]) ) / value(res["pi_33"][8]) * 100
rel_err_no_scale = 0.00026082519574522223 / value(res["pi_33"][8]) * 100

rel_err_paper = 0.0008 / 0.0560 * 100
rel_err_paper_no_scale = 0.00061 / 0.0560 * 100

axhline(5, 0.0, rel_err / 2.5, color="#0072B2", lw=3, alpha=0.5)
axhline(5, 0.0, rel_err_no_scale / 2.5, color="#0072B2", lw=5)
axhline(4.9, 0.0, rel_err_paper / 2.5, color="#009E73", lw=3, alpha=0.5)
axhline(4.9, 0.0, rel_err_paper_no_scale / 2.5, color="#009E73", lw=3)

# pi88
rel_err = ( err(res["pi_88"][8]) ) / value(res["pi_88"][8]) * 100
rel_err_no_scale = 0.0002517374410905558 / value(res["pi_88"][8]) * 100

rel_err_paper = 0.0008 / 0.0445 * 100
rel_err_paper_no_scale = 0.00062 / 0.0445 * 100
axhline(4, 0.0, rel_err / 2.5, color="#0072B2", lw=3, alpha=0.5)
axhline(4, 0.0, rel_err_no_scale / 2.5, color="#0072B2", lw=5)

axhline(3.9, 0.0, rel_err_paper / 2.5, color="#009E73", lw=3, alpha=0.5)
axhline(3.9, 0.0, rel_err_paper_no_scale / 2.5, color="#009E73", lw=5)

# pi 08
rel_err = ( err(res["pi_08"][8]) ) / value(res["pi_08"][8]) * 100
rel_err_no_scale = 7.851889552228499e-5 / value(res["pi_08"][8]) * 100

rel_err_paper = 0.00017 / 0.00701 * 100
rel_err_paper_no_scale = 0.00009 / 0.00701 * 100

axhline(3, 0.0, rel_err / 2.5, color="#0072B2", lw=3, alpha=0.5)
axhline(3, 0.0, rel_err_no_scale / 2.5, color="#0072B2", lw=5, zorder=1)

axhline(2.9, 0.0, rel_err_paper / 2.5, color="#009E73", lw=3, alpha=0.5)
axhline(2.9, 0.0, rel_err_paper_no_scale / 2.5, color="#009E73", lw=5, zorder=1)

# pi cc
# println("YOU ARE SUBTRACTING 1 FROM PICC ERR AT 5 GEV2")
rel_err = ( err(res["pi_cc"][8]) ) / value(res["pi_cc"][8]) * 100 
rel_err_no_scale =  0.00013200748487943054 / value(res["pi_cc"][8]) * 100 

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
rel_err = ( err(res["pi_33"][10]) ) / value(res["pi_33"][10]) * 100
rel_err_no_scale = 0.00025591981780368785 / value(res["pi_33"][10]) * 100

rel_err_paper = 0.0008 / 0.0608 * 100
rel_err_paper_no_scale = 0.00063 / 0.0608 * 100
axhline(5, 0.0, rel_err / 2.5, color="#0072B2", lw=3, alpha=0.5)
axhline(5, 0.0, rel_err_no_scale / 2.5, color="#0072B2", lw=5)

axhline(4.9, 0.0, rel_err_paper / 2.5, color="#009E73", lw=3, alpha=0.5)
axhline(4.9, 0.0, rel_err_paper_no_scale / 2.5, color="#009E73", lw=5)

# pi88
rel_err = ( err(res["pi_88"][10]) ) / value(res["pi_88"][10]) * 100
rel_err_no_scale = 0.0002752680320349962 / value(res["pi_88"][10]) * 100

rel_err_paper = 0.0008 / 0.0493 * 100
rel_err_paper_no_scale = 0.00064 / 0.0493 * 100

axhline(4, 0.0, rel_err / 2.5, color="#0072B2", lw=3, alpha=0.5)
axhline(4, 0.0, rel_err_no_scale / 2.5, color="#0072B2", lw=5)

axhline(3.9, 0.0, rel_err_paper / 2.5, color="#009E73", lw=3, alpha=0.5)
axhline(3.9, 0.0, rel_err_paper_no_scale / 2.5, color="#009E73", lw=5)

# pi 08
rel_err = ( err(res["pi_08"][10]) ) / value(res["pi_08"][10]) * 100
rel_err_no_scale = 7.819447774689946e-5 / value(res["pi_08"][10]) * 100

rel_err_paper = 0.00017 / 0.00704 * 100
rel_err_paper_no_scale = 0.00008 / 0.00704 * 100

axhline(3, 0.0, rel_err / 2.5, color="#0072B2", lw=3, alpha=0.5)
axhline(3, 0.0, rel_err_no_scale / 2.5, color="#0072B2", lw=5)

axhline(2.9, 0.0, rel_err_paper / 2.5, color="#009E73", lw=3, alpha=0.5)
axhline(2.9, 0.0, rel_err_paper_no_scale / 2.5, color="#009E73", lw=5)

# pi cc
# println("YOU ARE SUBTRACTING 1 FROM PICC ERR AT 7 GEV2")
rel_err = ( err(res["pi_cc"][10]) ) / value(res["pi_cc"][10]) * 100 
rel_err_no_scale = 0.00015078593060933852 / value(res["pi_cc"][10]) * 100 

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
sizes = [64.8, 18.1, 15.7]
push!(sizes, 100-sum(sizes))
explode = fill(0.03,4)
colors = ["#E69F00", "#56B4E9", "#009E73", "#CC79A7"]
pie(sizes, explode=explode, colors=colors)
legend(labels, fontsize=8, bbox_to_anchor=(0.5,0.9), loc="upper right")

subplot(132) # 5 GeV^2
title(L"$ \mathit{\bar{\Pi}}^{(3,3)}(Q^2 = 5 \ \mathrm{GeV}^2)$", fontsize=12)
labels = ("Systematics",  string("Stat. +",L"$\chi$", "-cont. limit"), "Scale setting", "Others")
sizes = [69.5, 17.1, 12.4]
push!(sizes, 100-sum(sizes))
explode = fill(0.03,4)
colors = ["#E69F00", "#56B4E9", "#009E73", "#CC79A7"]
pie(sizes, explode=explode, colors=colors)
legend(labels, fontsize=8, bbox_to_anchor=(0.5,0.9), loc="upper right")

subplot(133) # 9 GeV^2
title(L"$ \mathit{\bar{\Pi}}^{(3,3)}(Q^2 = 9 \ \mathrm{GeV}^2)$", fontsize=12)
labels = ("Systematics",  string("Stat. +",L"$\chi$", "-cont. limit"), "Scale setting", "Others")
sizes = [65.5, 21.5, 11.8]
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