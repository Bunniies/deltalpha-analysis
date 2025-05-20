using HVPobs
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




Qgev = [ 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0]
my_prec = [
    0.5704281186341132,
    0.5225064969006548,
    0.5140034666017413,
    0.5083762017222127,
    0.5040526668289174,
    0.49809619710809033,
    0.49181913431481356
]

dalpha_2022_val = [ 0.003864, 0.00521, 0.00605, 0.00666, 0.00716, 0.00757, 0.00793]
dalpha_2022_err = [ 0.000032, 0.00004, 0.00007, 0.00009, 0.00009, 0.00009, 0.00009]
dalpha_2022_err_myprec = dalpha_2022_val .* my_prec ./ 100

adPy_val = [0.023945, 0.022565, 0.021729, 0.021117, 0.020628, 0.020220, 0.019868]
adPy_err = [0.000479, 0.000147, 0.000085, 0.000062, 0.000051, 0.000044, 0.000040]

pQCD_val = [0.023926, 0.022489, 0.021638, 0.021018, 0.020525, 0.020114, 0.019760]
pQCD_err = [0.000223, 0.000148, 0.000128, 0.000116, 0.000106, 0.000099, 0.000092]


dalpha_2022_res = Vector{uwreal}(undef, length(Qgev))
dalpha_2022_res_myprec = Vector{uwreal}(undef, length(Qgev))
adPy_res =  Vector{uwreal}(undef, length(Qgev))
pQCD_res =  Vector{uwreal}(undef, length(Qgev))

for k in eachindex(Qgev)
    dalpha_2022_res[k] = uwreal([dalpha_2022_val[k], dalpha_2022_err[k]], "2022 res")
    dalpha_2022_res_myprec[k] = uwreal([dalpha_2022_val[k], dalpha_2022_err_myprec[k]], "2022 res my prec")
    adPy_res[k] = uwreal([adPy_val[k], adPy_err[k]], "Adpy")
    pQCD_res[k] = uwreal([pQCD_val[k], pQCD_err[k]], "pQCD")
end

pQCD_anal_continuation = uwreal([0.000045, 0.000002], "anal_continuation")

dalpha_fin_2022 = dalpha_2022_res .+ pQCD_res .+ pQCD_anal_continuation; uwerr.(dalpha_fin_2022); dalpha_fin_2022
dalpha_fin_myprec_pQCD = dalpha_2022_res_myprec .+ pQCD_res .+ pQCD_anal_continuation; uwerr.(dalpha_fin_myprec_pQCD)
dalpha_fin_myprec_adPy = dalpha_2022_res_myprec .+ adPy_res .+ pQCD_anal_continuation; uwerr.(dalpha_fin_myprec_adPy)

##
fig = figure()
fill_between(Qgev, value.(dalpha_fin_2022) .- err.(dalpha_fin_2022), value.(dalpha_fin_2022) .+ err.(dalpha_fin_2022), color="royalblue", alpha=0.5, label="Mainz 2022")
# fill_between(Qgev, value.(dalpha_fin_myprec_pQCD) .- err.(dalpha_fin_myprec_pQCD), value.(dalpha_fin_myprec_pQCD) .+ err.(dalpha_fin_myprec_pQCD), color="gold", alpha=0.2, label="2025 precision/PQCD")
# fill_between(Qgev, value.(dalpha_fin_myprec_adPy) .- err.(dalpha_fin_myprec_adPy), value.(dalpha_fin_myprec_adPy) .+ err.(dalpha_fin_myprec_adPy), color="tomato", alpha=0.2, label="2025 precision/AdPy")
# my precision with pQCD
plot(Qgev, value.(dalpha_fin_myprec_pQCD) .- err.(dalpha_fin_myprec_pQCD), ls="--", color="black", label="2025 precision/pQCD")
plot(Qgev, value.(dalpha_fin_myprec_pQCD) .+ err.(dalpha_fin_myprec_pQCD), ls="--", color="black")
# my precision with AdPy
plot(Qgev, value.(dalpha_fin_myprec_adPy) .- err.(dalpha_fin_myprec_adPy), ls="--", color="red", label="2025 precision/AdPy")
plot(Qgev, value.(dalpha_fin_myprec_adPy) .+ err.(dalpha_fin_myprec_adPy), ls="--", color="red")


legend()
ylabel(L"$\Delta\alpha_{\mathrm{had}}^{(5)}(M_Z^2)$")
xlabel(L"$Q^2 \ [\mathrm{GeV}^2]$")
tight_layout()
display(fig)
savefig(joinpath(path_plot, "my_precision_vs_2022_data.pdf"))
close("all")

## Pie Chart with  2022 results, 2025 precision, AdPy and pQCD for perturbative running
fig = figure(figsize=(10,8))
explode= fill(0.03,3)

labels = ("Lattice 2022", "Pert. pQCD", "Others")


subplot(331)# Q^2 = 3 GeV^2
title(L"$Q^2=3 \ \mathrm{GeV}^2$")
sizes = [23.02, 76.96]
push!(sizes, 100-sum(sizes))
colors = ["#E69F00", "#56B4E9", "#009E73"]
pie(sizes, explode=explode, colors=colors)
legend(labels, fontsize=8, bbox_to_anchor=(0.5,0.9), loc="upper right")

subplot(332)# Q^2 = 5 GeV^2
title(L"$Q^2=5 \ \mathrm{GeV}^2$")
sizes = [41.88, 58.10]
push!(sizes, 100-sum(sizes))
colors = ["#E69F00", "#56B4E9", "#009E73"]
pie(sizes, explode=explode, colors=colors)
legend(labels, fontsize=8, bbox_to_anchor=(0.5,0.9), loc="upper right")

subplot(333)# Q^2 = 7 GeV^2
title(L"$Q^2=7 \ \mathrm{GeV}^2$")
sizes = [48.89, 51.09]
push!(sizes, 100-sum(sizes))
colors = ["#E69F00", "#56B4E9", "#009E73"]
pie(sizes, explode=explode, colors=colors)
legend(labels, fontsize=8, bbox_to_anchor=(0.5,0.9), loc="upper right")


labels = ("Lattice 2025", "Pert. pQCD", "Others")
subplot(334) # Q^2 = 3 GeV^2
# title(L"$Q^2=3 \ \mathrm{GeV}^2$")
sizes = [5.57, 94.4]
push!(sizes, 100-sum(sizes))
colors = ["#E69F00", "#56B4E9", "#009E73"]
pie(sizes, explode=explode, colors=colors)
legend(labels, fontsize=8, bbox_to_anchor=(0.5,0.9), loc="upper right")

subplot(335) # Q^2 = 5 GeV^2
# title(L"$Q^2=5 \ \mathrm{GeV}^2$")
title(L"$\mathrm{Lattice\ 2025\ precision\ with\ pQCD}$")

sizes = [10.38, 89.58]
push!(sizes, 100-sum(sizes))
colors = ["#E69F00", "#56B4E9", "#009E73"]
pie(sizes, explode=explode, colors=colors)
legend(labels, fontsize=8, bbox_to_anchor=(0.5,0.9), loc="upper right")

subplot(336) # Q^2 = 7 GeV^2
# title(L"$Q^2=7 \ \mathrm{GeV}^2$")
sizes = [15.23, 84.73]
push!(sizes, 100-sum(sizes))
colors = ["#E69F00", "#56B4E9", "#009E73"]
pie(sizes, explode=explode, colors=colors)
legend(labels, fontsize=8, bbox_to_anchor=(0.5,0.9), loc="upper right")

labels = ("Lattice 2025", "Pert. AdPy", "Others")

subplot(337) # Q^2 = 3 GeV^2
# title(L"$Q^2=3 \ \mathrm{GeV}^2$")
sizes = [11.8, 88.15]
push!(sizes, 100-sum(sizes))
colors = ["#E69F00", "#56B4E9", "#009E73"]
pie(sizes, explode=explode, colors=colors)
legend(labels, fontsize=8, bbox_to_anchor=(0.5,0.9), loc="upper right")

subplot(338) # Q^2 = 5 GeV^2
# title(L"$Q^2=5 \ \mathrm{GeV}^2$")
title(L"$\mathrm{Lattice\ 2025\ precision\ with\ AdlerPy}$")
sizes = [33.33, 66.56]
push!(sizes, 100-sum(sizes))
colors = ["#E69F00", "#56B4E9", "#009E73"]
pie(sizes, explode=explode, colors=colors)
legend(labels, fontsize=8, bbox_to_anchor=(0.5,0.9), loc="upper right")

subplot(339)
# title(L"$Q^2=7 \ \mathrm{GeV}^2$")
sizes = [48.67, 51.20]
push!(sizes, 100-sum(sizes))
colors = ["#E69F00", "#56B4E9", "#009E73"]
pie(sizes, explode=explode, colors=colors)
legend(labels, fontsize=8, bbox_to_anchor=(0.5,0.9), loc="upper right")

suptitle(L"$\mathrm{Lattice\ 2022\ with\ pQCD}$", fontsize=18)
tight_layout()
display(fig)
savefig(joinpath(path_plot, "pie_chart_2022_2025_prec_pQCD_AdPy.pdf"))
close("all")