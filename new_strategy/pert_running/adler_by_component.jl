using HVPobs
using DelimitedFiles
using PyPlot
using LaTeXStrings
using ADerrors
import ADerrors:err
using Interpolations
using Statistics

include("../../Running/tools_running.jl")

#================ SET UP PLOT PARAMETERS ==================#

rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
rcParams["text.usetex"] =  true
rcParams["mathtext.fontset"]  = "cm"
rcParams["font.size"] =20
rcParams["axes.labelsize"] =26
rcParams["axes.titlesize"] = 22
plt.rc("text", usetex=true) # set to true if a LaTeX installation is present

include("tools_adler.jl")

alpha0=1/137.035999180;
Mz = 91.1876

path_plot = "/Users/alessandroconigli/MyDrive/postdoc-mainz/projects/deltalpha/plots/new_strategy/pert_running/adler_func_comparison/"

# path AdlerPy
path_data_mu_as_trunc = "/Users/alessandroconigli/MyDrive/postdoc-mainz/projects/deltalpha/running_jegerlhener/running_AdlerPy/adler_per_comp/mu_as_trunc_err/"
path_data_std_trunc = "/Users/alessandroconigli/MyDrive/postdoc-mainz/projects/deltalpha/running_jegerlhener/running_AdlerPy/adler_per_comp/std_trunc_err/"
path_data_rodolfo_err = "/Users/alessandroconigli/MyDrive/postdoc-mainz/projects/deltalpha/running_jegerlhener/running_AdlerPy/adler_per_comp/rodolfo_err"
name_mu_as_trunc = "adlerPy_mutrunc_alpha_old_err.dat"
name_std_trunc = "adlerPy_trunc_err_included_alpha_old_err_QED_false.dat"
name_rodolfo_err = "adlerPy_rodolfo_old_alpha_err.dat"

# pQCD1 modified
path_pQCD1_modified = "/Users/alessandroconigli/Desktop/postdoc-mainz/projects/deltalpha/running_jegerlhener/test_pQCDAdler/dalibor_pQCD1_v1_bfmom/pQCDAdler.dat53"
# pQCD1
path_pQCD1_bfmom = "/Users/alessandroconigli/MyDrive/postdoc-mainz/projects/deltalpha/running_jegerlhener/test_pQCDAdler/using_bf_mom_pqcd1/pQCDAdler.dat53"

## 
# reading AdPy mu_trunc_err 
qval, adl_mu, adc_mu, adb_mu, adozi_mu = read_AdlerPy_by_component(path_data_mu_as_trunc, name_mu_as_trunc)
uwerr.(adl_mu); uwerr.(adc_mu); uwerr.(adb_mu); uwerr.(adozi_mu)
adtot_mu = adl_mu + adc_mu + adb_mu + adozi_mu; uwerr.(adtot_mu)

# reading AdPy std_strunc_err
qval, adl_std, adc_std, adb_std, adozi_std = read_AdlerPy_by_component(path_data_std_trunc, name_std_trunc)
uwerr.(adl_std); uwerr.(adc_std); uwerr.(adb_std); uwerr.(adozi_std)
adtot_std = adl_std + adc_std + adb_std + adozi_std; 
uwerr.(adtot_std)

# reading AdPy rodolfo err
qval, adl_rod, adc_rod, adb_rod, adozi_rod = read_AdlerPy_by_component(path_data_rodolfo_err, name_rodolfo_err)
uwerr.(adl_rod); uwerr.(adc_rod); uwerr.(adb_rod); uwerr.(adozi_rod)
adtot_rod = adl_rod + adc_rod + adb_std + adozi_std; 
uwerr.(adtot_rod)

# reading pQCD1 modified
qval_pQCD_modified, adler_pQCD1_modified = read_pQCD1_datfile(path_pQCD1_modified)
uwerr.(adler_pQCD1_modified)
# reading pQCD1
qval_pQCD1_bfmom, adler_pQCD1_bfmom = read_pQCD1_datfile(path_pQCD1_bfmom)
uwerr.(adler_pQCD1_bfmom)

## plotting std_trunc_err
fig = figure(figsize=(10,8))
errorbar(qval, value.(adl_std), err.(adl_std), fmt="s", label="light")
errorbar(qval, value.(adc_std), err.(adc_std), fmt="^", label="charm")
errorbar(qval, value.(adb_std), err.(adb_std), fmt="<", label="bottom")
errorbar(qval, value.(adozi_std), err.(adozi_std), fmt="o", label="disc")
errorbar(qval, value.(adtot_std), err.(adtot_std), fmt="d", color="black", label="tot")

legend(loc="upper left", ncol=1)
xlim(2,100)
xlabel(L"$Q \ [\mathrm{GeV}]$")
ylabel(L"$D(Q)$")
tight_layout()
display(fig)
# savefig(joinpath(path_plot, "adler_comp_std_trunc_err.pdf"))
close("all")


## plot contribution to the variance

var_tot = err.(adtot_std).^2
var_ratio_l = err.(adl_std).^2 ./  var_tot
var_ratio_c = err.(adc_std).^2 ./ var_tot
var_ratio_b = err.(adb_std).^2 ./ var_tot
var_ratio_ozi = err.(adozi_std).^2 ./ var_tot

fig = figure(figsize=(10,8))
fill_between(qval, 0.0, var_ratio_l, label="light")
fill_between(qval, var_ratio_l, var_ratio_l+var_ratio_c, label="charm")
fill_between(qval, var_ratio_l+var_ratio_c, var_ratio_l+var_ratio_c+var_ratio_b, label="bottom")
fill_between(qval, var_ratio_l+var_ratio_c+var_ratio_b, var_ratio_l+var_ratio_c+var_ratio_b+var_ratio_ozi, label="disc")

xlabel(L"$Q \ [\mathrm{GeV}]$")
ylabel(L"$\sigma_i^2 / \sum\sigma_i^2$")

legend()
tight_layout()
# savefig(joinpath(path_plot, "variance_contrib_by_component_charm_trunc.pdf"))
display(fig)
close("all")

## plotting mu_trunc_err
fig = figure(figsize=(10,8))
errorbar(qval, value.(adl_mu), err.(adl_mu), fmt="s", label="light")
errorbar(qval, value.(adc_mu), err.(adc_mu), fmt="^", label="charm")
errorbar(qval, value.(adb_mu), err.(adb_mu), fmt="<", label="bottom")
errorbar(qval, value.(adozi_mu), err.(adozi_mu), fmt="o", label="disc")
errorbar(qval, value.(adtot_mu), err.(adtot_mu), fmt="d", color="black", label="tot")

legend()
xlim(2,100)
xlabel(L"$Q \ [\mathrm{GeV}]$")
ylabel(L"$D(Q)$")
tight_layout()
display(fig)
savefig(joinpath(path_plot, "adler_comp_mu_err.pdf"))

close("all")

## plot adler from lattice vs adler from AdPy std_trunc_err

qlist = collect(range(sqrt(1), sqrt(8), length=100))

ad_lat = adler_lattice.(qlist); uwerr.(ad_lat)
fig = figure(figsize=(10,8))

fill_between(qval_pQCD1_bfmom.^2, value.(adler_pQCD1_bfmom).- err.(adler_pQCD1_bfmom), value.(adler_pQCD1_bfmom).+ err.(adler_pQCD1_bfmom), alpha=0.5, label="pQCD1", color="forestgreen")
fill_between(qval_pQCD_modified.^2, value.(adler_pQCD1_modified).- err.(adler_pQCD1_modified), value.(adler_pQCD1_modified).+ err.(adler_pQCD1_modified), alpha=0.5, color="black", label="pQCD1 modified")
fill_between(qlist.^2, value.(ad_lat) .- err.(ad_lat), value.(ad_lat) .+ err.(ad_lat), alpha=0.5, label="lattice", color="red")
fill_between(qval.^2, value.(adtot_std).-err.(adtot_std), value.(adtot_std).+err.(adtot_std), alpha=0.5, label="AdlerPy", color="royalblue")
# fill_between(qval.^2, value.(adtot_rod).-err.(adtot_rod), value.(adtot_rod).+err.(adtot_rod), alpha=0.5, label="AdlerPy", color="royalblue")

# xlim(1.4, 6.5)
xlim(1.5, 14)
# ylim(2.4, 3.1)
ylim(2.4, 3.4)
xlabel(L"Q^2 \ [\mathrm{GeV}^2]")
ylabel(L"$D(Q)$")
legend(loc="lower right")
tight_layout()
display(fig)
# savefig(joinpath(path_plot, "adler_lattice_vs_AdlerPy.pdf"))
# savefig(joinpath(path_plot, "adler_lattice_vs_AdlerPy_vs_pQCD.pdf"))
close("all")

####################
## compute integral
####################

Emin, Emax = extrema(qval)
distq = mean(diff(log.(qval.^2)))
# distq = mean(diff(log.(qval_pQCD1_bfmom.^2)))

res_adtot_std = []
res_adtot_mu = []
res_adtot_rodolfo = []
store_res_pqcd1_bfmom = []
store_res_pqcd1_modified = []

# qval_pQCD1_bfmom[393] = Mz
for mom in range(1,9)
    ub=findall(x-> x - Mz >= 0, qval)[1] -1
    lb=findall(x-> x - sqrt(mom) >= 0, qval)[1] -1

    lb == 0 ? lb +=1 : lb
    println(qval[ub], " ", qval[lb], " ub=", ub)

    # adpy std_trunc_err
    mres1 = get_integral(adtot_std[lb:ub], distq) * (alpha0/(3*pi))
    if mom == 1
        push!(res_adtot_std, mres1)
    else
        mres2 = get_integral(adtot_std[lb-1:ub], distq) * (alpha0/(3*pi))
        mres = mres2 + (mres1 - mres2) / (qval[lb-1] - qval[lb]) * (qval[lb-1] - sqrt(mom))
        push!(res_adtot_std, mres)
    end
    # adpy mu_trunc_err
    mres1 = get_integral(adtot_mu[lb:ub], distq) * (alpha0 / (3pi))
    if mom == 1
        push!(res_adtot_mu, mres1)
    else
        mres2 = get_integral(adtot_mu[lb-1:ub], distq) * (alpha0/ (3pi))
        mres = mres2 + (mres1 - mres2) / (qval[lb-1] - qval[lb]) * (qval[lb-1] - sqrt(mom))
        push!(res_adtot_mu, mres)
    end
    # adpy rodolfo
    mres1 = get_integral(adtot_rod[lb:ub], distq) * (alpha0 / (3pi))
    if mom == 1
        push!(res_adtot_rodolfo, mres1)
    else
        mres2 = get_integral(adtot_rod[lb-1:ub], distq) * (alpha0/ (3pi))
        mres = mres2 + (mres1 - mres2) / (qval[lb-1] - qval[lb]) * (qval[lb-1] - sqrt(mom))
        push!(res_adtot_rodolfo, mres)
    end
    # pqcd1 BF-MOM
    mres1 = get_integral(adler_pQCD1_bfmom[lb:ub], distq) * (alpha0 / (3pi))
    if mom ==1 
        push!(store_res_pqcd1_bfmom, mres1)
    else
        mres2 = get_integral(adler_pQCD1_bfmom[lb-1:ub], distq) * (alpha0 / (3pi))
        mres = mres2 + (mres1 - mres2) / (qval_pQCD1_bfmom[lb-1] - qval_pQCD1_bfmom[lb]) * (qval_pQCD1_bfmom[lb-1] - sqrt(mom))
        push!(store_res_pqcd1_bfmom, mres)
    end
    # pqcd1 BF-MOM modified by Dalibor
    mres1 = get_integral(adler_pQCD1_modified[lb:ub], distq) * (alpha0 / (3pi))
    if mom ==1 
        push!(store_res_pqcd1_modified, mres1)
    else
        mres2 = get_integral(adler_pQCD1_modified[lb-1:ub], distq) * (alpha0 / (3pi))
        mres = mres2 + (mres1 - mres2) / (qval_pQCD1_bfmom[lb-1] - qval_pQCD1_bfmom[lb]) * (qval_pQCD1_bfmom[lb-1] - sqrt(mom))
        push!(store_res_pqcd1_modified, mres)
    end
end 
uwerr.(res_adtot_mu)
uwerr.(res_adtot_std)
uwerr.(store_res_pqcd1_bfmom)
uwerr.(store_res_pqcd1_modified)
uwerr.(res_adtot_rodolfo)


## PLOT COMPARISON
fig = figure(figsize=(10,8))
QQ = collect(1:9)
errorbar(QQ, value.(store_res_pqcd1_bfmom), err.(store_res_pqcd1_bfmom), fmt="s", ms=8, mfc="none", label="KM")
errorbar(QQ.+0.2, value.(store_res_pqcd1_modified), err.(store_res_pqcd1_modified), fmt="s", ms=8, mfc="none", label="pQCD1 modified")
errorbar(QQ.+0.4, value.(res_adtot_std), err.(res_adtot_std), fmt="s", mfc="none", ms=8, label="AdlerPy")
xlabel(L"$Q_0^2 \ \mathrm{[GeV]^2}$")
ylabel(L"$\Delta\alpha_{\mathrm{had}}^{(5)}(-M_Z^2) - \Delta\alpha_{\mathrm{had}}^{(5)}(-Q_0^2)$")
legend()
tight_layout()
display(fig)
savefig(joinpath(path_plot,"running_comparison.pdf"))
close("all")