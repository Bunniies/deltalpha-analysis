using PyPlot, LaTeXStrings
using Statistics
using ADerrors
import ADerrors:err
using OrderedCollections
using DelimitedFiles
#================ SET UP PLOT PARAMETERS ==================#

rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
rcParams["text.usetex"] =  true
rcParams["mathtext.fontset"]  = "cm"
rcParams["font.size"] =20
rcParams["axes.labelsize"] =26
rcParams["axes.titlesize"] = 22
plt.rc("text", usetex=true) # set to true if a LaTeX installation is present

path_plot = "/Users/alessandroconigli/MyDrive/postdoc-mainz/projects/deltalpha/plots/running"


path_kohtaroh = "/Users/alessandroconigli/MyDrive/postdoc-mainz/projects/deltalpha/running_jegerlhener/kohtaroh/pqcd_adler_jeg.dat"
path_pQCD_bfmom = "/Users/alessandroconigli/MyDrive/postdoc-mainz/projects/deltalpha/running_jegerlhener/test_pQCDAdler/using_bf-mom/pQCDAdler.dat53"
path_pQCD_msbar = "/Users/alessandroconigli/MyDrive/postdoc-mainz/projects/deltalpha/running_jegerlhener/test_pQCDAdler/using_msbar/pQCDAdler.dat53"
path_adlerpy = "/Users/alessandroconigli/MyDrive/postdoc-mainz/projects/deltalpha/running_jegerlhener/running_AdlerPy/using_my_input_values/adlerPy.dat"
path_adlerpy_mpole = "/Users/alessandroconigli/MyDrive/postdoc-mainz/projects/deltalpha/running_jegerlhener/running_AdlerPy/using_my_input_values/adlerPy_mpole.dat"


fname_ko = readdlm(path_kohtaroh, comments=true )
fname_pQCD_bfmom = readdlm(path_pQCD_bfmom, comments=true )
fname_pQCD_msbar = readdlm(path_pQCD_msbar, comments=true )
fname_adlerpy = readdlm(path_adlerpy, comments=true)
fname_adlerpy_mpole = readdlm(path_adlerpy_mpole, comments=true)
##

function read_pQCD_datfile(path; N=400, errinp=true, pl=false)
    fname = readdlm(path, comments=true)
    
    qval = (-fname[1:N,1])
    mean_val = (fname[1:N,3])

    err_npc = (fname[N+1:2N, 3] .- fname[2*N+1:end, 3]) #./2
    tot_err = err_npc

    if errinp
        path_alpha_err = joinpath(splitpath(path)[1:end-1])
        fname_pl = readdlm(joinpath(path_alpha_err,"data_err_alpha/pQCDAdler.dap52"), comments=true)
        fname_min = readdlm(joinpath(path_alpha_err,"data_err_alpha/pQCDAdler.dam52") , comments=true)

        err_input = (fname_pl[1:N,3] .- fname_min[1:N,3]) ./ 2
        tot_err = sqrt.(tot_err.^2 .+ err_input.^2) 
    end

    if pl
        errorbar(qval, mean_val, yerr=tot_err, fmt="s")
        
        display(gcf())
        close("all")
    end

    ad_res = []
    for k in eachindex(mean_val)
        push!(ad_res, uwreal([mean_val[k], tot_err[k]], "Dfunc Adpy"))
    end

    # return reverse(qval), reverse(ad_res)
    return qval, ad_res
end


# qval, val_mean, tot_err = read_pQCD_datfile(path_pQCD_msbar)
##
fig = figure(figsize=(10,8))

# Kohtaroh
Q_ko = fname_ko[:,2]
d_ko = fname_ko[:,3]
err_ko = fname_ko[:,6]
plot(Q_ko, d_ko, color="#FFC20A")
fill_between(Q_ko, d_ko .- err_ko, d_ko.+err_ko, color="#FFC20A", label="Kohtaroh", alpha=0.8)

# pQCD BF-MOM
Q_pqcd, d_pqcd_bfmom = read_pQCD_datfile(path_pQCD_bfmom)
uwerr.(d_pqcd_bfmom)
plot(Q_pqcd, value.(d_pqcd_bfmom),  color="royalblue")# color="#005AB5")
fill_between(Q_pqcd, value.(d_pqcd_bfmom) .- err.(d_pqcd_bfmom), value.(d_pqcd_bfmom) .+ err.(d_pqcd_bfmom), label="pQCD (BF-MOM)", color="royalblue", alpha=0.5)

# pQCD MSbar
Q_pqcd, d_pqcd_msbar = read_pQCD_datfile(path_pQCD_msbar)
uwerr.(d_pqcd_msbar)
# d_pqcd_msbar = fname_pQCD_msbar[1:400,2]
plot(Q_pqcd, value.(d_pqcd_msbar), color="black")# color="#005AB5")
fill_between(Q_pqcd, value.(d_pqcd_msbar) .- err.(d_pqcd_msbar) , value.(d_pqcd_msbar) .+ err.(d_pqcd_msbar),  label="pQCD (MSbar)", color="black", alpha=0.7)

# AdlerPy MSbar
Q_adlerpy = fname_adlerpy[:,2]
d_adlerpy = fname_adlerpy[:,3]
err_adlerpy = fname_adlerpy[:,4]
plot(Q_adlerpy, d_adlerpy, color="#D55E00")
fill_between(Q_adlerpy, d_adlerpy .- err_adlerpy, d_adlerpy.+err_adlerpy, color="#D55E00", label="AdlerPy (MSbar) ", alpha=0.8)

# AdlerPy MPOLE
d_adlerpy_mpole = fname_adlerpy_mpole[:,3]
err_adlerpy_mpole = fname_adlerpy_mpole[:,4]
plot(Q_adlerpy, d_adlerpy_mpole, color="#009E73")
fill_between(Q_adlerpy, d_adlerpy_mpole .- err_adlerpy_mpole, d_adlerpy_mpole .+ err_adlerpy_mpole, color="#009E73", label="AdlerPy (pole mass)", alpha=0.8)

# xlim(2.0,100) # 0
xlim(2,20) # 1
# xlim(15,100) # 2
# ylim(3.7,3.85) # 2 
# ylim(2.8,3.85) # 0 

legend(loc="lower right", ncol=2)
xlabel(L"$Q \ \mathrm{[GeV]}$")
ylabel(L"$D(Q)$")

# Cropped plot
fig.add_axes([0.3,0.42,0.657,0.25])

plot(Q_pqcd, value.(d_pqcd_bfmom), label="pQCD (BF-MOM)", color="royalblue")
fill_between(Q_pqcd, value.(d_pqcd_bfmom) .- err.(d_pqcd_bfmom), value.(d_pqcd_bfmom) .+ err.(d_pqcd_bfmom), color="royalblue", alpha=0.5)

plot(Q_pqcd, value.(d_pqcd_msbar), label="pQCD (MSbar)", color="black")# color="#005AB5")
fill_between(Q_pqcd, value.(d_pqcd_msbar) .- err.(d_pqcd_msbar), value.(d_pqcd_msbar) .+ err.(d_pqcd_msbar), color="black", alpha=0.5)

fill_between(Q_ko, d_ko .- err_ko, d_ko.+err_ko, color="#FFC20A", label="Kohtaroh", alpha=0.8)
plot(Q_adlerpy, d_adlerpy, color="#D55E00")
fill_between(Q_adlerpy, d_adlerpy .- err_adlerpy, d_adlerpy.+err_adlerpy, color="#D55E00", label="AdlerPy (MSbar)", alpha=0.8)
fill_between(Q_adlerpy, d_adlerpy_mpole .- err_adlerpy_mpole, d_adlerpy_mpole.+err_adlerpy_mpole, color="#009E73", label="AdlerPy (pole masses)", alpha=0.8)
xlim(6,20)
ylim(3.5,3.84)
# 

tight_layout()
display(fig)
savefig(joinpath(path_plot,"Adler_func_running_1.pdf"))
close("all")


## COMPARE RUNNING RESULTS (Tab 7 of 2203.08676 vs my results)
qval_km = [1,2,3,4,5,6,7]
KM_val = [0.023928, 0.022492, 0.021640, 0.021020, 0.020528, 0.020117, 0.019763]
KM_err = [223, 149, 129, 116, 107, 99, 93] .* 1e-6

qval_ac = [1,2,3,4,5,6,7,8,9]
AC_val_msbar = [0.024441, 0.022848, 0.021879, 0.021273, 0.020767, 0.020339, 0.019986, 0.019655, 0.019428]
AC_err_msbar = [194, 51, 35, 26, 23, 19, 18, 17, 17] .* 1e-6

fig = figure(figsize=(10,8))

errorbar(qval_km, KM_val, KM_err, fmt="s", label="KM", capsize=2, mfc="none", ms=8, color="royalblue")
errorbar(qval_ac, AC_val_msbar, AC_err_msbar, fmt="s", label="AC MSbar", capsize=2, mfc="none", ms=8, color="tomato")
xlabel(L"$Q_0^2 \ \mathrm{[GeV]^2}$")
ylabel(L"$\Delta\alpha_{\mathrm{had}}^{(5)}(-M_Z^2) - \Delta\alpha_{\mathrm{had}}^{(5)}(-Q_0^2)$")
legend()
tight_layout()
display(fig)
savefig(joinpath(path_plot,"running_comparison.pdf"))

close("all")
