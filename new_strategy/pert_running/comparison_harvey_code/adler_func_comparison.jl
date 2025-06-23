using HVPobs
using DelimitedFiles
using PyPlot
using LaTeXStrings
using ADerrors
import ADerrors:err
using Interpolations
using Statistics

include("../../../Running/tools_running.jl")


#================ SET UP PLOT PARAMETERS ==================#

rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
rcParams["text.usetex"] =  true
rcParams["mathtext.fontset"]  = "cm"
rcParams["font.size"] =20
rcParams["axes.labelsize"] =26
rcParams["axes.titlesize"] = 22
plt.rc("text", usetex=true) # set to true if a LaTeX installation is present

path_plot = "/Users/alessandroconigli/MyDrive/postdoc-mainz/projects/deltalpha/plots/new_strategy/pert_running/adler_func_comparison"

# pQCD
path_pQCD1_bfmom = "/Users/alessandroconigli/MyDrive/postdoc-mainz/projects/deltalpha/running_jegerlhener/test_pQCDAdler/using_bf_mom_pqcd1/pQCDAdler.dat53"
path_pQCD_msbar = "/Users/alessandroconigli/MyDrive/postdoc-mainz/projects/deltalpha/running_jegerlhener/test_pQCDAdler/using_msbar/pQCDAdler.dat53"

# pQCD modified by Dalibor
path_pQCD1_modified = "/Users/alessandroconigli/Desktop/postdoc-mainz/projects/deltalpha/running_jegerlhener/test_pQCDAdler/dalibor_pQCD1_v1_bfmom/pQCDAdler.dat53"

# AdlerPy my input
path_adpy_my_input = "/Users/alessandroconigli/MyDrive/postdoc-mainz/projects/deltalpha/running_jegerlhener/running_AdlerPy/using_my_input_values/log_spaced_mom_correct_quark_cond_sign/adlerPy.dat" # with bootstrap
path_adpy_my_input_4loops = "/Users/alessandroconigli/MyDrive/postdoc-mainz/projects/deltalpha/running_jegerlhener/running_AdlerPy/using_my_input_values/log_spaced_mom_correct_quark_cond_sign/adlerPy_4loops.dat" # with bootstrap
# path_adpy_my_input = "/Users/alessandroconigli/MyDrive/postdoc-mainz/projects/deltalpha/running_jegerlhener/running_AdlerPy/using_my_input_values/log_spaced_mom_no_bootstrap/adlerPy.dat" # no bootstrap

# AdlerPy km input
path_adpy_km_input = "/Users/alessandroconigli/MyDrive/postdoc-mainz/projects/deltalpha/running_jegerlhener/running_AdlerPy/using_km_input_values/adlerPy.dat"
path_adpy_km_input_4loops = "/Users/alessandroconigli/MyDrive/postdoc-mainz/projects/deltalpha/running_jegerlhener/running_AdlerPy/using_km_input_values/adlerPy_4loops.dat"

# Kohtaroh's data
path_adl_km = "/Users/alessandroconigli/MyDrive/postdoc-mainz/projects/deltalpha/running_jegerlhener/kohtaroh/pqcd_adler_jeg.dat"

## functions 
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
        push!(ad_res, uwreal([mean_val[k], tot_err[k]], "Dfunc pQCD"))
    end

    return reverse(qval), reverse(ad_res)
end

function read_pQCD1_datfile(path; N=400)
    fname = readdlm(path, comments=true)
    println(path)
    qval= -fname[1:N,1]
    mean_val = fname[1:N,3]

    path_err = joinpath(splitpath(path)[1:end-1])
    fname_pl = readdlm(joinpath(path_err,"data_err_alpha/pQCDAdler.dap52"), comments=true)
    fname_min = readdlm(joinpath(path_err,"data_err_alpha/pQCDAdler.dam52") , comments=true)
    err_input = (fname_pl[1:N,2] .- fname_min[1:N,2]) ./ 2

    fname_pl = readdlm(joinpath(path_err,"data_err_ope/pQCDAdler.dat52p"), comments=true)
    fname_min = readdlm(joinpath(path_err,"data_err_ope/pQCDAdler.dat52m") , comments=true)
    err_ope = (fname_pl[1:N,3] .- fname_min[1:N,3]) ./ 2
    
    tot_err = sqrt.(err_input.^2 + err_ope.^2)
    
    ad_res = []
    for k in eachindex(mean_val)
        push!(ad_res, uwreal([mean_val[k], tot_err[k]], "Dfunc pQCD1"))
    end
    return reverse(qval), reverse(ad_res)
end


function read_adlerpy_datfile(path; path_4loops::Union{String,Nothing}=nothing)
    fname = readdlm(path, comments=true)

    qval = fname[:,2]
    d_val = fname[:,3]
    d_err = fname[:,4]

    if !isnothing(path_4loops)
        fname = readdlm(path_4loops, comments=true)
        d_val_4loops = fname[:,3]
        d_err_4loops = d_val .- d_val_4loops
        d_err = sqrt.(d_err.^2 .+ d_err_4loops.^2)
    end

    ad_res = []
    for k in eachindex(qval)
        push!(ad_res, uwreal([d_val[k], d_err[k]], "Dfunc Adpy"))
    end
    return qval, ad_res
end


## reading data

qval_pQCD1_bfmom, adler_pQCD1_bfmom = read_pQCD1_datfile(path_pQCD1_bfmom)
qval_pQCD_msbar, adler_pQCD_msbar = read_pQCD_datfile(path_pQCD_msbar)

qval_adpy_my_input, adler_adpy_my_input = read_adlerpy_datfile(path_adpy_my_input, path_4loops=path_adpy_my_input_4loops)
qval_adpy_km_input, adler_adpy_km_input = read_adlerpy_datfile(path_adpy_km_input, path_4loops=path_adpy_km_input_4loops)

qval_pQCD_modified, adler_pQCD1_modified = read_pQCD1_datfile(path_pQCD1_modified)

# reading Kohtaroh's data
fname_adl_km= readdlm(path_adl_km, comments=true)
Q_val_adl_km = reverse(fname_adl_km[:,2])
dfunc_val_adl_km = reverse(fname_adl_km[:,3])
dfunc_err_adl_km = reverse(fname_adl_km[:,6])

dfunc_adl_km = Vector{uwreal}(undef, length(Q_val_adl_km))
for k in eachindex(Q_val_adl_km)
    dfunc_adl_km[k] =uwreal([dfunc_val_adl_km[k], dfunc_err_adl_km[k]], "Dfunc AdPy")
end

## PLOT ADLER FUNCTION

fig = figure(figsize=(10,8))

# AdlerPy KM input
#uwerr.(adler_adpy_km_input)
#fill_between(qval_adpy_km_input, value.(adler_adpy_km_input) .- err.(adler_adpy_km_input),  value.(adler_adpy_km_input) .+ err.(adler_adpy_km_input), alpha=0.5, color="forestgreen", label="AdlerPy (KM input)")

# pQCD1 BF-MOM
uwerr.(adler_pQCD1_bfmom)
fill_between(qval_pQCD1_bfmom, value.(adler_pQCD1_bfmom) .- err.(adler_pQCD1_bfmom),  value.(adler_pQCD1_bfmom) .+ err.(adler_pQCD1_bfmom), alpha=0.5, color="royalblue", label="pQCD1")

# pQCD MSbar
#uwerr.(adler_pQCD_msbar)
#fill_between(qval_pQCD_msbar, value.(adler_pQCD_msbar) .- err.(adler_pQCD_msbar),  value.(adler_pQCD_msbar) .+ err.(adler_pQCD_msbar), alpha=0.8, color="black", label="pQCD (MSbar)")

# AdlerPy my input
uwerr.(adler_adpy_my_input)
fill_between(qval_adpy_my_input, value.(adler_adpy_my_input) .- err.(adler_adpy_my_input),  value.(adler_adpy_my_input) .+ err.(adler_adpy_my_input), alpha=0.5, color="orange", label="AdlerPy ")

# pQCD1 modified by dalibot
uwerr.(adler_pQCD1_modified)
fill_between(qval_pQCD_modified, value.(adler_pQCD1_modified) .- err.(adler_pQCD1_modified),  value.(adler_pQCD1_modified) .+ err.(adler_pQCD1_modified), alpha=0.5, color="forestgreen", label="pQCD1 modified")

# Kohtaroh
uwerr.(dfunc_adl_km)
fill_between(Q_val_adl_km, value.(dfunc_adl_km) .- err.(dfunc_adl_km),  value.(dfunc_adl_km) .+ err.(dfunc_adl_km), alpha=0.5, color="red", label="Kohtaroh")


legend(loc="lower right", ncol=2)
xlabel(L"$Q \ \mathrm{[GeV]}$")
ylabel(L"$D(Q)$")

ylim(3.7,3.85) # 0 
# ylim(2.5,3.8)
xlim(10,100)

tight_layout()
display(fig)
savefig(joinpath(path_plot,"Adler_func_running_10_100.pdf"))
close("all")