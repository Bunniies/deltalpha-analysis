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


alpha0=1/137.035999180;
Mz = 91.1876

path_plot = "/Users/alessandroconigli/MyDrive/postdoc-mainz/projects/deltalpha/plots/new_strategy/pert_running"

# pQCD
path_pQCD1_bfmom = "/Users/alessandroconigli/MyDrive/postdoc-mainz/projects/deltalpha/running_jegerlhener/test_pQCDAdler/using_bf_mom_pqcd1/pQCDAdler.dat53"

# pQCD modified
path_pQCD1_modified = "/Users/alessandroconigli/Desktop/postdoc-mainz/projects/deltalpha/running_jegerlhener/test_pQCDAdler/dalibor_pQCD1_v1_bfmom/pQCDAdler.dat53"

# AdlerPy my input
# path_adpy_my_input = "/Users/alessandroconigli/MyDrive/postdoc-mainz/projects/deltalpha/running_jegerlhener/running_AdlerPy/using_my_input_values/log_spaced_mom_correct_quark_cond_sign/adlerPy.dat" # with bootstrap
# path_adpy_my_input_4loops = "/Users/alessandroconigli/MyDrive/postdoc-mainz/projects/deltalpha/running_jegerlhener/running_AdlerPy/using_my_input_values/log_spaced_mom_correct_quark_cond_sign/adlerPy_4loops.dat" # with bootstrap
# path_adpy_my_input = "/Users/alessandroconigli/MyDrive/postdoc-mainz/projects/deltalpha/running_jegerlhener/running_AdlerPy/using_my_input_values/log_spaced_mom_no_bootstrap/adlerPy.dat" # no boostrap, error in quadrature
path_adpy_my_input = "/Users/alessandroconigli/MyDrive/postdoc-mainz/projects/deltalpha/running_jegerlhener/running_AdlerPy/using_my_input_values/log_spaced_mom_noBoot_maxErr/adlerPy.dat" # no bootstrap, maxErr
# Kohtaroh's values
path_adl_km = "/Users/alessandroconigli/MyDrive/postdoc-mainz/projects/deltalpha/running_jegerlhener/kohtaroh/pqcd_adler_jeg.dat"

## functions 
function read_pQCD_datfile(path; N=400, errinp=true, pl=false)
    fname = readdlm(path, comments=true)
    
    qval = (-fname[1:N,1])
    mean_val = (fname[1:N,3])

    err_npc = (fname[N+1:2N, 3] .- fname[2*N+1:end, 3]) #./2
    # err_npc = (mean_val .- fname[N+1:2*N, 3])
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
qvaltest, ad_res_test = read_pQCD1_datfile(path_pQCD1_bfmom)


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

function get_interp(qval, adler)
    itp = interpolate((qval,),  adler , Gridded(Linear()))
    return itp
end

function get_integral(data, dl) # apply trapezodial rule assuming data to be equally dl distant 
    return sum((data[1:end-1] + data[2:end])*0.5*dl)
end


function compute_integral(itp, a, b)
    
    int_range = collect(range(a,b,length=10000))
    int_res = sum(itp(int_range)[1:end] .* push!(diff(int_range), 0.0) ./ int_range[1:end])
    #int_res = sum(itp(int_range)[2:end] .* diff(int_range) ./ int_range[2:end])
    # int_res = sum((itp(int_range)[1:end-1]./int_range[1:end-1] + itp(int_range)[2:end]./int_range[2:end])/2 .* diff(int_range))

    return int_res 
end

## reading data

# pQCD1
qval_pQCD1_bfmom, adler_pQCD1_bfmom = read_pQCD1_datfile(path_pQCD1_bfmom)
# pQCD1 modified
qval_pQCD_modified, adler_pQCD1_modified = read_pQCD1_datfile(path_pQCD1_modified)
# AdPy 
# qval_adpy_my_input, adler_adpy_my_input = read_adlerpy_datfile(path_adpy_my_input, path_4loops=path_adpy_my_input_4loops) # with bootstrap
qval_adpy_my_input, adler_adpy_my_input = read_adlerpy_datfile(path_adpy_my_input) # no bootstrap

# Kohtaroh's values
fname_adl_km= readdlm(path_adl_km, comments=true)
Q_val_adl_km = reverse(fname_adl_km[:,2])
dfunc_val_adl_km = reverse(fname_adl_km[:,3])
dfunc_err_adl_km = reverse(fname_adl_km[:,6])

dfunc_adl_km = Vector{uwreal}(undef, length(Q_val_adl_km))
for k in eachindex(Q_val_adl_km)
    dfunc_adl_km[k] =uwreal([dfunc_val_adl_km[k], dfunc_err_adl_km[k]], "Dfunc AdPy")
end
##

## NEW INTEGRAL EVALUATION
store_res_km = []
store_res_adpy = []
store_res_pqcd1_bfmom = []
store_res_pqcd1_modified = []

Emin, Emax = extrema(Q_val_adl_km)
distq = mean(diff(log.(qval_pQCD1_bfmom.^2)))

for mom in range(1,9)# loop in q^2 GeV^2
    # ub=findall(x-> x - Mz >= 0, Q_val_adl_km)[1] -1
    # lb=findall(x-> x - sqrt(mom) >= 0, Q_val_adl_km)[1] -1
    ub=findall(x-> x - Mz >= 0, Q_val_adl_km)[1] -1
    lb=findall(x-> x - sqrt(mom) >= 0, Q_val_adl_km)[1] -1

    lb == 0 ? lb +=1 : lb
    println(Q_val_adl_km[ub], " ", Q_val_adl_km[lb])
    # KM data
    mres1 = get_integral(dfunc_adl_km[lb:ub], distq) * (alpha0/(3*pi))
    if mom == 1
        push!(store_res_km, mres1)
    else
        mres2 = get_integral(dfunc_adl_km[lb-1:ub], distq) * (alpha0/(3*pi))
        mres = mres2 + (mres1 - mres2) / (Q_val_adl_km[lb-1] - Q_val_adl_km[lb]) * (Q_val_adl_km[lb-1] - sqrt(mom))
        push!(store_res_km, mres )
    end
    # AdPy my input
    mres1 = get_integral(adler_adpy_my_input[lb:ub], distq) * (alpha0/(3*pi))
    if mom <= 1
        push!(store_res_adpy, mres1)
    else
        # mres2 = get_integral(itp_adpy_msbar[lb-1:ub], distq) * (alpha0/(3*pi))    
        mres2 = get_integral(adler_adpy_my_input[lb-1:ub], distq) * (alpha0/(3*pi))    
        mres = mres2 + (mres1 - mres2) / (Q_val_adl_km[lb-1] - Q_val_adl_km[lb]) * (Q_val_adl_km[lb-1] - sqrt(mom))
        push!(store_res_adpy, mres)
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
uwerr.(store_res_km)
uwerr.(store_res_pqcd1_bfmom)
uwerr.(store_res_pqcd1_modified)
uwerr.(store_res_adpy)


##
fig = figure(figsize=(10,8))
QQ = collect(1:9)
errorbar(QQ, value.(store_res_km), err.(store_res_km), fmt="s", ms=8, mfc="none", label="KM")
errorbar(QQ .+ 0.15, value.(store_res_pqcd1_bfmom), err.(store_res_pqcd1_bfmom), fmt="s", ms=8, mfc="none", label="pQCD1")
errorbar(QQ .+ 0.3, value.(store_res_adpy), err.(store_res_adpy), fmt="s", ms=8, mfc="none", label="AdPy (my input)")
errorbar(QQ .+ 0.45, value.(store_res_pqcd1_modified), err.(store_res_pqcd1_modified), fmt="s", ms=8, mfc="none", label="pQCD1 modified")
# errorbar(QQ .+ 0.6, value.(store_res_pqcd1_my_input), err.(store_res_pqcd1_my_input), fmt="s", ms=8, mfc="none", label="pQCD1(my input)")

xlabel(L"$Q_0^2 \ \mathrm{[GeV]^2}$")
ylabel(L"$\Delta\alpha_{\mathrm{had}}^{(5)}(-M_Z^2) - \Delta\alpha_{\mathrm{had}}^{(5)}(-Q_0^2)$")
legend()
tight_layout()
display(fig)
# savefig(joinpath(path_plot,"running_comparison.pdf"))
close("all")


## PLOT dalpha at Z pole using multiple codes for perturbative running

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

adPy_val = value.(store_res_adpy)
adPy_err = err.(store_res_adpy)

pQCD1_val = value.(store_res_pqcd1_bfmom)
pQCD1_err = err.(store_res_pqcd1_bfmom)

pQCD1_modified_val = value.(store_res_pqcd1_modified)
pQCD1_modified_err = err.(store_res_pqcd1_modified)

dalpha_2022_res = Vector{uwreal}(undef, length(Qgev))
dalpha_2022_res_myprec = Vector{uwreal}(undef, length(Qgev))
adPy_res =  Vector{uwreal}(undef, length(Qgev))
pQCD_res =  Vector{uwreal}(undef, length(Qgev))
pQCD_modified_res =  Vector{uwreal}(undef, length(Qgev))

for k in eachindex(Qgev)
    dalpha_2022_res[k] = uwreal([dalpha_2022_val[k], dalpha_2022_err[k]], "2022 res")
    dalpha_2022_res_myprec[k] = uwreal([dalpha_2022_val[k], dalpha_2022_err_myprec[k]], "2022 res my prec")
    adPy_res[k] = uwreal([adPy_val[k], adPy_err[k]], "Adpy")
    pQCD_res[k] = uwreal([pQCD1_val[k], pQCD1_err[k]], "pQCD")
    pQCD_modified_res[k] = uwreal([pQCD1_modified_val[k], pQCD1_modified_err[k]], "pQCD modified")
end

pQCD_anal_continuation = uwreal([0.000045, 0.000002], "anal_continuation")

dalpha_fin_2022 = dalpha_2022_res .+ pQCD_res .+ pQCD_anal_continuation; uwerr.(dalpha_fin_2022); dalpha_fin_2022
dalpha_fin_myprec_pQCD = dalpha_2022_res_myprec .+ pQCD_res .+ pQCD_anal_continuation; uwerr.(dalpha_fin_myprec_pQCD)
dalpha_fin_myprec_adPy = dalpha_2022_res_myprec .+ adPy_res .+ pQCD_anal_continuation; uwerr.(dalpha_fin_myprec_adPy)
dalpha_fin_myprec_pQCD_modified = dalpha_2022_res_myprec .+ pQCD_modified_res .+ pQCD_anal_continuation; uwerr.(dalpha_fin_myprec_pQCD_modified)

fig = figure(figsize=(10,8))
fill_between(Qgev, value.(dalpha_fin_2022) .- err.(dalpha_fin_2022), value.(dalpha_fin_2022) .+ err.(dalpha_fin_2022), color="royalblue", alpha=0.5, label="Mainz 2022")
# fill_between(Qgev, value.(dalpha_fin_myprec_pQCD) .- err.(dalpha_fin_myprec_pQCD), value.(dalpha_fin_myprec_pQCD) .+ err.(dalpha_fin_myprec_pQCD), color="gold", alpha=0.2, label="2025 precision/PQCD")
# fill_between(Qgev, value.(dalpha_fin_myprec_adPy) .- err.(dalpha_fin_myprec_adPy), value.(dalpha_fin_myprec_adPy) .+ err.(dalpha_fin_myprec_adPy), color="tomato", alpha=0.2, label="2025 precision/AdPy")
# my precision with pQCD
plot(Qgev, value.(dalpha_fin_myprec_pQCD) .- err.(dalpha_fin_myprec_pQCD), ls="--", color="black", label="2025 precision/pQCD")
plot(Qgev, value.(dalpha_fin_myprec_pQCD) .+ err.(dalpha_fin_myprec_pQCD), ls="--", color="black")
# my precision with AdPy
plot(Qgev, value.(dalpha_fin_myprec_adPy) .- err.(dalpha_fin_myprec_adPy), ls="-.", color="red", label="2025 precision/AdPy")
plot(Qgev, value.(dalpha_fin_myprec_adPy) .+ err.(dalpha_fin_myprec_adPy), ls="-.", color="red")
# my precision with pQCD modified
plot(Qgev, value.(dalpha_fin_myprec_pQCD_modified) .- err.(dalpha_fin_myprec_pQCD_modified), ls="-", color="forestgreen", label="2025 precision/pQCD modifed")
plot(Qgev, value.(dalpha_fin_myprec_pQCD_modified) .+ err.(dalpha_fin_myprec_pQCD_modified), ls="-", color="forestgreen")

legend()
ylabel(L"$\Delta\alpha_{\mathrm{had}}^{(5)}(M_Z^2)$")
xlabel(L"$Q^2 \ [\mathrm{GeV}^2]$")
xticks()
yticks()
tight_layout()
display(fig)
# savefig(joinpath(path_plot, "my_precision_vs_2022_data.pdf"))
close("all")


## saving adler function with error from AdlerPy for Dalibor
uwerr.(adler_adpy_my_input)
using DelimitedFiles

open("/Users/alessandroconigli/Desktop/adlerpy_adlerfunc.txt", "w") do io
    writedlm(io, ["# q     val      err"])
    [writedlm(io, [round(qval_adpy_my_input[k],digits=4) round(value(adler_adpy_my_input[k]),digits=5) err(adler_adpy_my_input[k])]) for k in eachindex(adler_adpy_my_input)]
end