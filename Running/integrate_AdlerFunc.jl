using HVPobs
using DelimitedFiles
using PyPlot
using LaTeXStrings
using ADerrors
import ADerrors:err
using Interpolations
using Statistics

include("./tools_running.jl")

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

path_plot = "/Users/alessandroconigli/MyDrive/postdoc-mainz/projects/deltalpha/plots/running"

path_adpy_msbar = "/Users/alessandroconigli/MyDrive/postdoc-mainz/projects/deltalpha/running_jegerlhener/running_AdlerPy/using_my_input_values/log_spaced_mom/adlerPy.dat"
path_adpy_msbar_4loops = "/Users/alessandroconigli/MyDrive/postdoc-mainz/projects/deltalpha/running_jegerlhener/running_AdlerPy/using_my_input_values/log_spaced_mom/adlerPy_4loops.dat"
path_adpy_msbar_km_input = "/Users/alessandroconigli/MyDrive/postdoc-mainz/projects/deltalpha/running_jegerlhener/running_AdlerPy/using_km_input_values/adlerPy.dat"
path_adpy_msbar_4loops_km_input = "/Users/alessandroconigli/MyDrive/postdoc-mainz/projects/deltalpha/running_jegerlhener/running_AdlerPy/using_km_input_values/adlerPy_4loops.dat"

path_adl_km = "/Users/alessandroconigli/MyDrive/postdoc-mainz/projects/deltalpha/running_jegerlhener/kohtaroh/pqcd_adler_jeg.dat"
path_pQCD_bfmom = "/Users/alessandroconigli/MyDrive/postdoc-mainz/projects/deltalpha/running_jegerlhener/test_pQCDAdler/using_bf-mom/pQCDAdler.dat53"
path_pQCD_msbar = "/Users/alessandroconigli/MyDrive/postdoc-mainz/projects/deltalpha/running_jegerlhener/test_pQCDAdler/using_msbar/pQCDAdler.dat53"
path_pQCD1_bfmom = "/Users/alessandroconigli/MyDrive/postdoc-mainz/projects/deltalpha/running_jegerlhener/test_pQCDAdler/using_bf_mom_pqcd1/pQCDAdler.dat53"
path_pQCD1_my_input = "/Users/alessandroconigli/MyDrive/postdoc-mainz/projects/deltalpha/running_jegerlhener/test_pQCDAdler/using_bf_mom_pqcd1_newInput/pQCDAdler.dat53"
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

qval_adpy, dfunc_adpy_msbar = read_adlerpy_datfile(path_adpy_msbar, path_4loops=path_adpy_msbar_4loops)  
qval_adpy, dfunc_adpy_msbar_km_input = read_adlerpy_datfile(path_adpy_msbar_km_input, path_4loops=path_adpy_msbar_4loops_km_input)  
qval_pqcd, dfunc_pqcd_bfom = read_pQCD_datfile(path_pQCD_bfmom)
qval_pqcd, dfunc_pqcd_msbar = read_pQCD_datfile(path_pQCD_msbar)
qval_pqcd, dfunc_pqcd1_bfom = read_pQCD1_datfile(path_pQCD1_bfmom)
qval_pqcd, dfunc_pqcd1_my_input = read_pQCD1_datfile(path_pQCD1_my_input)

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
store_res_km_myval = []
store_res_pqcd_msbar = []
store_res_pqcd_bfmom = []
store_res_adpy_msbar = []
store_res_adpy_msbar_km_input = []
store_res_pqcd1_bfmom = []
store_res_pqcd1_my_input = []

Emin, Emax = extrema(Q_val_adl_km)
distq = mean(diff(log.(qval_pqcd.^2)))

itp_adpy_msbar = get_interp(log.(qval_adpy.^2), dfunc_adpy_msbar)(log.(Q_val_adl_km.^2))
itp_adpy_msbar_km_input = get_interp(qval_adpy, dfunc_adpy_msbar_km_input)(Q_val_adl_km)
itp_km_myval = get_interp(log.(Q_val_adl_km.^2), dfunc_adl_km)

for mom in range(1,9)
    # ub=findall(x-> x - Mz >= 0, Q_val_adl_km)[1] -1
    # lb=findall(x-> x - sqrt(mom) >= 0, Q_val_adl_km)[1] -1
    ub=findall(x-> x - Mz >= 0, Q_val_adl_km)[1] -1
    lb=findall(x-> x - sqrt(mom) >= 0, Q_val_adl_km)[1] -1

    lb == 0 ? lb +=1 : lb
    # KM data
    mres1 = get_integral(dfunc_adl_km[lb:ub], distq) * (alpha0/(3*pi))
    if mom == 1
        push!(store_res_km, mres1)
    else
        mres2 = get_integral(dfunc_adl_km[lb-1:ub], distq) * (alpha0/(3*pi))
        mres = mres2 + (mres1 - mres2) / (Q_val_adl_km[lb-1] - Q_val_adl_km[lb]) * (Q_val_adl_km[lb-1] - sqrt(mom))
        push!(store_res_km, mres )
    end
    # pqcd BF-MOM
    mres1 = get_integral(dfunc_pqcd_bfom[lb:ub], distq) * (alpha0/(3*pi))
    if mom == 1
        push!(store_res_pqcd_bfmom, mres1)
    else
        mres2 = get_integral(dfunc_pqcd_bfom[lb-1:ub], distq) * (alpha0/(3*pi))
        mres = mres2 + (mres1 - mres2) / (Q_val_adl_km[lb-1] - Q_val_adl_km[lb]) * (Q_val_adl_km[lb-1] - sqrt(mom))
        push!(store_res_pqcd_bfmom, mres )
    end
    # pqcd MSBAR
    mres1 = get_integral(dfunc_pqcd_msbar[lb:ub], distq) * (alpha0/(3*pi))
    if mom ==1
        push!(store_res_pqcd_msbar, mres1)
    else
        mres2 = get_integral(dfunc_pqcd_msbar[lb-1:ub], distq) * (alpha0/(3*pi))
        mres = mres2 + (mres1 - mres2) / (Q_val_adl_km[lb-1] - Q_val_adl_km[lb]) * (Q_val_adl_km[lb-1] - sqrt(mom))
        push!(store_res_pqcd_msbar, mres)
    end
    # AdPy MSBAR
    # mres1 = get_integral(itp_adpy_msbar[lb:ub], distq) * (alpha0/(3*pi))
    mres1 = get_integral(dfunc_adpy_msbar[lb:ub], distq) * (alpha0/(3*pi))
    if mom <= 1
        push!(store_res_adpy_msbar, mres1)
    else
        # mres2 = get_integral(itp_adpy_msbar[lb-1:ub], distq) * (alpha0/(3*pi))    
        mres2 = get_integral(dfunc_adpy_msbar[lb-1:ub], distq) * (alpha0/(3*pi))    
        mres = mres2 + (mres1 - mres2) / (Q_val_adl_km[lb-1] - Q_val_adl_km[lb]) * (Q_val_adl_km[lb-1] - sqrt(mom))
        push!(store_res_adpy_msbar, mres)
    end
    # AdPy MSBAR KM input
    push!(store_res_adpy_msbar_km_input, get_integral(itp_adpy_msbar_km_input[lb:ub], distq) * (alpha0/(3*pi)))
    # pqcd1 BF-MOM
    push!(store_res_pqcd1_bfmom, get_integral(dfunc_pqcd1_bfom[lb:ub], distq) * (alpha0/(3*pi)))
    # pqcd1 BF-MOM my input
    push!(store_res_pqcd1_my_input, get_integral(dfunc_pqcd1_my_input[lb:ub], distq)* (alpha0/(3*pi)))
    # test my method with KM data
    val_to_int = itp_km_myval(range(log(mom), log(Mz^2), step=distq))
    push!(store_res_km_myval, get_integral(val_to_int, distq) * (alpha0/(3*pi)))
end
uwerr.(store_res_km)
uwerr.(store_res_pqcd_bfmom)
uwerr.(store_res_pqcd_msbar)
uwerr.(store_res_adpy_msbar)
uwerr.(store_res_adpy_msbar_km_input)
uwerr.(store_res_pqcd1_bfmom)
uwerr.(store_res_km_myval)
uwerr.(store_res_pqcd1_my_input)


##
fig = figure(figsize=(10,8))
QQ = collect(1:9)
errorbar(QQ, value.(store_res_km), err.(store_res_km), fmt="s", ms=8, mfc="none", label="KM")
errorbar(QQ .+ 0.15, value.(store_res_adpy_msbar), err.(store_res_adpy_msbar), fmt="s", ms=8, mfc="none", label="AdPy (my input)")
errorbar(QQ .+ 0.3 , value.(store_res_adpy_msbar_km_input), err.(store_res_adpy_msbar_km_input), fmt="s", ms=8, mfc="none", label="AdPy (KM input)")
#errorbar(QQ , value.(store_res_pqcd_msbar), err.(store_res_pqcd_msbar), fmt="s", ms=8, mfc="none", label="pQCD(MSbar)")
#errorbar(QQ .+ 0.3, value.(store_res_pqcd_bfmom), err.(store_res_pqcd_bfmom), fmt="s", ms=8, mfc="none", label="pQCD(BF-MOM)")
errorbar(QQ .+ 0.45, value.(store_res_pqcd1_bfmom), err.(store_res_pqcd1_bfmom), fmt="s", ms=8, mfc="none", label="pQCD1(BF-MOM)")
# errorbar(QQ .+ 0.6, value.(store_res_pqcd1_my_input), err.(store_res_pqcd1_my_input), fmt="s", ms=8, mfc="none", label="pQCD1(my input)")

xlabel(L"$Q_0^2 \ \mathrm{[GeV]^2}$")
ylabel(L"$\Delta\alpha_{\mathrm{had}}^{(5)}(-M_Z^2) - \Delta\alpha_{\mathrm{had}}^{(5)}(-Q_0^2)$")
legend()
tight_layout()
display(fig)
savefig(joinpath(path_plot,"running_comparison.pdf"))
close("all")

## TEST FOLLOWING  RGE eq in https://arxiv.org/pdf/2308.05740


function rge_dalpha_running(dalpha_n_q2, q2; n=3)

    delta_mc = uwreal([1.278, 0.006], "mc(mc)") - 1.270
    delta_mb = uwreal([4.17, 0.02], "mc(mb)") - 4.180
    delta_alphas = uwreal([0.1183, 0.0007], "dalha_s") - 0.1185
    gluon_cond = uwreal([0.01, 0.01], "gluon_cond")  
    quark_cond = uwreal([-0.003, 0.003], "quark_cond")  
    
    parametric_formula = (218.32 - 17.3*log(q2/4)
    +350*delta_alphas - 20*delta_mc -1.3*delta_mb
    -(1.8-48/q2^2)*gluon_cond + 219/q2^2*quark_cond
    ) * 1e-4

    parametric_formula_charm = (135.85 - 17.3*log(q2/4)
    +143*delta_alphas  -1.3*delta_mb
    -(-48/q2^2)*gluon_cond + 219/q2^2*quark_cond
    ) * 1e-4


    if n == 3
        return dalpha_n_q2 + parametric_formula
    elseif n==4
        return dalpha_n_q2 +  parametric_formula_charm
    else
        error("The value of n has to be 3 or 4, depending whether the charm contribution 
        is computed non-perturbatively (n=4) or not (n=3)  ")
    end
    
end

alpha0=1/137.035999180
isovec_mainz_4q2 = uwreal([0.0529,0.0009], "isovec_mainz") 
isoscal_mainz_4q2 = uwreal([0.0414,0.0007], "isoscal_mainz")
charm_mainz_4q2 = uwreal([0.01348,0.00024], "charm_mainz") 

dalpha_3_mainz =  4*pi*alpha0 *(isovec_mainz_4q2 + 1/3*isoscal_mainz_4q2)
dalpha_4_mainz = 4*pi*alpha0 *(isovec_mainz_4q2 + 1/3*isoscal_mainz_4q2 + 4/9*charm_mainz_4q2)

## comparison of dalpha using different perturbative evaluation and Mainz 22 lattice data
dalpha_zpole_RGE = rge_dalpha_running(dalpha_3_mainz, 4, n=3)
dalpha_zpole_km = dalpha_4_mainz + store_res_km[4] + uwreal([0.000045,0.000002], "run_Mz_minus_Mz")
dalpha_zpole_adpy_msbar = dalpha_4_mainz + store_res_adpy_msbar[4] + uwreal([0.000045,0.000002], "run_Mz_minus_Mz")
dalpha_zpole_adpy_msbar_km_in = dalpha_4_mainz + store_res_adpy_msbar_km_input[4] + uwreal([0.000045,0.000002], "run_Mz_minus_Mz")
dalpha_zpole_pqcd_bfmom = dalpha_4_mainz + store_res_pqcd_bfmom[4] + uwreal([0.000045,0.000002], "run_Mz_minus_Mz")
dalpha_zpole_pqcd_msbar = dalpha_4_mainz + store_res_pqcd_msbar[4] + uwreal([0.000045,0.000002], "run_Mz_minus_Mz")
dalpha_mainz_pub = uwreal([0.02773, 0.00015], "Mainz 22")
dalpha_zpole_pqcd1_bfmom = dalpha_4_mainz + store_res_pqcd1_bfmom[4] + uwreal([0.000045,0.000002], "run_Mz_minus_Mz")

uwerr(dalpha_zpole_RGE)
uwerr(dalpha_zpole_km)
uwerr(dalpha_zpole_adpy_msbar)
uwerr(dalpha_zpole_adpy_msbar_km_in)
uwerr(dalpha_zpole_pqcd_bfmom)
uwerr(dalpha_zpole_pqcd_msbar)
uwerr(dalpha_mainz_pub)
uwerr(dalpha_zpole_pqcd1_bfmom)


##
fig = figure(figsize=(10,8))

errorbar(value(dalpha_mainz_pub), 7, xerr=err(dalpha_mainz_pub), fmt="s", ms=10, mfc="none", lw=2, color="forestgreen", capsize=4)
fill_betweenx(collect(1:7), value(dalpha_mainz_pub) .- err(dalpha_mainz_pub), value(dalpha_mainz_pub) .+ err(dalpha_mainz_pub), alpha=0.3, color="forestgreen")

errorbar(value(dalpha_zpole_km), 6, xerr=err(dalpha_zpole_km), fmt="s", ms=10, mfc="none", color="black", capsize=2)
errorbar(value(dalpha_zpole_adpy_msbar), 5, xerr=err(dalpha_zpole_adpy_msbar), fmt="s", ms=10, mfc="none", color="black", capsize=2)
errorbar(value(dalpha_zpole_adpy_msbar_km_in), 4, xerr=err(dalpha_zpole_adpy_msbar_km_in), fmt="s", ms=10, mfc="none", color="black", capsize=2)
errorbar(value(dalpha_zpole_pqcd_bfmom), 3, xerr=err(dalpha_zpole_pqcd_bfmom), fmt="s", ms=10, mfc="none", color="black", capsize=2)
errorbar(value(dalpha_zpole_pqcd1_bfmom), 2, xerr=err(dalpha_zpole_pqcd1_bfmom), fmt="s", ms=10, mfc="none", color="black", capsize=2)

# errorbar(value(dalpha_zpole_pqcd_msbar), 2, xerr=err(dalpha_zpole_pqcd_msbar), fmt="s", ms=10, mfc="none", color="black", capsize=2)
errorbar(value(dalpha_zpole_RGE), 1, xerr=err(dalpha_zpole_RGE), fmt="s", ms=10, mfc="none", color="black", capsize=2)


xlabel(L"$\Delta\alpha_{\mathrm{had}}^{(5)}(M_Z^2)$")

# yticks(reverse(collect(1:7)), ["Mainz 22", "KM data", "AdPy (MSbar)", "AdPy (pole mass)", "pQCD (BF-MOM)", "pQCD (MSbar)", "RGE"])
yticks(reverse(collect(1:7)), ["Mainz 22", "KM data", "AdPy (my input)", "AdPy (KM input)", "pQCD (BF-MOM)", "pQCD1 (BF-MOM)", "RGE"])

tight_layout()
display(fig)
# savefig(joinpath(path_plot, "dalpha_zpole_from_Mainz22.pdf"))
close("all")

## TEST LOG INTEGRAL
function get_interp_new(qval, adler)
    itp = interpolate((log.(qval.^2),), (alpha0) / (3*pi)  * adler , Gridded(Linear()))
    return itp
end
function compute_integral_new(itp, a, b)
    
    int_range = collect(range(a,b,length=1000))
    int_res = sum((itp(int_range)[1:end-1] + itp(int_range)[2:end])/2 .* diff(int_range))
    return int_res 
end

function get_integral(data, dl)
    return sum((data[1:end-1] + data[2:end])*0.5*dl)
end

get_integral(reverse(dfunc_adl_km)[10:340], 0.023083559829514246)*(alpha0/(3*pi))
get_integral(dfunc_adl_km[61:391], 0.023083559829514246)*(alpha0/(3*pi))
get_integral(dfunc_adl_km[62:392], 0.023083559829514246)*(alpha0/(3*pi))


Emin, Emax = extrema(Q_val_adl_km)

using Statistics
distq = mean(diff(log.(Q_val_adl_km.^2)))
res = []
for mom in range(1,9)
    ub=findall(x-> x - Mz >= 0, Q_val_adl_km)[1] -1
    lb=findall(x-> x - sqrt(mom) >= 0, Q_val_adl_km)[1] -1
    lb == 0 ? lb +=1 : lb
    aux = get_integral(dfunc_adl_km[lb:ub], distq)*(alpha0/(3*pi))
    push!(res, aux)
    println(ub, " ", lb)
end
ref_paper = 0.021020
itp_km = get_interp(Q_val_adl_km, dfunc_adl_km)
int_km = compute_integral(itp_km, sqrt(4), Mz)

itp_km_new = get_interp_new(Q_val_adl_km, dfunc_adl_km)
int_km_new = (compute_integral_new(itp_km_new, log(4), log(Mz^2)))
uwerr.(int_km_new); int_km_new

## figure integral comparison
alpha0=1/137.035999180
isovec_mainz_4q2 = uwreal([0.0529,0.0009], "isovec_mainz") 
isoscal_mainz_4q2 = uwreal([0.0414,0.0007], "isoscal_mainz")
charm_mainz_4q2 = uwreal([0.01348,0.00024], "charm_mainz") 

dalpha_4_mainz = 4*pi*alpha0 *(isovec_mainz_4q2 + 1/3*isoscal_mainz_4q2 + 4/9*charm_mainz_4q2)

dalpha_mainz_pub = uwreal([0.02773, 0.00015], "Mainz 22")
dalpha_mainz_pub_4q2 = dalpha_4_mainz + uwreal([0.000045,0.000002], "run_Mz_minus_Mz") + uwreal([0.021020,0.000116], "test")
dalpha_km_int1 = dalpha_4_mainz + int_km +  uwreal([0.000045,0.000002], "run_Mz_minus_Mz")
dalpha_km_int2 = dalpha_4_mainz + int_km_new +  uwreal([0.000045,0.000002], "run_Mz_minus_Mz")

uwerr(dalpha_mainz_pub)
uwerr(dalpha_km_int1)
uwerr(dalpha_km_int2)

fig = figure(figsize=(10,8))
title("Integral comparison using KM data")

errorbar(value(dalpha_mainz_pub), 7, xerr=err(dalpha_mainz_pub), fmt="s", ms=10, mfc="none", lw=2, color="forestgreen", capsize=4)
fill_betweenx(collect(1:7), value(dalpha_mainz_pub) .- err(dalpha_mainz_pub), value(dalpha_mainz_pub) .+ err(dalpha_mainz_pub), alpha=0.3, color="forestgreen")
errorbar(value(dalpha_km_int1), 6, xerr=err(dalpha_km_int1), fmt="s", ms=10, mfc="none", color="black", capsize=2)
errorbar(value(dalpha_km_int2), 5, xerr=err(dalpha_km_int2), fmt="s", ms=10, mfc="none", color="black", capsize=2)

xlabel(L"$\Delta\alpha_{\mathrm{had}}^{(5)}(M_Z^2)$")
yticks(reverse(collect(5:7)), ["Mainz 22", "Standard int", "Log int"])
ylim(4,7.2)

tight_layout()
display(fig)
savefig(joinpath(path_plot, "dalpha_zpole_int_comparison.pdf"))
close("all")



## TEST FITTING TO PERFORM INTEGRAL

@. func(x,p) = p[1] + p[2]*x^2 + p[3]*x^3 + p[4]*x^4 + p[5]*x^5 + p[6]*x^6 + p[7]*x^7 + p[8]*x^8 + p[9]*x^9 

function make_model(;N=10)
    ff = (x,p) ->  sum([p[k].* x.^(k) for k in 1:N])
    return ff
end

ydata = (2*alpha0) / (3*pi) * dfunc_adl_km ./ Q_val_adl_km
uwerr.(ydata)

funz = make_model(N=10)
fit_test = fit_routine(funz, Q_val_adl_km[1:end], ydata[1:end], 10)

test = int_error(funz, sqrt(4), Mz, fit_test.param )

##
errorbar(Q_val_adl_km, value.(ydata), err.(ydata), fmt="s")
yy = funz(Q_val_adl_km, fit_test.param)
plot(Q_val_adl_km, value.(yy), lw=10)
display(gcf())
close("all")




## OLD INTEGRAL EVALUATION
import Base.zero
Base.zero(::Type{uwreal}) = uwreal(0.0)

itp_km = get_interp(Q_val_adl_km, dfunc_adl_km)
itp_pqcd_msbar = get_interp(qval_pqcd, dfunc_pqcd_msbar)
itp_pqcd_bfmom = get_interp(qval_pqcd, dfunc_pqcd_bfom)
itp_adpy_msbar = get_interp(qval_pqcd, dfunc_adpy_msbar)
itp_adpy_mpole = get_interp(qval_pqcd, dfunc_adpy_mpole)
itp_pqcd1_bfmom = get_interp(qval_pqcd, dfunc_pqcd1_bfom)

store_res_km = []
store_res_pqcd_msbar = []
store_res_pqcd_bfmom = []
store_res_adpy_msbar = []
store_res_adpy_mpole = []
store_res_pqcd1_bfmom = []

for mom in collect(1:9)
    if mom !=4
        push!(store_res_km, uwreal(0.0))
        push!(store_res_pqcd_msbar, uwreal(0.0))
        push!(store_res_pqcd_bfmom, uwreal(0.0))
        push!(store_res_adpy_mpole, uwreal(0.0))
        push!(store_res_adpy_msbar, uwreal(0.0))
        push!(store_res_pqcd1_bfmom,uwreal(0.0))
        continue
    end
    # KM data
    res_aux = compute_integral(itp_km, sqrt(mom), Mz)
    push!(store_res_km, res_aux)
    # pqcd BF-MOM
    res_aux = compute_integral(itp_pqcd_bfmom, sqrt(mom), Mz)
    push!(store_res_pqcd_bfmom, res_aux)
    # pqcd MSBAR
    res_aux = compute_integral(itp_pqcd_msbar, sqrt(mom), Mz)
    push!(store_res_pqcd_msbar, res_aux)
    # AdPy MSBAR
    res_aux = compute_integral(itp_adpy_msbar, sqrt(mom), Mz)
    push!(store_res_adpy_msbar, res_aux)
    # AdPy mpole
    res_aux = compute_integral(itp_adpy_mpole, sqrt(mom), Mz)
    push!(store_res_adpy_mpole, res_aux)
    # pqcd1 BF-MOM
    res_aux = compute_integral(itp_pqcd1_bfmom, sqrt(mom), Mz)
    push!(store_res_pqcd1_bfmom, res_aux)
end

uwerr.(store_res_km)
uwerr.(store_res_pqcd_bfmom)
uwerr.(store_res_pqcd_msbar)
uwerr.(store_res_adpy_msbar)
uwerr.(store_res_adpy_mpole)
uwerr.(store_res_pqcd1_bfmom)

## test with central value integral
quadgk(value.(itp_km),sqrt(2), Mz )


