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



alpha0=1/137.035999180;
Mz = 91.1876

path_adpy_light = "/Users/alessandroconigli/MyDrive/postdoc-mainz/projects/deltalpha/deltalpha-analysis/new_strategy/pert_running/comparison_harvey_code/data/adlerPy/adler_light_4loop.dat"
path_pQCD1 = "/Users/alessandroconigli/Desktop/postdoc-mainz/projects/deltalpha/deltalpha-analysis/new_strategy/pert_running/comparison_harvey_code/data/pQCD1/pQCDAdler.dat33"

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

function get_integral(data, dl) # apply trapezodial rule assuming data to be equally dl distant 
    return sum((data[1:end-1] + data[2:end])*0.5*dl)
end
function get_interp(qval, adler)
    itp = interpolate((qval,),  adler , Gridded(Linear()))
    return itp
end

qval_pqcd, adler_pqcd = read_pQCD_datfile(path_pQCD1, errinp=false)
qval_adpy, adler_adpy = read_adlerpy_datfile(path_adpy_light)

distq = mean(diff(log.(qval_adpy.^2)))

itp_adpy = get_interp(log.(qval_adpy.^2), adler_adpy)
itp_pqcd = get_interp(log.(qval_pqcd.^2), adler_pqcd)
##
store_res1 = []
store_res2 = []
store_res_pqcd = []
store_res_adpy = []

for mom in [5,10,15,20].^2
    ub = findall(x->x - sqrt(mom) >=0, qval_adpy)[1]-1
    lb = findall(x->x - sqrt(mom/4) >=0, qval_adpy)[1]-1
    #pQCD
    xx = collect(range(log(mom/4), log(mom), step=distq))
    mres = get_integral(itp_pqcd(xx), distq) / (12*pi) ./2
    push!(store_res_pqcd, mres)
    #adPy
    # println(qval_adpy[ub], " ", qval_adpy[lb])

    mres1 = get_integral(adler_adpy[lb:ub], distq) / (12*pi) ./2
    mres2 = get_integral(adler_adpy[lb-1:ub], distq) / (12*pi) ./ 2
    mres = mres2 + (mres1 - mres2) / (qval_adpy[lb-1] - qval_adpy[lb]) * (qval_adpy[lb-1] - sqrt(mom/4))
    xx = collect(range(log(mom/4), log(mom), step=distq))
    # mres = get_integral(itp_adpy(xx), distq) / (12*pi) ./ 2
    # push!(store_res1, mres1)
    # push!(store_res2, mres2)
    push!(store_res_adpy, mres)
end