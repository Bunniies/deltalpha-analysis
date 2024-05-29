using Revise
using HVPobs, ADerrors
using PyPlot, LaTeXStrings
using TimerOutputs

#============= SET UP  VARIABLES ===========#

rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
rcParams["text.usetex"] =  true
rcParams["mathtext.fontset"]  = "cm"
rcParams["font.size"] =12
rcParams["axes.labelsize"] =20
rcParams["axes.titlesize"] = 18

plt.rc("text", usetex=true) # set to true once you install latex
#===========================================#

include("tools.jl")

path_data = "/Users/alessandroconigli/Lattice/data/HVP/2ptdata"
path_rw   = "/Users/alessandroconigli/Lattice/data/aux_obs_data/rwf"
path_ms   = "/Users/alessandroconigli/Lattice/data/aux_obs_data/wilson"
path_plot = "/Users/alessandroconigli/MyDrive/postdoc-mainz/projects/deltalpha/plots"
improve = true
renorm  = true

enslist = ["H101"]

@time begin

    println("Reading data...")
    v1v1 = get_corr(path_data, enslist[1], "light", "V1V1", path_rw=path_rw, frw_bcwd=false)#, L=CLS_db[enslist[1]]["L"])
    v2v2 = get_corr(path_data, enslist[1], "light", "V2V2", path_rw=path_rw, frw_bcwd=false)#, L=CLS_db[enslist[1]]["L"])
    v3v3 = get_corr(path_data, enslist[1], "light", "V3V3", path_rw=path_rw, frw_bcwd=false)#, L=CLS_db[enslist[1]]["L"])
    
    v1t10 = get_corr(path_data, enslist[1], "light", "V1T10", path_rw=path_rw, frw_bcwd=false)#, L=CLS_db[enslist[1]]["L"])
    v2t20 = get_corr(path_data, enslist[1], "light", "V2T20", path_rw=path_rw, frw_bcwd=false)#, L=CLS_db[enslist[1]]["L"])
    v3t30 = get_corr(path_data, enslist[1], "light", "V3T30", path_rw=path_rw, frw_bcwd=false)#, L=CLS_db[enslist[1]]["L"])
    
    t0_ens = get_t0(path_ms, enslist[1], path_rw=path_rw, pl=true)

    
    # cv = uwreal([-0.069, 0.032], "test")
    if improve
        beta = CLS_db[enslist[1]]["beta"]
        cv = cv_loc(beta)
        improve_corr_vkvk!(v1v1, v1t10, cv)
        improve_corr_vkvk!(v2v2, v2t20, cv)
        improve_corr_vkvk!(v3v3, v3t30, cv)
    end

    g33 = Corr(0.5 .* (v1v1.obs .+ v2v2.obs .+ v3v3.obs)./3, v1v1.id, v1v1.gamma)
    
        
    if renorm
        beta = CLS_db[enslist[1]]["beta"]
        ml = 0.5 * (1 / CLS_db[enslist[1]]["kappa_l"] - 1 / CLS_kappa_crit[beta] )
        ms = 0.5 * (1 / CLS_db[enslist[1]]["kappa_s"] - 1 / CLS_kappa_crit[beta] )
        trmq = 2*ml + ms
        
        g33.obs .*= (ZV.(beta) .* (1. .+ bv_bar.(beta).*trmq .+ bv.(beta).*ml  )).^2
    end
    frwd_bckwrd_symm!(g33)
end

## Save G33 correlator for comparison

path_corr = "/Users/alessandroconigli/Desktop/g33_H101.txt"

open(path_corr, "w") do f
    write(f, "#t        val                     err \n")
    for k in eachindex(g33.obs)
        write(f, "$(k-1)   $(value(g33.obs[k]))      $(err(g33.obs[k]))  \n")
    end
end
##
wpmm = Dict{String, Vector{Float64}}()
wpmm["H101"] = [-1.0, 4.0, -1.0, -1.0]

uwerr.(v1v1.obs)
uwerr.(v2v2.obs)
uwerr.(v3v3.obs)

g33.obs[2:end] = (g33.obs[2:end] + reverse(g33.obs[2:end])) ./2.

[uwerr(g33.obs[k], wpmm) for k in eachindex(g33.obs)]


##
uwerr.(v1t10.obs)
errorbar(collect(1:length(v1t10.obs)), value.(v1t10.obs), err.(v1t10.obs), fmt="s")
display(gcf())
close("all")
##

#=====================#
# Check Kernel
#=====================#

##
Qgev = [3., 5., 9.] # Q^2
Qgev = [0.1, 1., 3.] # Q^2

#3 * t0sqrt_ph^2 / t0_ens / hc^2 *1e6

#Qlat = Qgev .* value.(t0sqrt_ph.^2 ./ t0_ens ./hc^2 *1e6)

#Qlat_half = Qgev ./2 .* value.(t0sqrt_ph.^2 ./ t0_ens ./hc^2 *1e6)
t = collect(range(0,96,length=1000)) .* value.(t0sqrt_ph ./ sqrt.(t0_ens))
t0_ens = obs[1]["t0"]
for q in Qgev
    qlat = q*1e6/hc^2 #*  value(t0sqrt_ph^2 / t0_ens /hc^2 *1e6)

    k = [kernel_deltalpha(t[k], qlat) for k in eachindex(t)] ./ t.^3#  .* value.(sqrt.(t0_ens) ./ t0sqrt_ph)
    # k = [krnl_dα_qhalf(t[k], qlat) for k in eachindex(t)] ./ t.^3  .* value.(sqrt.(t0_ens) ./ t0sqrt_ph)
    replace!(k, NaN=>0)
    plot(t, k, label=L"$Q^2=-$"*"$(q) GeV" )
    legend()
    ylim(0.0,1.6)
end
xlabel(L"$t/a$")
ylabel(L"$K(t, Q^2)/t^3 \ [\mathrm{fm}^{-1}]$")
display(gcf())
# savefig(joinpath(path_plot, "kernel.pdf"))
close("all")
##
function tmr_integrand(obs::Vector{uwreal}, Q2::Union{Float64,uwreal})#, t0ens::uwreal)
    tvals = length(obs)
    tfm = Float64.(collect(0:tvals-1)) #.* value.(t0sqrt_ph ./ sqrt.(t0ens))
    # k = [kernel_deltalpha(t, Q2) for t in tfm] # .* value.(sqrt.(t0ens) ./ t0sqrt_ph ).^(-2)
    k = [krnl_dα_qhalf(t, Q2) for t in tfm] # .* value.(sqrt.(t0ens) ./ t0sqrt_ph ).^(-2)
    replace!(k, NaN=>0)
    #plot(tfm, k)
    #display(gcf())
    #close()
    integrand = obs .* k
    return integrand
end

##
# integrand = tmr_integrand((g33.obs[2:end] .+ reverse(g33.obs[2:end])) ./2, 27.05)

Qgev = [3, 5, 9] # Q^2
Qlat = Qgev .* (t0sqrt_ph.^2 ./ t0_ens ./hc^2 *1e6)

QGEV = Qgev[3]
QQhalf = Qlat_half[3]

QQ = Qlat[3]
# itegrand_full = tmr_integrand(g33.obs  ,QQ) #, t0_ens)
# integrand_half = tmr_integrand(g33.obs  ,QQhalf) #, t0_ens)
# integrand = (integrand_full .- integrand_half) #./ sqrt.(8 .* t0_ens)  # value.(t0sqrt_ph ./ sqrt.(t0_ens)).^(-2)

figure(figsize=(8,8))
cc = ["gold", "purple", "lightskyblue"]
for (k, q) in enumerate(Qlat)
    integrand= tmr_integrand(g33.obs, q ) 
    sum_int = sum(integrand[1:50]); uwerr(sum_int)
    #integrand = integrand ./ sum_int .*  value.(sqrt.(t0_ens) ./ t0sqrt_ph)
    println(sum_int)
    integrand = integrand ./ sqrt.( 8 .* t0_ens)
    uwerr.(integrand)
    xx = collect(1:length(integrand)) .* value.(t0sqrt_ph ./ sqrt.(t0_ens))
    errorbar(xx, abs.(value.(integrand)), err.(integrand), fmt="o", ms=8, mfc=cc[k], color=cc[k], label=L"$Q^2=-$"*"$(Qgev[k]) " *L"$\mathrm{GeV}^2$")
    fill_between(xx, abs.(value.(integrand)) - err.(integrand), abs.(value.(integrand)) + err.(integrand), alpha=0.5, color=cc[k]  )

end    
ylim(0.0, .0015)
xlim(0., 2.)
# xlabel(L"$t/a$")
xlabel(L"$t\ [\mathrm{fm}]$")
ylabel(L"$G^{33}(t) \cdot K(t, Q^2)/ \sqrt{8t_0}$")

# ylabel(L"$G_{33}(t) \cdot K(t, Q^2)/ \int G_{33}(t) \cdot K(t, Q^2) \ [\mathrm{fm}^{-1}]$")
legend()
#twiny()
#xx_fm =  xx[1:48] #.* value.(t0sqrt_ph ./ sqrt.(t0_ens))
#fm_time = round.((xx_fm[1:4:end] ) .* value.(t0sqrt_ph ./ sqrt.(t0_ens)) , digits=1)
#xticks(xx_fm[1:4:end], fm_time)
#tight_layout()
display(gcf())
tt = string("G33_kernel_", QGEV, "_GeV.pdf")
#savefig(joinpath(path_plot, tt))

close()

## check PT prediction

piq_piqhalf = 0.5 * (3/(12*pi^2)) * log(4) * (1+0.3/pi)
sum(abs.(value.(integrand[1:40])))
sum(integrand[1:40])