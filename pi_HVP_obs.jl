######################################################################################
# This file was created by Alessandro Conigli 
# Here we compute the finite-a HVP required in the analysis of the running of δα
# The observables are then stored in a BDIO file, following the order 
########################################################################################

using Revise
using HVPobs
using PyPlot, LaTeXStrings
using OrderedCollections
using TimerOutputs
using BDIO, ADerrors

#================ SET UP VARIABLES ===============#

rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
rcParams["text.usetex"] =  true
rcParams["mathtext.fontset"]  = "cm"
rcParams["font.size"] =13
rcParams["axes.labelsize"] =22
rcParams["axes.titlesize"] = 18

plt.rc("text", usetex=true) # set to true if a LaTeX installation is present

#======== INCLUDE, PATHS AND CONSTANTS ==========#
# Madrid scale setting
const t0sqrt_ph = uwreal([0.1439, 0.0006], "sqrtt0 [fm]") 

include("./pi_HVP_types.jl")
include("tools.jl")
include("data_management.jl")
include("plot_utils.jl")

path_data = "/Users/alessandroconigli/Lattice/data/HVP/2ptdata"
path_rw   = "/Users/alessandroconigli/Lattice/data/HVP/rwf_deflated"
path_ms   = "/Users/alessandroconigli/Lattice/data/HVP/ms_t0_dat"
path_fvc  = "/Users/alessandroconigli/Lattice/data/HVP/FVC"

path_plot = "/Users/alessandroconigli/MyDrive/postdoc-mainz/projects/deltalpha/plots"
path_bdio = "/Users/alessandroconigli/MyDrive/postdoc-mainz/projects/deltalpha/data"

IMPR      = true
IMPR_SET  = "1" # either "1" or "2"
RENORM    = true
STD_DERIV = true

# enslist = ["H101", "B450", "N202", "N300", "J500"]
enslist = ["H101"]
ensinfo = EnsInfo.(enslist)

Qgev = [3., 5., 9.] # Q^2
Qmgev = 9.0 # Qm^2


#============ OBSERVABLE ALLOCATIONS ============#

g33_ll = Vector{Corr}(undef, length(ensinfo))
g33_lc = Vector{Corr}(undef, length(ensinfo))

g88_ll_conn = Vector{Corr}(undef, length(ensinfo))
g88_lc_conn = Vector{Corr}(undef, length(ensinfo))

g08_ll_conn = Vector{Corr}(undef, length(ensinfo))
g08_lc_conn = Vector{Corr}(undef, length(ensinfo))

obs = Vector(undef, length(ensinfo))


#================ READING DATA ====================#

@info("Reading data...")
@time begin

    for (k, ens) in enumerate(ensinfo)
        
        println("\n - Ensemble: $(ens.id)")

        println("   - G33 ll and lc correlator")
        g33_ll[k], g33_lc[k] = corr33(path_data, ens, sector="light", path_rw=path_rw, impr=IMPR, impr_set=IMPR_SET, cons=true, frw_bcwd=true, std=STD_DERIV)
                
        println("   - G88 connected ll and lc correlator")
        g88_ll_conn[k], g88_lc_conn[k] = corr88_conn(path_data, ens, g33_ll[k], g33_lc=g33_lc[k], path_rw=path_rw, impr=IMPR, impr_set=IMPR_SET, cons=true, frw_bcwd=true, std=STD_DERIV)
        
        println("   - G08 connected ll and lc correlator")
        g08_ll_conn[k], g08_lc_conn[k] = corr08_conn(g33_ll[k], g88_ll_conn[k], g33_lc=g33_lc[k], g88_lc=g88_lc_conn[k])

        if RENORM

            Z3 = get_Z3(ens, impr_set=IMPR_SET)
            renormalize!(g33_ll[k], Z3^2)
            renormalize!(g33_lc[k], Z3)
            Z8 = get_Z8(ens, impr_set=IMPR_SET)
            renormalize!(g88_ll_conn[k], Z8^2)
            renormalize!(g88_lc_conn[k], Z8)
            Z08 = get_Z08(ens, impr_set=IMPR_SET)
            renormalize!(g08_ll_conn[k], Z8*Z08)
            renormalize!(g08_lc_conn[k], Z08)
        end

        println("   - Gradient flow t0")
        obs[k] = OrderedDict()
        obs[k]["t0"] = get_t0(path_ms, ens, path_rw=path_rw, pl=false)
    end
end

#=========== COMPUTE FVC  ==============#
##
for (k, ens) in enumerate(ensinfo)
    
    Qlat = Qgev .* value.(t0sqrt_ph.^2 ./ obs[k]["t0"]) ./hc^2 *1e6
    # Qlat = 1.0 .* 0.1467.^2 ./ obs[k]["t0"] ./hc^2 *1e6

    KRNL = krnl_dα_qhalf_sub

    fvc_raw = get_fvc(path_fvc, ens.id)
    fvc = vcat(sum(fvc_raw, dims=1)...) 
    Thalf = Float64.(collect(1:length(fvc)))

    aux = Vector{uwreal}(undef, 0)

    for (j,q) in enumerate(Qlat) 
        try
            krnl = KRNL.(Thalf, q)
            push!(aux, sum(fvc .* krnl))
        catch
            println("- subtracted kernel")
            qmlat = Qmgev * value.(t0sqrt_ph^2 / obs[k]["t0"]) ./hc^2 *1e6
            krnl = KRNL.(Thalf, q, qmlat)
            push!(aux, sum(fvc .* krnl))
        end
    end
    if ens.kappa_l == ens.kappa_s
        println("Symmetric point ensemble")
        mult = -1.5
    else
        mult = -1
    end
    obs[k]["fvc"] =  mult .* aux
end

## test FVC and plotting at different accumulating wrapping
wpmm = Dict{String, Vector{Float64}}()
wpmm["H101"]     = [5.0, -2.0, -1.0, -1.0]
wpmm["H102r002"] = [5.0, -2.0, -1.0, -1.0]
wpmm["H400"]     = [5.0, -1.5, -1.0, -1.0]
wpmm["N202"]     = [5.0, -2.0, -1.0, -1.0]
wpmm["N200"]     = [5.0, -2.0, -1.0, -1.0]
wpmm["N203"]     = [5.0, -2.0, -1.0, -1.0]
wpmm["N300"]     = [5.0, -1.5, -1.0, -1.0]
wpmm["J303"]     = [5.0, -2.0, -1.0, -1.0]


fvc = Vector{Matrix{uwreal}}(undef, length(ensinfo))
fvc_acc = Vector{Vector{Vector{uwreal}}}(undef, length(ensinfo))
for (k, ens) in enumerate(ensinfo)
    fvc[k] = get_fvc(path_fvc, ens.id)
    fvc_acc[k] = Vector{Vector{uwreal}}(undef, 0)

    KRNL = krnl_dα_qhalf_sub
    qmlat = Qmgev * value(t0sqrt_ph^2 / obs[k]["t0"]) /hc^2 *1e6



        qref = 3.0 * value(t0sqrt_ph^2) / obs[k]["t0"] /hc^2 *1e6
        T = 0
        for (idx,n) in enumerate([1, 2, 6])
            push!(fvc_acc[k],  vcat(sum(fvc[k][1:n,:], dims=1)...))

            T = Float64.(collect(1:length(fvc_acc[k][idx])))
            tmr = -1. *fvc_acc[k][idx] .*  KRNL.(T, qref, qmlat)  
            tmp = sum(tmr); 
            # [uwerr(tmp[ll], wpmm) for ll in eachindex(tmp)]
            uwerr(tmp, wpmm) 
            #println(tmp )
            tmr .=  tmr ./ value.(t0sqrt_ph) .* sqrt.(obs[k]["t0"]) # convert to phys units
            uwerr.(tmr)
            # plot(T, value.(tmr),  label=L"$\vec{n}^2 = $"* "$(n)")
            fill_between(T, value.(tmr)-err.(tmr), value.(tmr)+err.(tmr), alpha=0.8, lw=50, label=L"$F_\pi(-Q^2),\ \vec{n}^2 \leq$"* " $(n)" )
        end

        yy = tmr_integrand(g33_ll[k], qref, qmlat, KRNL, pl=false) ./ value.(t0sqrt_ph) .* sqrt.(obs[k]["t0"])
        # wpm=Dict{String, Vector{Float64}}()
        # wpm["J501"] =[-10.0, -1., 100., -1.] 
        # [uwerr(yy[k], wpm) for k in eachindex(yy)]
        [uwerr(yy[ll], wpmm) for ll in eachindex(yy)]
        fill_between(T, 0, err.(yy), alpha=0.2, color="gray")
        legend()
        xlabel(L"$t/a$")
        ylabel(L"$\delta G^{33}(t)K(t,Q^2) \ [\mathrm{fm}^{-1}]$")
        tight_layout()
        display(gcf())
        t = "$(ens.id)_fvc_Q2_3.pdf"
        # savefig(joinpath(path_plot,t))
        close("all")

end



#=========== COMPUTE TMR ==============#
##
for (k, ens) in enumerate(ensinfo)

    Qlat = Qgev .* value.(t0sqrt_ph.^2) ./ obs[k]["t0"] ./hc^2 *1e6
    qmlat = Qmgev * value(t0sqrt_ph^2) / obs[k]["t0"] /hc^2 *1e6

    # Qlat = 1.0 .* 0.1467.^2 ./ obs[k]["t0"] ./hc^2 *1e6

    pi_33_ll, pi_33_lc, pi_88_ll_conn, pi_88_lc_conn, pi_08_ll_conn, pi_08_lc_conn = [Vector{uwreal}(undef, 0) for _ = 1:6]
    
    KRNL = krnl_dα_qhalf_sub

    for (j, q) in enumerate(Qlat)

        try
            push!(pi_33_ll, sum(tmr_integrand(g33_ll[k], q,  KRNL, pl=true)))
            push!(pi_33_lc, sum(tmr_integrand(g33_lc[k], q,  KRNL, pl=false)))
            
            push!(pi_88_ll_conn, sum(tmr_integrand(g88_ll_conn[k], q,  KRNL, pl=false)))
            push!(pi_88_lc_conn, sum(tmr_integrand(g88_lc_conn[k], q,  KRNL, pl=false)))
        catch
            println("- subtracted kernel")
            push!(pi_33_ll, sum(tmr_integrand(g33_ll[k], q, qmlat, KRNL, pl=true, t0ens=obs[k]["t0"])))
            push!(pi_33_lc, sum(tmr_integrand(g33_lc[k], q, qmlat, KRNL, pl=false)))
            
            push!(pi_88_ll_conn, sum(tmr_integrand(g88_ll_conn[k], q, qmlat, KRNL, pl=false)))
            push!(pi_88_lc_conn, sum(tmr_integrand(g88_lc_conn[k], q, qmlat, KRNL, pl=false)))
        end
        
        #push!(pi_08_ll_conn, sum(tmr_integrand(g08_ll_conn[k], q, KRNL, pl=false)))
        #push!(pi_08_lc_conn, sum(tmr_integrand(g08_lc_conn[k], q, KRNL, pl=false)))
    end

    obs[k]["pi_33_ll"] = pi_33_ll #.+ obs[k]["fvc"]
    obs[k]["pi_33_lc"] = pi_33_lc #.+ obs[k]["fvc"]

    obs[k]["pi_88_ll_conn"] = pi_88_ll_conn #.+ obs[k]["fvc"]
    obs[k]["pi_88_lc_conn"] = pi_88_lc_conn #.+ obs[k]["fvc"]

    obs[k]["pi_08_ll_conn"] = sqrt(3)/2 * (pi_33_ll - pi_88_ll_conn) # pi_08_ll_conn #.+ (all(iszero.(value.(pi_08_ll_conn))) ? 0.0 : obs[k]["fvc"] ) 
    obs[k]["pi_08_lc_conn"] = sqrt(3)/2 * (pi_33_lc - pi_88_lc_conn) #pi_08_lc_conn #.+ (all(iszero.(value.(pi_08_lc_conn))) ? 0.0 : obs[k]["fvc"] )
end

#========== ADD FVC ==========#

for k in eachindex(ensinfo)
    obs[k]["pi_33_ll"] .+= obs[k]["fvc"]
    obs[k]["pi_33_lc"] .+= obs[k]["fvc"]

    obs[k]["pi_88_ll_conn"] .+= obs[k]["fvc"]
    obs[k]["pi_88_lc_conn"] .+= obs[k]["fvc"]

    obs[k]["pi_08_ll_conn"] .+= (all(iszero.(value.(obs[k]["pi_08_ll_conn"]))) ? 0.0 : obs[k]["fvc"] ) 
    obs[k]["pi_08_lc_conn"] .+= (all(iszero.(value.(obs[k]["pi_08_lc_conn"]))) ? 0.0 : obs[k]["fvc"] )
end
##

#======== SAVE TO BDIO =========#
##
for (k, o) in enumerate(obs)

    ens = ensinfo[k].id
    p = joinpath(path_bdio, ens,  string(ens, "_PI_hvp_set$(IMPR_SET).bdio"))
    fb = BDIO_open(p, "w", ens)
    uinfo = 0 

    for (key, value) in o
        if isa(value, uwreal)
            write_uwreal(value, fb, uinfo)
        elseif isa(value, Vector{uwreal})
            [write_uwreal(v, fb, uinfo) for v in value]
        end
        uinfo +=1
    end

    BDIO_close!(fb)
end
##

test_pp = "/Users/alessandroconigli/MyDrive/postdoc-mainz/projects/deltalpha/data/B450/B450_PI_hvp.bdio"
cc = read_BDIO(test_pp, "dalpha", "pi_33_ll")




## OLDER CHECKS
#======== check continuum limit =======#
# global cl for local and  conserved
Res = []
xdata = [1 / (8 * obs[k]["t0"]) for k in eachindex(obs)]
idx=2:5
uwerr.(xdata)
@. cl_loc(x,p) = p[1] + p[2]*x + p[3]*x^(3/2.) #+ p[4]*x^2
@. cl_cons(x,p) = p[1] + p[4]*x + p[5]*x^(3/2.) #+ p[4]*x^2
col = ["plum", "gold", "lightskyblue"]
for (k,q) in enumerate(Qgev)
    # local local
    ydata_loc = [obs[j]["pi_33_ll"][k] for j in eachindex(ensinfo)]; uwerr.(ydata_loc)
    ydata_cons = [obs[j]["pi_33_lc"][k] for j in eachindex(ensinfo)]; uwerr.(ydata_cons)
    fit = fit_routine([cl_loc, cl_cons], [value.(xdata)[idx], value.(xdata)[idx]], [ydata_loc[idx], ydata_cons[idx]], 5)
    cl_val = fit.param[1]; uwerr(cl_val)
    push!(Res, cl_val)
    xarr = Float64.(range(0.0, maximum(value.(xdata[idx])), length=100))
    yarr = cl_loc(xarr, fit.param); uwerr.(yarr)
    
    errorbar(value.(xdata), xerr=err.(xdata), value.(ydata_loc), yerr=err.(ydata_loc), fmt="s", capsize=2, label=L"$Q^2= $"*"$(Qgev[k]) "*L"$\mathrm{GeV}^2$", color=col[k], mfc="none" )
    errorbar(0.0, value(cl_val), err(cl_val), fmt="s", mfc="none", capsize=2, color=col[k])

    fill_between(xarr, value.(yarr)-err.(yarr), value.(yarr)+err.(yarr), alpha=0.2, color=col[k] )

    # local conserved
    cl_val = fit.param[1]; uwerr(cl_val)
    xarr = Float64.(range(0.0, maximum(value.(xdata[idx])), length=100))
    yarr = cl_cons(xarr, fit.param); uwerr.(yarr)
    
    errorbar(value.(xdata), xerr=err.(xdata), value.(ydata_cons), yerr=err.(ydata_cons), fmt="d", capsize=2, color=col[k], mfc="none" )

    fill_between(xarr, value.(yarr)-err.(yarr), value.(yarr)+err.(yarr), alpha=0.2, color=col[k] )

end
axvline(ls="dashed", color="black", lw=0.2, alpha=0.7) 
xlabel(L"$a^2/8t_0$")
ylabel(L"$\bar{\Pi}^{33}(-Q^2)$")
legend(loc="upper left")
ylim(0.0135,0.02)
# ylim(0.0175,0.0245)
xlim(-0.001, 0.05)
tight_layout()
#savefig(joinpath(path_plot,"pi33_cont_lim.pdf"))
display(gcf())
close()
##

#======== check continuum limit =======#
# independent cl for local and  conserved
xdata = [1 / (8 * obs[k]["t0"]) for k in eachindex(obs)]
uwerr.(xdata)
@. cl(x,p) = p[1] + p[2]*x + p[3]*x^(3/2.) #+ p[4]*x^2
col = ["plum", "gold", "lightskyblue"]
idx = findall(x->x>3.4, getfield.(ensinfo, :beta))
for (k,q) in enumerate(Qgev)
    # local local
    ydata = [obs[j]["pi_33_ll"][k] for j in eachindex(ensinfo)]
    uwerr.(ydata)
    fit = fit_routine(cl, value.(xdata)[idx], ydata[idx], 3)
    cl_val = fit.param[1]; uwerr(cl_val)
    xarr = Float64.(range(0.0, maximum(value.(xdata)[idx]), length=100))
    yarr = cl(xarr, fit.param); uwerr.(yarr)
    
    errorbar(value.(xdata), xerr=err.(xdata), value.(ydata), yerr=err.(ydata), fmt="s", capsize=2, label=L"$Q^2= $"*"$(Qgev[k]) "*L"$\mathrm{GeV}^2$", color=col[k], mfc="none" )
    #errorbar(0.0, value(cl_val), err(cl_val), fmt="s", mfc="none", capsize=2, color=col[k])

    fill_between(xarr, value.(yarr)-err.(yarr), value.(yarr)+err.(yarr), alpha=0.2, color=col[k] )

    # local conserved
    ydata = [obs[j]["pi_33_lc"][k] for j in eachindex(ensinfo)]
    uwerr.(ydata)
    fit = fit_routine(cl, value.(xdata)[idx], ydata[idx], 3)
    cl_val = fit.param[1]; uwerr(cl_val)
    xarr = Float64.(range(0.0, maximum(value.(xdata)[idx]), length=100))
    yarr = cl(xarr, fit.param); uwerr.(yarr)
    
    errorbar(value.(xdata), xerr=err.(xdata), value.(ydata), yerr=err.(ydata), fmt="d", capsize=2, color=col[k], mfc="none" )
    #errorbar(-0.001, value(cl_val), err(cl_val), fmt="d", mfc="none", capsize=2, color=col[k])

    fill_between(xarr, value.(yarr)-err.(yarr), value.(yarr)+err.(yarr), alpha=0.2, color=col[k] )

end
axvline(ls="dashed", color="black", lw=0.2, alpha=0.7) 
xlabel(L"$a^2/8t_0$")
ylabel(L"$\bar{\Pi}^{33}(-Q^2)$")
# ylabel(L"$\bar{\Pi}^{33,\mathrm{sub}}(-Q^2)$")
legend(loc="lower right")
ylim(0.013,0.023)
# ylim(0.0115,0.02)
tight_layout()
xlim(-0.002, 0.05)
display(gcf())
savefig(joinpath(path_plot, "su3_ens_cl.pdf"))
close()

##
# CHECK KERNEL
Qgev = [3., 5., 9.] # Q^2
Qm = 8.
# Qgev = [0.1, 1., 3.] # Q^2

#3 * t0sqrt_ph^2 / t0_ens / hc^2 *1e6

#Qlat = Qgev .* value.(t0sqrt_ph.^2 ./ t0_ens ./hc^2 *1e6)

#Qlat_half = Qgev ./2 .* value.(t0sqrt_ph.^2 ./ t0_ens ./hc^2 *1e6)
t = collect(range(0,96,length=1000)) .* value.(t0sqrt_ph ./ sqrt.(obs[1]["t0"]))
t0_ens = obs[1]["t0"]
kk=1
for q in Qgev
    plid = 310 + kk
    subplot(plid)
    qlat = q*1e6/hc^2 #*  value(t0sqrt_ph^2 / t0_ens /hc^2 *1e6)
    qmlat = Qm*1e6/hc^2
    k1 = [krnl_dα_qhalf_sub(t[k], qlat, qmlat) for k in eachindex(t)] ./ t.^3#  .* value.(sqrt.(t0_ens) ./ t0sqrt_ph)
    k = [krnl_dα_qhalf(t[k], qlat) for k in eachindex(t)] ./ t.^3#  .* value.(sqrt.(t0_ens) ./ t0sqrt_ph)
    # k = [krnl_dα_qhalf(t[k], qlat) for k in eachindex(t)] ./ t.^3  .* value.(sqrt.(t0_ens) ./ t0sqrt_ph)
    replace!(k, NaN=>0)
    replace!(k1, NaN=>0)
    plot(t, k)#, label=L"$Q^2=-$"*"$(q) "* L"$\mathrm{GeV}^2$" )
    plot(t, k1, ls="dashed")#, label=L"$Q^2=-$"*"$(q) GeV", ls="dashed" )
    ax = gca()
    text(0.7,0.6,L"$Q^2=\ $"*"$(q) "* L"$\mathrm{GeV}^2$", transform=ax.transAxes )
    #legend()
    #ylim(0.0,1.6)
    if kk != 3
        setp(ax.get_xticklabels(),visible=false) # Disable x tick labels
    end
    if kk ==2
        ylabel(L"$K(t, Q^2)/t^3 \ [\mathrm{fm}^{-1}]$")
    end
    kk+=1
    xlim(0,4)
end
xlabel(L"$t\ [\mathrm{fm}]$")
display(gcf())
savefig(joinpath(path_plot, "kernel_comparison.pdf"))
close("all")








##

3 * t0sqrt_ph^2 / t0_ens / hc^2 *1e6

Qlat = Qgev .* value.(t0sqrt_ph.^2 ./ t0_ens ./hc^2 *1e6)

Qlat_half = Qgev ./2 .* value.(t0sqrt_ph.^2 ./ t0_ens ./hc^2 *1e6)
t = collect(range(0,96,length=1000)) #.* value.(t0sqrt_ph ./ sqrt.(t0_ens))

for q in Qgev
    qlat = q *  value(t0sqrt_ph^2 / t0_ens /hc^2 *1e6)

    # k = [kernel_deltalpha(t[k], qlat) for k in eachindex(t)] ./ t.^3  .* value.(sqrt.(t0_ens) ./ t0sqrt_ph)
    k = [krnl_dα_qhalf(t[k], qlat) for k in eachindex(t)] ./ t.^3  .* value.(sqrt.(t0_ens) ./ t0sqrt_ph)
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
savefig(joinpath(path_plot, tt))

close()

## check PT prediction

piq_piqhalf = 0.5 * (3/(12*pi^2)) * log(4) * (1+0.3/pi)
sum(abs.(value.(integrand[1:40])))
sum(integrand[1:40])