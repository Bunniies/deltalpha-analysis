######################################################################################
# This file was created by Alessandro Conigli 
# Here we load the finite-a HVP  stored in BDIO and perform
# the chiral-continuum extrapolations
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


include("./pi_HVP_types.jl")
include("./pi_HVP_func_comb.jl")
include("./tools.jl")
include("./data_management.jl")
include("./plot_utils.jl")
include("./IO_BDIO.jl")

path_bdio = "/Users/alessandroconigli/MyDrive/postdoc-mainz/projects/deltalpha/data"
path_plot = "/Users/alessandroconigli/MyDrive/postdoc-mainz/projects/deltalpha/plots"

#======= PHYSICAL CONSTANTS ====================#
const t0sqrt_ph = uwreal([0.1439, 0.0006], "sqrtt0 [fm]") 
const MPI_ph = uwreal([134.9768, 0.0005], "mpi phys")
const MK_ph  = uwreal([495.011, 0.01], "mk phys")

const phi2_ph = (sqrt(8)*t0sqrt_ph * MPI_ph / hc)^2
const phi4_ph = (sqrt(8)*t0sqrt_ph)^2 * ((MK_ph/hc)^2 + 0.5*(MPI_ph/hc)^2)

const TREELEVEL = true
const Qgev = [3., 5., 9.]


#============== READ BDIO FILES =================#

dir_path = filter(isdir, readdir(path_bdio, join=true))
hvp_path_s1 = vcat(filter(!isempty, [filter(x-> occursin("hvp_set1.bdio", x)  , readdir(dir_path[k], join=true)) for k in eachindex(dir_path)])...)
hvp_path_s2 = vcat(filter(!isempty, [filter(x-> occursin("hvp_set2.bdio", x)  , readdir(dir_path[k], join=true)) for k in eachindex(dir_path)])...)

spectrum_path = vcat(filter(!isempty, [filter(x-> occursin("spectrum.bdio", x)  , readdir(dir_path[k], join=true)) for k in eachindex(dir_path)])...)

treelevel_path_s1 = vcat(filter(!isempty, [filter(x-> occursin("tree_level_set1.bdio", x)  , readdir(dir_path[k], join=true)) for k in eachindex(dir_path)])...)
treelevel_path_s2 = vcat(filter(!isempty, [filter(x-> occursin("tree_level_set2.bdio", x)  , readdir(dir_path[k], join=true)) for k in eachindex(dir_path)])...)


enslist = [basename(hvp_path_s1[k])[1:4] for k in eachindex(hvp_path_s1)]
ensinfo = EnsInfo.(enslist)
NENS = length(ensinfo)

# set 1 improvement coefficients
pi_33_ll_s1 = Vector{Vector{uwreal}}(undef, NENS)
pi_33_lc_s1 = Vector{Vector{uwreal}}(undef, NENS)

# set 2 improvement coefficients
pi_33_ll_s2 = Vector{Vector{uwreal}}(undef, NENS)
pi_33_lc_s2 = Vector{Vector{uwreal}}(undef, NENS)

phi2 = Vector{uwreal}(undef, NENS)
phi4 = Vector{uwreal}(undef, NENS)
a28t0 = Vector{uwreal}(undef, NENS)

# reading tree-level prediction in the continuum
pt_pred = []
p = joinpath(path_bdio, string("tree_level_cl_phys_val.bdio"))
fb = BDIO_open(p, "r")
BDIO_seek!(fb)
if BDIO_get_uinfo(fb) == 0
    push!(pt_pred, read_uwreal(fb))
end
while BDIO_seek!(fb, 2)
    if BDIO_get_uinfo(fb) == 0 
        push!(pt_pred, read_uwreal(fb))
    end
end
BDIO_close!(fb)

for (k, ens) in enumerate(ensinfo)
    println("- Ensemble ", ens.id)
    if TREELEVEL
        #pt_pred = @. (1-(sqrt(Qgev)/(2*3))^2 ) * 1 / (4pi^2) * log(2)

        hvp_3l_ll_s1 = value.(pt_pred) .- read_BDIO(treelevel_path_s1[k], "3level", "3l_33_ll")
        hvp_3l_lc_s1 = value.(pt_pred) .- read_BDIO(treelevel_path_s1[k], "3level", "3l_33_lc")
        hvp_3l_ll_s2 = value.(pt_pred) .- read_BDIO(treelevel_path_s2[k], "3level", "3l_33_ll")
        hvp_3l_lc_s2 = value.(pt_pred) .- read_BDIO(treelevel_path_s2[k], "3level", "3l_33_lc")
    else
        hvp_3l_ll_s1, hvp_3l_lc_s1, hvp_3l_ll_s2, hvp_3l_lc_s2 = [fill(1.0, length(Qgev)) for _ in 1:4]
    end

    # set 1 improvement coefficients
    pi_33_ll_s1[k] = read_BDIO(hvp_path_s1[k], "dalpha", "pi_33_ll" ) .+ hvp_3l_ll_s1
    pi_33_lc_s1[k] = read_BDIO(hvp_path_s1[k], "dalpha", "pi_33_lc" ) .+ hvp_3l_lc_s1

    # set 2 improvement coefficients
    pi_33_ll_s2[k] = read_BDIO(hvp_path_s2[k], "dalpha", "pi_33_ll" ) .+ hvp_3l_ll_s2
    pi_33_lc_s2[k] = read_BDIO(hvp_path_s2[k], "dalpha", "pi_33_lc" ) .+ hvp_3l_lc_s2

    # spectrum
    m_pi  = read_BDIO(spectrum_path[k], "spectrum", "mpi")[1]
    m_k   = read_BDIO(spectrum_path[k], "spectrum", "mk")[1]
    t0ens = read_BDIO(spectrum_path[k], "spectrum", "t0")[1]

    phi2[k]  = 8 * t0ens * m_pi^2
    phi4[k]  = 8 * t0ens * (m_k^2 + 0.5*m_pi^2 )
    a28t0[k] =  1 / (8*t0ens)
end

##############################
## CREATE FIT CATEGORIES
##############################
NMOM = length(pi_33_ll_s1[1])
# set 1 pi33
fitcat_pi33_ll_s1 = [Vector{FitCat}(undef,0) for i=1:NMOM]
fitcat_pi33_lc_s1 = [Vector{FitCat}(undef,0) for i=1:NMOM]
# set 2 pi33
fitcat_pi33_ll_s2 = [Vector{FitCat}(undef,0) for i=1:NMOM]
fitcat_pi33_lc_s2 = [Vector{FitCat}(undef,0) for i=1:NMOM]

## no cuts in data
for s in 1:2
    xdata = [a28t0 phi2 phi4]
    if s == 1
        for q in 1:length(pi_33_ll_s1[1])
            str = "all_data_set$(s)_q$(q)"
            # pi 33
            push!(fitcat_pi33_ll_s1[q], FitCat(xdata, getindex.(pi_33_ll_s1, q), str))
            push!(fitcat_pi33_lc_s1[q], FitCat(xdata, getindex.(pi_33_lc_s1, q), str))
        end
    elseif s == 2
        for q in 1:length(pi_33_ll_s2[1])
            str = "all_data_set$(s)_q$(q)"
            # pi 33
            push!(fitcat_pi33_ll_s2[q], FitCat(xdata, getindex.(pi_33_ll_s2, q), str))
            push!(fitcat_pi33_lc_s2[q], FitCat(xdata, getindex.(pi_33_lc_s2, q), str))
        end

    end
end

## cuts in phi2<0.6 
i_cutphi2 = findall(x->x<0.6, value.(phi2))
for s in 1:2
    xdata = [a28t0[i_cutphi2] phi2[i_cutphi2] phi4[i_cutphi2]]
    if s == 1
        for q in 1:length(pi_33_ll_s1[1])
            str = "all_data_set$(s)_q$(q)"
            # pi 33
            push!(fitcat_pi33_ll_s1[q], FitCat(xdata, getindex.(pi_33_ll_s1, q)[i_cutphi2], str))
            push!(fitcat_pi33_lc_s1[q], FitCat(xdata, getindex.(pi_33_lc_s1, q)[i_cutphi2], str))
        end
    elseif s == 2
        for q in 1:length(pi_33_ll_s2[1])
            str = "all_data_set$(s)_q$(q)"
            # pi 33
            push!(fitcat_pi33_ll_s2[q], FitCat(xdata, getindex.(pi_33_ll_s2, q)[i_cutphi2], str))
            push!(fitcat_pi33_lc_s2[q], FitCat(xdata, getindex.(pi_33_lc_s2, q)[i_cutphi2], str))
        end

    end
end

## cuts in phi2<0.6 and beta<3.4
i_cutphi2 = findall(x->x<0.6, value.(phi2))
i_cutbeta = findall(x->x>3.4, getfield.(ensinfo, :beta) )
i_cutall = intersect(i_cutphi2, i_cutbeta)
for s in 1:2
    xdata = [a28t0[i_cutall] phi2[i_cutall] phi4[i_cutall]]
    if s == 1
        for q in 1:length(pi_33_ll_s1[1])
            str = "phi2_and_beta_cuts_set$(s)_q$(q)"
            # pi 33
            push!(fitcat_pi33_ll_s1[q], FitCat(xdata, getindex.(pi_33_ll_s1, q)[i_cutall], str))
            push!(fitcat_pi33_lc_s1[q], FitCat(xdata, getindex.(pi_33_lc_s1, q)[i_cutall], str))
        end
    elseif s == 2
        for q in 1:length(pi_33_ll_s2[1])
            str = "phi2_and_beta_cuts_set$(s)_q$(q)"
            # pi 33
            push!(fitcat_pi33_ll_s2[q], FitCat(xdata, getindex.(pi_33_ll_s2, q)[i_cutall], str))
            push!(fitcat_pi33_lc_s2[q], FitCat(xdata, getindex.(pi_33_lc_s2, q)[i_cutall], str))
        end

    end
end


#######################################
## CHIRAL-CONTINUUM EXTRAPOLATION
#######################################

# pi 33 set 1
for q in 1:NMOM

    for k_cat in eachindex(fitcat_pi33_ll_s1[q])
        xdata = fitcat_pi33_ll_s1[q][k_cat].xdata
        ydata_ll_s1 = fitcat_pi33_ll_s1[q][k_cat].ydata
        ydata_ll_s2 = fitcat_pi33_ll_s2[q][k_cat].ydata
        ydata_lc_s1 = fitcat_pi33_lc_s1[q][k_cat].ydata
        ydata_lc_s2 = fitcat_pi33_lc_s2[q][k_cat].ydata

        for (k_mod, model) in enumerate(f_tot_isov)
            println(k_mod)
            fit_ll_s1 = fit_routine(model, value.(xdata), ydata_ll_s1, n_par_tot_isov[k_mod])
            fit_ll_s2 = fit_routine(model, value.(xdata), ydata_ll_s2, n_par_tot_isov[k_mod])
            fit_lc_s1 = fit_routine(model, value.(xdata), ydata_lc_s1, n_par_tot_isov[k_mod])
            fit_lc_s2 = fit_routine(model, value.(xdata), ydata_lc_s2, n_par_tot_isov[k_mod])
            push!(fitcat_pi33_ll_s1[q][k_cat].fit, fit_ll_s1)
            push!(fitcat_pi33_ll_s2[q][k_cat].fit, fit_ll_s2)
            push!(fitcat_pi33_lc_s1[q][k_cat].fit, fit_lc_s1)
            push!(fitcat_pi33_lc_s2[q][k_cat].fit, fit_lc_s2)

        end
    end
end



## plot attempt - CL
for q in 1:NMOM

    for k_cat in eachindex(fitcat_pi33_ll_s1[q])
        xdata = fitcat_pi33_ll_s1[q][k_cat].xdata
        ydata_ll_s1 = fitcat_pi33_ll_s1[q][k_cat].ydata
        ydata_ll_s2 = fitcat_pi33_ll_s2[q][k_cat].ydata
        ydata_lc_s1 = fitcat_pi33_lc_s1[q][k_cat].ydata
        ydata_lc_s2 = fitcat_pi33_lc_s2[q][k_cat].ydata

    
        ww_ll_s1 = get_w_from_fitcat(fitcat_pi33_ll_s1[q], norm=true)
        ww_ll_s2 = get_w_from_fitcat(fitcat_pi33_ll_s2[q], norm=true)
        ww_lc_s1 = get_w_from_fitcat(fitcat_pi33_lc_s1[q], norm=true)
        ww_lc_s2 = get_w_from_fitcat(fitcat_pi33_lc_s2[q], norm=true)

        for (k_mod, model) in enumerate(f_tot_isov)
            xarr = [Float64.(range(1e-8, 0.05, length=100)) fill(phi2_ph, 100) fill(phi4_ph, 100)]
            
            fit_param_ll_s1 = fitcat_pi33_ll_s1[q][k_cat].fit[k_mod].param
            fit_param_ll_s2 = fitcat_pi33_ll_s2[q][k_cat].fit[k_mod].param
            fit_param_lc_s1 = fitcat_pi33_lc_s1[q][k_cat].fit[k_mod].param
            fit_param_lc_s2 = fitcat_pi33_lc_s2[q][k_cat].fit[k_mod].param

            yy_ll_s1  = model(xarr, fit_param_ll_s1)
            yy_ll_s2  = model(xarr, fit_param_ll_s2)
            yy_lc_s1  = model(xarr, fit_param_lc_s1)
            yy_lc_s2  = model(xarr, fit_param_lc_s2)

            plot(xarr[:,1], value.(yy_ll_s1),alpha=ww_ll_s1[k_mod], color="forestgreen")
            plot(xarr[:,1], value.(yy_ll_s2),alpha=ww_ll_s2[k_mod], color="royalblue")
            plot(xarr[:,1], value.(yy_lc_s1),alpha=ww_lc_s1[k_mod], color="purple")
            plot(xarr[:,1], value.(yy_lc_s2),alpha=ww_lc_s2[k_mod], color="gold")

        end
    end
    L2D = PyPlot.matplotlib.lines.Line2D
    custom_lines = [L2D([0], [0], color="forestgreen", lw=2),
                L2D([0], [0], color="royalblue", lw=2),
                L2D([0], [0], color="purple", lw=2),
                L2D([0], [0], color="gold", lw=2)]
    legend(custom_lines,  [L"$\mathrm{LL,\ set\ 1}$", L"$\mathrm{LL,\ set\ 2}$", L"$\mathrm{LC,\ set\ 1}$", L"$\mathrm{LC,\ set\ 2}$"])
    if q == 1
        ylim(0.0158, 0.0212)
    elseif q ==2 
        ylim(0.014, 0.021)
    elseif q==3
        ylim(0.011, 0.019)
    end
    xlabel(L"$a^2/(8t_0)$")
    ylabel(L"$\bar{\Pi}^{33,\mathrm{sub}}(-Q^2)$")
    tight_layout()
    display(gcf())
    fname = string("cont_lim_3limp_q", q, ".pdf")
    savefig(joinpath(path_plot, fname))
    close("all")
end

## plot attemp - chiral limit
using Statistics

#i_vec_cuts = [collect(1:length(ensinfo)), i_cutphi2]
for q in 1:NMOM
    tt = ["set", "2", "LL"]
    for k_cat in eachindex(fitcat_pi33_lc_s2[q])

        xdata = fitcat_pi33_ll_s2[q][k_cat].xdata; uwerr.(xdata)
        ydata = fitcat_pi33_ll_s2[q][k_cat].ydata
        
        wtot = get_w_from_fitcat(fitcat_pi33_ll_s2[q], norm=false)
        w, widx = findmax(wtot)
        model = f_tot_isov[widx]
        
        fit_par = fitcat_pi33_ll_s2[q][k_cat].fit[widx].param
        xdata_proj = [xdata[:,1] xdata[:,2] fill(phi4_ph, length(xdata[:,1]))]
        ydata_proj = ydata - model(xdata, fit_par) + model(xdata_proj, fit_par)

        uwerr.(ydata_proj)

        # phys res
        ph_res = model([0.0 phi2_ph phi4_ph], fit_par)[1]; uwerr(ph_res)
        println(ph_res)
        errorbar(value(phi2_ph), value(ph_res), err(ph_res), fmt="o", capsize=2, color="red")
        axvline(value(phi2_ph), ls="dashed", color="black", lw=0.2, alpha=0.7) 

        # plot data points
        a_sp = [0.087, 0.077, 0.065, 0.05, 0.039]
        color = ["tomato", "forestgreen", "violet", "orange", "navy"]

        fmttot = ["d", "s", "^", "h", "8"]
        for (k, b) in enumerate(sort(unique(getfield.(ensinfo, :beta))))
            # n_ = findall(x->x.beta == b, ensinfo[i_vec_cuts[k_cat]])
            n_ = findall(x->x.beta == b, ensinfo)
            a2_aux  = mean(value.(xdata[:,1][n_]))
            errorbar(value.(xdata[n_,2]), value.(ydata_proj[n_]), xerr=err.(xdata[n_,2]), yerr=err.(ydata_proj[n_]), fmt=fmttot[k], label=string(L"$\beta = \ $", b), color=color[k], capsize=2 )
            # dashed lines
            xxx = [fill(a2_aux, 100) Float64.(range(0.0, 0.8, length=100)) fill(phi4_ph, 100)]
            yy_ll = model(xxx, fit_par)
            replace!(value.(yy_ll), Inf=>0.0)
            plot(xxx[:,2], value.(yy_ll), ls="--", color=color[k], lw=0.5)
        end

        # cont lim band
        xarr = [fill(1e-8, 100) Float64.(range(0.01, 0.8, length=100)) fill(phi4_ph,100)]
        yarr = model(xarr, fit_par); uwerr.(yarr)
        fill_between(xarr[:,2], value.(yarr) .- err.(yarr), value.(yarr) .+ err.(yarr), alpha=0.2, color="gray")



    end
    legend(ncol=2)
    xlabel(L"$\phi_2$")
    ylabel(L"$\bar{\Pi}^{33,\mathrm{sub}}(-Q^2)$")
    title(join(tt, " "))
    display(gcf())
    tight_layout()
    fname = string("chiral_lim_pi33_", join(tt,"_"), "q_$(q).pdf") 
    savefig(joinpath(path_plot, fname ))
    close("all")
end


###################################
## RESULTS
##################################


for q in 1:NMOM
    @info "Momentum no. $(q)"
    fitcat_pi33_tot = vcat(vcat(fitcat_pi33_ll_s1[q],
                #fitcat_pi33_ll_s2[q],
                fitcat_pi33_lc_s1[q],
                fitcat_pi33_lc_s2[q])...)

    ww_tot = get_w_from_fitcat(fitcat_pi33_tot)

    w, widx  =  findmax(ww_tot)
  
    model_idx = mod(widx, length(f_tot_isov))
    println("   wmax: ", w)
    model = f_tot_isov[model_idx]
    
    
    cat_idx = Int((widx - model_idx ) / length(f_tot_isov))+1
    if cat_idx < 0
        cat_idx +=1
    end
    println("   best χ2/χ2exp: ", fitcat_pi33_tot[cat_idx].fit[model_idx].chi2 / fitcat_pi33_tot[cat_idx].fit[model_idx].chi2exp)
    println("   widx: ", widx, " model_idx: ", model_idx, " catidx: ", cat_idx)
    ## Best Res
    best_mod = f_tot_isov[model_idx]
    xdata = fitcat_pi33_tot[cat_idx].xdata
    param = fitcat_pi33_tot[cat_idx].fit[model_idx].param

    ph_res_best = best_mod([0.0 phi2_ph phi4_ph], param)[1]; uwerr(ph_res_best)
    println("   best res: ", ph_res_best )
    ## histogram

    all_res = Vector{uwreal}()
    for (k, cat) in enumerate(fitcat_pi33_tot)
        for (j, mod) in enumerate(f_tot_isov)
            push!(all_res, mod([0.0 phi2_ph phi4_ph], cat.fit[j].param)[1])
        end
    end

    final_res, syst = model_average(all_res, ww_tot); uwerr(final_res)
    println("   Model ave:  ", final_res)
    println("   systematic: ", syst)


    #hist(value.(all_res), bins=800, histtype="stepfilled", alpha=0.5, ec="k", color="navy", weights=ww_tot)
    #display(gcf())
    #close()

end





## OLD TESTS ∼ 11/2023
## functional forms

chiral_dep_pi33_ll(x, p) = p[1] .+ p[2] .* (x[:,2] .- p[3]) .+ p[4] .* log.(x[:,2] ./ p[3])
chiral_dep_pi33_lc(x, p) = p[1] .+ p[2] .* (x[:,2] .- p[3]) .+ p[4] .* log.(x[:,2] ./ p[3])

cl_dep_pi33_ll(x, p) = p[5] .* x[:,1] .+ p[6] .* x[:,1].^(3/2)
cl_dep_pi33_lc(x, p) = p[7] .* x[:,1] .+ p[8] .* x[:,1].^(3/2)

fit_pi33_ll = (x,p) -> chiral_dep_pi33_ll(x,p) + cl_dep_pi33_ll(x,p)
fit_pi33_lc = (x,p) -> chiral_dep_pi33_lc(x,p) + cl_dep_pi33_lc(x,p)

## test fit
fix_mom_pi33_ll = getindex.(pi_33_ll, 3)
fix_mom_pi33_lc = getindex.(pi_33_lc, 3)

fit = fit_routine([fit_pi33_ll, fit_pi33_lc], [value.(xfit), value.(xfit)], [fix_mom_pi33_ll, fix_mom_pi33_lc], 8 )


ph_res_ll = fit_pi33_ll([0.0 phi2_ph], fit.param)
ph_res_lc = fit_pi33_lc([0.0 phi2_ph], fit.param)

## plot cl limit
xin = xfit
xproj = [xfit[:,1] fill(phi2_ph, length(xfit[:,1]))]

yproj_ll = fix_mom_pi33_ll - fit_pi33_ll(xin, fit.param) + fit_pi33_ll(xproj , fit.param)
yproj_lc = fix_mom_pi33_lc - fit_pi33_ll(xin, fit.param) + fit_pi33_ll(xproj , fit.param)

uwerr.(yproj_ll)
uwerr.(yproj_lc)

errorbar(value.(xfit[:,1]), value.(yproj_ll), err.(yproj_ll), fmt="d", capsize=2, mfc="none", color="forestgreen", label=L"$\overline{\Pi}^{33, LL}(-Q^2)$")
errorbar(value.(xfit[:,1]), value.(yproj_lc), err.(yproj_lc), fmt="^", capsize=2, mfc="none", color="royalblue", label=L"$\overline{\Pi}^{33, LC}(-Q^2)$")

# CL result
ph_res = fit_pi33_ll([0.0 phi2_ph], fit.param)[1]; uwerr(ph_res)
errorbar(0.0, value(ph_res), err(ph_res), fmt="o", capsize=2, color="red")

# CL Bands
xarr = [Float64.(range(1e-8, 0.05, length=100)) fill(phi2_ph, 100)]
yarr_ll = fit_pi33_ll(xarr, fit.param); uwerr.(yarr_ll)
yarr_lc = fit_pi33_lc(xarr, fit.param); uwerr.(yarr_lc)

fill_between(xarr[:,1], value.(yarr_ll) .- err.(yarr_ll), value.(yarr_ll) .+ err.(yarr_ll), alpha=0.2, color="forestgreen")
fill_between(xarr[:,1], value.(yarr_lc) .- err.(yarr_lc), value.(yarr_lc) .+ err.(yarr_lc), alpha=0.2, color="royalblue")

axvline(ls="dashed", color="black", lw=0.2, alpha=0.7) 
xlabel(L"$a^2/8t_0$")
ylabel(L"$\bar{\Pi}^{33}(-Q^2)$")
ylim(0.018, 0.0235)
legend()
tight_layout()
display(gcf())
tt = "cl_proj_Q2_3.pdf"
savefig(joinpath(path_plot, tt))
close("all")

## PLOT PHI2 limit
xin = xfit
xproj = [xfit[:,1] fill(phi2_ph, length(xfit[:,1]))]

yproj_ll = fix_mom_pi33_ll #- fit_pi33_ll(xin, fit.param) + fit_pi33_ll(xproj , fit.param)
yproj_lc = fix_mom_pi33_lc #- fit_pi33_ll(xin, fit.param) + fit_pi33_ll(xproj , fit.param)

uwerr.(yproj_ll)
uwerr.(yproj_lc)

# errorbar(value.(xfit[:,2]), value.(yproj_ll), err.(yproj_ll), fmt="d", capsize=2, mfc="none", color="forestgreen", label=L"$\overline{\Pi}^{33, LL}(-Q^2)$")
# errorbar(value.(xfit[:,2]), value.(yproj_lc), err.(yproj_lc), fmt="^", capsize=2, mfc="none", color="royalblue", label=L"$\overline{\Pi}^{33, LC}(-Q^2)$")

# CL result
ph_res = fit_pi33_ll([0.0 phi2_ph], fit.param)[1]; uwerr(ph_res)
errorbar(value(phi2_ph), value(ph_res), err(ph_res), fmt="o", capsize=2, color="red")

# CL Bands
xarr = [fill(1e-8, 100) Float64.(range(0.01, 0.8, length=100))]
yarr_ll = fit_pi33_ll(xarr, fit.param); uwerr.(yarr_ll)
yarr_lc = fit_pi33_lc(xarr, fit.param); uwerr.(yarr_lc)

fill_between(xarr[:,2], value.(yarr_ll) .- err.(yarr_ll), value.(yarr_ll) .+ err.(yarr_ll), alpha=0.2, color="forestgreen")
fill_between(xarr[:,2], value.(yarr_lc) .- err.(yarr_lc), value.(yarr_lc) .+ err.(yarr_lc), alpha=0.2, color="royalblue")

# Plotting data points
a_sp = [0.087, 0.077, 0.065, 0.05, 0.039]
color = ["orange", "darkgreen", "magenta", "navy", "royalblue"]

for (k,b) in enumerate(sort(unique(getfield.(ensinfo, :beta))))
    n_ = findall(x->x.beta == b, ensinfo)
    a2_aux = mean(value.(xfit[:,1][n_]))
    uwerr.(xfit)
    errorbar(value.(xfit[n_,2]), value.(yproj_ll[n_]), err.(yproj_ll[n_]), fmt="d",label=string(L"$a \approx$", a_sp[k], " fm"), color=color[k], mfc="none", capsize=2 )
    errorbar(value.(xfit[n_,2]), value.(yproj_lc[n_]), err.(yproj_lc[n_]), fmt="s", color=color[k], mfc="none", capsize=2 )
    # dashed lines
    xxx = [fill(a2_aux, 100) Float64.(range(0.0, 0.8, length=100))]
    yy_ll = fit_pi33_ll(xxx, fit.param)
    println(yy_ll)
    replace!(value.(yy_ll), Inf=>0.0)
    yy_lc = fit_pi33_lc(xxx, fit.param)
    replace!(value.(yy_lc), Inf=>0.0)
    plot(xxx[:,2], value.(yy_ll), ls="--", color=color[k])#dashed lines 
    plot(xxx[:,2], value.(yy_lc), ls="-.", color=color[k])#dashed lines 

end
axvline(ls="dashed", color="black", lw=0.2, alpha=0.7) 
xlabel(L"$\phi_2$")
ylabel(L"$\bar{\Pi}^{33}(-Q^2)$")
ylim(0.018, 0.0235)
legend(ncol=2)
tight_layout()
display(gcf())
tt = "phi2_proj_Q2_9.pdf"
savefig(joinpath(path_plot, tt))
close("all")

## test only cont lim



cl_dep_pi33_ll(x, p) = p[1] .+ p[2] .* x[:,1] .+ p[3] .* x[:,1].^(3/2)
cl_dep_pi33_lc(x, p) = p[1] .+ p[4] .* x[:,1] .+ p[5] .* x[:,1].^(3/2)

fit_tot_pi33_ll1 = (x,p) -> chiral_dep_pi33_ll(x,p) #.+ cl_dep_pi33_ll.(x,p)
fit_tot_pi33_lc1 = (x,p) -> chiral_dep_pi33_lc(x,p) #.+ cl_dep_pi33_lc.(x,p)

## test fit
fix_mom_pi33_ll = getindex.(pi_33_ll, 1)
fix_mom_pi33_lc = getindex.(pi_33_lc, 1)

fit = fit_routine([fit_tot_pi33_ll1, fit_tot_pi33_lc1], [value.(xfit), value.(xfit)], [fix_mom_pi33_ll, fix_mom_pi33_lc], 5 )



