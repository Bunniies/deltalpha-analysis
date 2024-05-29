using Revise, HVPobs
using PyPlot, LaTeXStrings
using BDIO, ADerrors, ALPHAio
import ADerrors: err

#================ SET UP PLOT PARAMETERS ==================#

rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
rcParams["text.usetex"] =  true
rcParams["mathtext.fontset"]  = "cm"
rcParams["font.size"] =13
rcParams["axes.labelsize"] =22
rcParams["axes.titlesize"] = 18
plt.rc("text", usetex=true) # set to true if a LaTeX installation is present


include("../utils/const.jl")
include("../utils/types.jl")
include("../utils/plot_utils.jl")
include("../utils/IO_BDIO.jl")
include("../utils/tools.jl")
include("func_comb_PI33.jl")

path_bdio_obs = "/Users/alessandroconigli/MyDrive/postdoc-mainz/projects/deltalpha/data"
path_store_pi = "/Users/alessandroconigli/MyDrive/postdoc-mainz/projects/deltalpha/PIdata/impr_deriv"
path_plot = "/Users/alessandroconigli/MyDrive/postdoc-mainz/projects/deltalpha/plots/isovector/ILD"
path_phys_res = "/Users/alessandroconigli/MyDrive/postdoc-mainz/projects/deltalpha/physical_results"

#======= PHYSICAL CONSTANTS ====================#
const MPI_ph = uwreal([134.9768, 0.0005], "mpi phys")
const MK_ph  = uwreal([495.011, 0.01], "mk phys")

const phi2_ph = (sqrt(8)*t0sqrt_ph * MPI_ph / hc)^2
const phi4_ph = (sqrt(8)*t0sqrt_ph)^2 * ((MK_ph/hc)^2 + 0.5*(MPI_ph/hc)^2)

const Qgev = [3., 5., 9.] # Q^2
const Qmgev = 9.0 # Qm^2


#============== READ CORRELATORS FROM BDIO FILES =================#

enslist = sort([ "H101", "H102", "N101", "C101", "C102", "D150",
          "B450", "N451", "D450", "D451", "D452",
         "N202", "N203", "N200", "D200", "D201", "E250",
          "N300", "J303", "E300",
         "J500", "J501"])

ensinfo = EnsInfo.(enslist)

#============== READ PS MASSES AND t0 FROM BDIO FILES =================#
dir_path = filter(x-> basename(x) in enslist, readdir(path_bdio_obs, join=true))
spectrum_path = vcat(filter(!isempty, [filter(x-> occursin("spectrum.bdio", x)  , readdir(dir_path[k], join=true)) for k in eachindex(dir_path)])...)

phi2 = Vector{uwreal}(undef, length(enslist))
phi4 = Vector{uwreal}(undef, length(enslist))
a28t0 = Vector{uwreal}(undef, length(enslist))
t0ens = Vector{uwreal}(undef, length(enslist))

for (k, ens) in enumerate(ensinfo)
    m_pi  = read_BDIO(spectrum_path[k], "spectrum", "mpi")[1]
    m_k   = read_BDIO(spectrum_path[k], "spectrum", "mk")[1]
    t0ens_aux = read_BDIO(spectrum_path[k], "spectrum", "t0")[1]

    phi2[k]  = 8 * t0ens_aux * m_pi^2
    phi4[k]  = 8 * t0ens_aux * (m_k^2 + 0.5*m_pi^2 )
    a28t0[k] =  a(ens.beta)^2 #1 / (8*t0ens_aux)
    t0ens[k] = t0ens_aux
end

#============= READ PI 33  FROM BDIO =================#
# set 1 improvement coefficients
pi_33_ll_ILD_s1 = Vector{Vector{uwreal}}(undef, 0)
pi_33_lc_ILD_s1 = Vector{Vector{uwreal}}(undef, 0)

# set 2 improvement coefficients
pi_33_ll_ILD_s2 = Vector{Vector{uwreal}}(undef, 0)
pi_33_lc_ILD_s2 = Vector{Vector{uwreal}}(undef, 0)

fb = BDIO_open(joinpath(path_store_pi, "PI_33.bdio"), "r")
res = Dict()
count=0
while ALPHAdobs_next_p(fb)
    count+=1
    d = ALPHAdobs_read_parameters(fb)
    sz = tuple(d["size"]...)
    extra = d["extra"]
    ks = collect(d["keys"])
    res = ALPHAdobs_read_next(fb, size=sz, keys=ks)
    if extra["Ens"] ∉ enslist # != enslist[count]
        @info("Mismatch with EnsID in current BDIO uinfo2 ")
        println(extra["Ens"], " ", enslist)
        continue
    end

    push!(pi_33_ll_ILD_s1, res["pi33_ll_ILD_s1"])
    push!(pi_33_lc_ILD_s1, res["pi33_lc_ILD_s1"])
    push!(pi_33_ll_ILD_s2, res["pi33_ll_ILD_s2"])
    push!(pi_33_lc_ILD_s2, res["pi33_lc_ILD_s2"])

end
BDIO_close!(fb)

##

##############################
## CREATE FIT CATEGORIES
##############################
NMOM = length(Qgev)
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
        for q in 1:NMOM
            str = "all_data_set$(s)_q$(q)"
            # pi 33
            push!(fitcat_pi33_ll_s1[q], FitCat(xdata, getindex.(pi_33_ll_ILD_s1, q), str))
            push!(fitcat_pi33_lc_s1[q], FitCat(xdata, getindex.(pi_33_lc_ILD_s1, q), str))
        end
    elseif s == 2
        for q in 1:NMOM
            str = "all_data_set$(s)_q$(q)"
            # pi 33
            push!(fitcat_pi33_ll_s2[q], FitCat(xdata, getindex.(pi_33_ll_ILD_s2, q), str))
            push!(fitcat_pi33_lc_s2[q], FitCat(xdata, getindex.(pi_33_lc_ILD_s2, q), str))
        end

    end
end
## cuts in phi2<0.6 
i_cutphi2 = findall(x->x<0.6, value.(phi2))
for s in 1:2
    xdata = [a28t0[i_cutphi2] phi2[i_cutphi2] phi4[i_cutphi2]]
    if s == 1
        for q in 1:length(NMOM)
            str = "all_data_set$(s)_q$(q)"
            # pi 33
            push!(fitcat_pi33_ll_s1[q], FitCat(xdata, getindex.(pi_33_ll_ILD_s1, q)[i_cutphi2], str))
            push!(fitcat_pi33_lc_s1[q], FitCat(xdata, getindex.(pi_33_lc_ILD_s1, q)[i_cutphi2], str))
        end
    elseif s == 2
        for q in 1:length(NMOM)
            str = "all_data_set$(s)_q$(q)"
            # pi 33
            push!(fitcat_pi33_ll_s2[q], FitCat(xdata, getindex.(pi_33_ll_ILD_s2, q)[i_cutphi2], str))
            push!(fitcat_pi33_lc_s2[q], FitCat(xdata, getindex.(pi_33_lc_ILD_s2, q)[i_cutphi2], str))
        end

    end
end

##
#================= FITTING ====================#
# pi 33 
for q in 1:NMOM

    for k_cat in eachindex(fitcat_pi33_ll_s1[q])
        xdata = fitcat_pi33_ll_s1[q][k_cat].xdata
        ydata_ll_s1 = fitcat_pi33_ll_s1[q][k_cat].ydata
        ydata_ll_s2 = fitcat_pi33_ll_s2[q][k_cat].ydata
        ydata_lc_s1 = fitcat_pi33_lc_s1[q][k_cat].ydata
        ydata_lc_s2 = fitcat_pi33_lc_s2[q][k_cat].ydata

        for (k_mod, model) in enumerate(f_tot_isov)
            println(k_mod)
            fit_ll_s1 = fit_routine(model, value.(xdata), ydata_ll_s1, n_par_tot_isov[k_mod], pval=true)
            fit_ll_s2 = fit_routine(model, value.(xdata), ydata_ll_s2, n_par_tot_isov[k_mod], pval=true)
            fit_lc_s1 = fit_routine(model, value.(xdata), ydata_lc_s1, n_par_tot_isov[k_mod], pval=true)
            fit_lc_s2 = fit_routine(model, value.(xdata), ydata_lc_s2, n_par_tot_isov[k_mod], pval=true)
            push!(fitcat_pi33_ll_s1[q][k_cat].fit, fit_ll_s1)
            push!(fitcat_pi33_ll_s2[q][k_cat].fit, fit_ll_s2)
            push!(fitcat_pi33_lc_s1[q][k_cat].fit, fit_lc_s1)
            push!(fitcat_pi33_lc_s2[q][k_cat].fit, fit_lc_s2)

        end
    end
end

##########################
## PLOTS
#########################
using Statistics
plot_cl_all_set(fitcat_pi33_ll_s1, fitcat_pi33_ll_s2, fitcat_pi33_lc_s1, fitcat_pi33_lc_s2, path_plot=path_plot, ylab=L"$(\Delta\alpha^{3,3})_{\mathrm{sub}}^{\mathrm{ILD}}$", f_tot_isov=f_tot_isov)
plot_chiral_best_fit(fitcat_pi33_ll_s2, path_plot=path_plot, tt=["Set", "2", "LL"], f_tot_isov=f_tot_isov, ylab=L"$(\Delta\alpha^{3,3})_{\mathrm{sub}}^{\mathrm{ILD}}$")
plot_cl_best_fit(fitcat_pi33_lc_s2, path_plot=path_plot, tt=["Set", "2", "LC"], f_tot_isov=f_tot_isov, ylab=L"$(\Delta\alpha^{3,3})_{\mathrm{sub}}^{\mathrm{ILD}}$")

cattot = [vcat(fitcat_pi33_ll_s1[k], fitcat_pi33_lc_s1[k], fitcat_pi33_ll_s2[k],fitcat_pi33_lc_s2[k]...) for k in eachindex(fitcat_pi33_lc_s2)]
plot_mAve_summary(cattot, xlab=vcat(label_tot_isov,label_tot_isov,label_tot_isov,label_tot_isov), charge_factor=1., ylab=L"$(\Delta\alpha^{3,3})_{\mathrm{sub}}^{\mathrm{ILD}}$", tt=["Set", "1-2", "LL-LC"], path_plot=path_plot)
plot_mAve_summary(fitcat_pi33_lc_s1, xlab=label_tot_isov, charge_factor=1., ylab=L"$(\Delta\alpha^{3,3})_{\mathrm{sub}}^{\mathrm{ILD}}$")


###################################
## RESULTS
##################################

RES = []
SYST = []
for q in 1:NMOM
    @info "Momentum no. $(q)"
    fitcat_pi33_tot = vcat(vcat(fitcat_pi33_ll_s1[q],
                fitcat_pi33_ll_s2[q],
                fitcat_pi33_lc_s1[q],
                fitcat_pi33_lc_s2[q])...)

    ww_tot = get_w_from_fitcat(fitcat_pi33_tot)

    w, widx  =  findmax(ww_tot)
  
    model_idx = mod(widx, length(f_tot_isov))
    if model_idx == 0 
        model_idx = length(f_tot_isov)
    end
    println("   wmax: ", w, " model_idx: ", model_idx)
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

    ph_res_best = best_mod([0.0 value(phi2_ph) value(phi4_ph)], param)[1]; uwerr(ph_res_best)
    println("   best res: ", ph_res_best )
    ## histogram

    all_res = Vector{uwreal}()
    for (k, cat) in enumerate(fitcat_pi33_tot)
        for (j, mod) in enumerate(f_tot_isov)
            push!(all_res, mod([0.0 value(phi2_ph) value(phi4_ph)], cat.fit[j].param)[1])
        end
    end

    final_res, syst =  model_average(all_res, ww_tot); uwerr(final_res)
    push!(RES, final_res)
    push!(SYST, syst)
    println("   Model ave:  ", final_res)
    println("   systematic: ", syst)


    #hist(value.(all_res), bins=800, histtype="stepfilled", alpha=0.5, ec="k", color="navy", weights=ww_tot)
    #display(gcf())
    #close()

end

## saving physical results in BDIO

# aaa  = RES[1] + uwreal([0.0, SYST[1]], "SYST")

io = IOBuffer()
write(io, "PI33 ILD physical results")
fb = ALPHAdobs_create(joinpath(path_phys_res, "PI33_ILD_physRes.bdio"), io)
for k in eachindex(RES)
    aux = RES[k] + uwreal([0.0, SYST[k]], "Syst Pi33 ILD")
    ALPHAdobs_write(fb, aux)
end
ALPHAdobs_close(fb)

## test reading
fb = BDIO_open(joinpath(path_phys_res, "PI33_ILD_physRes.bdio"), "r")
res = []
while ALPHAdobs_next_p(fb)
    d = ALPHAdobs_read_parameters(fb)
    push!(res, ALPHAdobs_read_next(fb))
end
BDIO_close!(fb)