######################################################################################
# This file was created by Alessandro Conigli 
# Here we load the finite-a HVP for the isoscalar component stored in BDIO and perform
# the chiral-continuum extrapolations
########################################################################################

using Revise
using HVPobs
using PyPlot, LaTeXStrings
using OrderedCollections
using TimerOutputs
using BDIO, ADerrors, QuadGK, Statistics

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
path_plot = "/Users/alessandroconigli/MyDrive/postdoc-mainz/projects/deltalpha/plots/isoscalar"

#======= PHYSICAL CONSTANTS ====================#
const t0sqrt_ph = uwreal([0.1443, 0.0007], "sqrtt0 [fm]") 
const MPI_ph = uwreal([134.9768, 0.0005], "mpi phys")
const MK_ph  = uwreal([495.011, 0.01], "mk phys")

const phi2_ph = (sqrt(8)*t0sqrt_ph * MPI_ph / hc)^2
const phi4_ph = (sqrt(8)*t0sqrt_ph)^2 * ((MK_ph/hc)^2 + 0.5*(MPI_ph/hc)^2)

const Qgev = [3., 5., 9.] # Q^2
const Qmgev = 9.0 # Qm^2

#============== READ BDIO FILES =================#
enslist = sort([ "H102", "N101", "C101", "C102", "D150",
         "N451", "D450", "D451", "D452",
         "N203", "N200", "D200", "D201", "E250",
         "J303", "E300",
         ])

dir_path = filter(x-> basename(x) in enslist, readdir(path_bdio, join=true))

##

# deltaiso G33- G88 is stored both in SD and ILD files. No windows was usued in both cases
hvp_path_s1 = vcat(filter(!isempty, [filter(x-> occursin("hvp_set1_SD", x)  , readdir(dir_path[k], join=true)) for k in eachindex(dir_path)])...)
hvp_path_s2 = vcat(filter(!isempty, [filter(x-> occursin("hvp_set2_SD", x)  , readdir(dir_path[k], join=true)) for k in eachindex(dir_path)])...)

spectrum_path = vcat(filter(!isempty, [filter(x-> occursin("spectrum.bdio", x)  , readdir(dir_path[k], join=true)) for k in eachindex(dir_path)])...)

ensinfo = EnsInfo.(enslist)
NENS = length(ensinfo)

# set 1 improvement coefficients
pi_dltiso_ll_s1 = Vector{Vector{uwreal}}(undef, NENS)
pi_dltiso_lc_s1 = Vector{Vector{uwreal}}(undef, NENS)

# set 2 improvement coefficients
pi_dltiso_ll_s2 = Vector{Vector{uwreal}}(undef, NENS)
pi_dltiso_lc_s2 = Vector{Vector{uwreal}}(undef, NENS)

phi2 = Vector{uwreal}(undef, NENS)
phi4 = Vector{uwreal}(undef, NENS)
a28t0 = Vector{uwreal}(undef, NENS)

for (k, ens) in enumerate(ensinfo)
    println("- Ensemble ", ens.id)

    # set 1 improvement coefficients
    pi_dltiso_ll_s1[k] = read_BDIO(hvp_path_s1[k], "dalpha", "pi_88_ll_conn") 
    pi_dltiso_lc_s1[k] = read_BDIO(hvp_path_s1[k], "dalpha", "pi_88_lc_conn") 

    # set 2 improvement coefficients
    pi_dltiso_ll_s2[k] = read_BDIO(hvp_path_s2[k], "dalpha", "pi_88_ll_conn" ) 
    pi_dltiso_lc_s2[k] = read_BDIO(hvp_path_s2[k], "dalpha", "pi_88_lc_conn" ) 

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
NMOM = length(pi_dltiso_ll_s1[1])
# set 1 pi88
fitcat_dltiso_ll_s1 = [Vector{FitCat}(undef,0) for i=1:NMOM]
fitcat_dltiso_lc_s1 = [Vector{FitCat}(undef,0) for i=1:NMOM]
# set 2 pi88
fitcat_dltiso_ll_s2 = [Vector{FitCat}(undef,0) for i=1:NMOM]
fitcat_dltiso_lc_s2 = [Vector{FitCat}(undef,0) for i=1:NMOM]

## no cuts in data
for s in 1:2
    xdata = [a28t0 phi2 phi4]
    if s == 1
        for q in 1:length(pi_dltiso_ll_s1[1])
            str = "all_data_set$(s)_q$(q)"
            # pi 33
            push!(fitcat_dltiso_ll_s1[q], FitCat(xdata, getindex.(pi_dltiso_ll_s1, q), str))
            push!(fitcat_dltiso_lc_s1[q], FitCat(xdata, getindex.(pi_dltiso_lc_s1, q), str))
        end
    elseif s == 2
        for q in 1:length(pi_dltiso_ll_s2[1])
            str = "all_data_set$(s)_q$(q)"
            # pi 33
            push!(fitcat_dltiso_ll_s2[q], FitCat(xdata, getindex.(pi_dltiso_ll_s2, q), str))
            push!(fitcat_dltiso_lc_s2[q], FitCat(xdata, getindex.(pi_dltiso_lc_s2, q), str))
        end
    end
end


# pi 88 
for q in 1:NMOM

    for k_cat in eachindex(fitcat_dltiso_ll_s1[q])
        xdata = fitcat_dltiso_ll_s1[q][k_cat].xdata
        ydata_ll_s1 = fitcat_dltiso_ll_s1[q][k_cat].ydata
        ydata_ll_s2 = fitcat_dltiso_ll_s2[q][k_cat].ydata
        ydata_lc_s1 = fitcat_dltiso_lc_s1[q][k_cat].ydata
        ydata_lc_s2 = fitcat_dltiso_lc_s2[q][k_cat].ydata

        for (k_mod, model) in enumerate(f_tot_dltiso)
            println(k_mod)
            fit_ll_s1 = fit_routine(model, value.(xdata), ydata_ll_s1, n_par_tot_dltiso[k_mod])
            fit_ll_s2 = fit_routine(model, value.(xdata), ydata_ll_s2, n_par_tot_dltiso[k_mod])
            fit_lc_s1 = fit_routine(model, value.(xdata), ydata_lc_s1, n_par_tot_dltiso[k_mod])
            fit_lc_s2 = fit_routine(model, value.(xdata), ydata_lc_s2, n_par_tot_dltiso[k_mod])
            push!(fitcat_dltiso_ll_s1[q][k_cat].fit, fit_ll_s1)
            push!(fitcat_dltiso_ll_s2[q][k_cat].fit, fit_ll_s2)
            push!(fitcat_dltiso_lc_s1[q][k_cat].fit, fit_lc_s1)
            push!(fitcat_dltiso_lc_s2[q][k_cat].fit, fit_lc_s2)

        end
    end
end

## plot all cl limits
plot_cl_all_set(fitcat_dltiso_ll_s1, fitcat_dltiso_ll_s2, fitcat_dltiso_lc_s1, fitcat_dltiso_lc_s2, path_plot=path_plot, ylab=L"$\Delta_{ls}(\Delta\alpha)$", f_tot_isov=f_tot_dltiso)
plot_chiral_best_fit(fitcat_dltiso_lc_s2, path_plot=path_plot, tt=["Set", "2", "LC"], f_tot_isov=f_tot_dltiso, ylab=L"$\Delta_{ls}(\Delta\alpha)$")
plot_cl_best_fit(fitcat_dltiso_lc_s1, path_plot=path_plot, tt=["Set", "1", "LC"], f_tot_isov=f_tot_dltiso, ylab=L"$\Delta_{ls}(\Delta\alpha)$")



###################################
## RESULTS
##################################


for q in 1:NMOM
    @info "Momentum no. $(q)"
    fitcat_pi88_tot = vcat(vcat(fitcat_dltiso_ll_s1[q],
                fitcat_dltiso_ll_s2[q],
                fitcat_dltiso_lc_s1[q],
                fitcat_dltiso_lc_s2[q])...)

    ww_tot = get_w_from_fitcat(fitcat_pi88_tot)

    w, widx  =  findmax(ww_tot)
  
    model_idx = mod(widx, length(f_tot_dltiso))
    if model_idx == 0 
        model_idx = length(f_tot_dltiso)
    end
    println("   wmax: ", w, "widx: ", widx)
    model = f_tot_dltiso[model_idx]
    
    
    cat_idx = Int((widx - model_idx ) / length(f_tot_dltiso))+1
    if cat_idx < 0
        cat_idx +=1
    end
    println("   best χ2/χ2exp: ", fitcat_pi88_tot[cat_idx].fit[model_idx].chi2 / fitcat_pi88_tot[cat_idx].fit[model_idx].chi2exp)
    println("   widx: ", widx, " model_idx: ", model_idx, " catidx: ", cat_idx)
    ## Best Res
    best_mod = f_tot_dltiso[model_idx]
    xdata = fitcat_pi88_tot[cat_idx].xdata
    param = fitcat_pi88_tot[cat_idx].fit[model_idx].param

    ph_res_best = best_mod([0.0 phi2_ph phi4_ph], param)[1]; uwerr(ph_res_best)
    println("   best res: ", ph_res_best )
    ## histogram

    all_res = Vector{uwreal}()
    for (k, cat) in enumerate(fitcat_pi88_tot)
        for (j, mod) in enumerate(f_tot_dltiso)
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
