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
include("./func_comb.jl")

path_bdio_obs = "/Users/alessandroconigli/MyDrive/postdoc-mainz/projects/deltalpha/data"
path_store_pi = "/Users/alessandroconigli/MyDrive/postdoc-mainz/projects/deltalpha/PIdata/impr_deriv/high_q_kernel/std_mom/"
path_plot = "/Users/alessandroconigli/MyDrive/postdoc-mainz/projects/deltalpha/plots/c8Disc"


#======= PHYSICAL CONSTANTS ====================#
const MPI_ph = uwreal([134.9768, 0.0005], "mpi phys")
const MK_ph  = uwreal([495.011, 0.01], "mk phys")

const phi2_ph = (sqrt(8)*t0sqrt_ph * MPI_ph / hc)^2
const phi4_ph = (sqrt(8)*t0sqrt_ph)^2 * ((MK_ph/hc)^2 + 0.5*(MPI_ph/hc)^2)

const Qgev = [3., 5., 9.] # Q^2
const Qmgev = 9.0 # Qm^2


enslist = sort([  "H102", "N101", "C101", "C102", "D150",
        "N451", "D451", "D452",
        "N203", "N200", "D200", "D201", "E250",
        "J303", "E300",
        "J501"])


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

#============= READ PI CHARM 8 DISC FROM BDIO =================#
# set 1 improvement coefficients -> SD window included 
pi_c8_cc_s1 = Vector{Vector{uwreal}}(undef, 0)

# set 2 improvement coefficients -> No window included
pi_c8_cc_s2 = Vector{Vector{uwreal}}(undef, 0)

fb = BDIO_open(joinpath(path_store_pi, "PIc8_disc.bdio"), "r")
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
    push!(pi_c8_cc_s1, res["pic8_cc_s1"])
    push!(pi_c8_cc_s2, res["pic8_cc_s2"])

end
BDIO_close!(fb)


##############################
## CREATE FIT CATEGORIES
##############################
NMOM = length(Qgev)
# set 1 picc
fitcat_c8_cc_s1 = [Vector{FitCat}(undef,0) for i=1:NMOM]
# set 2 picc
fitcat_c8_cc_s2 = [Vector{FitCat}(undef,0) for i=1:NMOM]

## no cuts in data
for s in 1:2
    xdata = [a28t0 phi2 phi4]
    if s == 1
        for q in 1:length(pi_c8_cc_s1[1])
            str = "all_data_set$(s)_q$(q)"
            # pi 33
            push!(fitcat_c8_cc_s1[q], FitCat(xdata, getindex.(pi_c8_cc_s1, q), str))
        end
    elseif s == 2
        for q in 1:length(pi_c8_cc_s2[1])
            str = "all_data_set$(s)_q$(q)"
            # pi 33
            push!(fitcat_c8_cc_s2[q], FitCat(xdata, getindex.(pi_c8_cc_s2, q), str))
        end
    end
end
##

#================= FITTING ====================#
# pi cc connected
for q in 1:NMOM

    for k_cat in eachindex(fitcat_c8_cc_s1[q])
        xdata = fitcat_c8_cc_s1[q][k_cat].xdata
        ydata_cc_s1 = fitcat_c8_cc_s1[q][k_cat].ydata
        ydata_cc_s2 = fitcat_c8_cc_s2[q][k_cat].ydata

        for (k_mod, model) in enumerate(f_tot_dltiso)
            println(k_mod)
            fit_cc_s1 = fit_routine(model, value.(xdata), ydata_cc_s1, n_par_tot_dltiso[k_mod], pval=true)
            fit_cc_s2 = fit_routine(model, value.(xdata), ydata_cc_s2, n_par_tot_dltiso[k_mod], pval=true)
            push!(fitcat_c8_cc_s1[q][k_cat].fit, fit_cc_s1)
            push!(fitcat_c8_cc_s2[q][k_cat].fit, fit_cc_s2)

        end
    end
end


##########################
## PLOTS
#########################
using Statistics
plot_cl_all_set(fitcat_cc_ll_s1, fitcat_cc_ll_s2, fitcat_cc_lc_s1, fitcat_cc_lc_s2, path_plot=path_plot, ylab=L"$(\Delta\alpha^{c,8})_{\mathrm{sub}}$", f_tot_isov=f_tot_charm)
plot_chiral_best_fit(fitcat_c8_cc_s1, path_plot=path_plot, tt=[""], f_tot_isov=f_tot_dltiso, ylab=L"$(\Delta\alpha^{c,8})_{\mathrm{disc}}$")
plot_cl_best_fit(fitcat_c8_cc_s2, path_plot=path_plot, tt=["No Wind", "CC"], f_tot_isov=f_tot_dltiso, ylab=L"$(\Delta\alpha^{c,8})_{\mathrm{disc}}$")

cattot = [vcat(fitcat_cc_lc_s1[k], fitcat_cc_lc_s2[k]...) for k in eachindex(fitcat_cc_lc_s1)]
plot_mAve_summary(cattot, xlab=vcat(label_tot_charm,label_tot_charm), charge_factor=2/(3*sqrt(3)), ylab=L"$\frac{2}{3\sqrt{3}}(\Delta\alpha^{c,8})_{\mathrm{sub}}$", tt=["Set", "1-2", "LC"], path_plot=path_plot)
plot_mAve_summary(fitcat_c8_cc_s1, xlab=label_tot_dltiso, charge_factor=2/(3*sqrt(3)), tt=["SD Wind", "CC"],  ylab=L"$\frac{2}{3\sqrt{3}}(\Delta\alpha^{c,8})_{\mathrm{disc}}$", path_plot=path_plot)


###################################
## RESULTS
##################################

RES = []
SYST = []
for q in 1:NMOM
    @info "Momentum no. $(q)"
    fitcat_c8_tot = vcat(vcat(#fitcat_cc_ll_s1[q],
                #fitcat_cc_ll_s2[q],
                fitcat_c8_cc_s1[q])...)

    ww_tot = get_w_from_fitcat(fitcat_c8_tot)

    w, widx  =  findmax(ww_tot)
  
    model_idx = mod(widx, length(f_tot_dltiso))
    if model_idx == 0 
        model_idx = length(f_tot_dltiso)
    end
    println("   wmax: ", w, " model_idx: ", model_idx)
    model = f_tot_dltiso[model_idx]
    
    
    cat_idx = Int((widx - model_idx ) / length(f_tot_dltiso))+1
    if cat_idx < 0
        cat_idx +=1
    end
    println("   best χ2/χ2exp: ", fitcat_c8_tot[cat_idx].fit[model_idx].chi2 / fitcat_c8_tot[cat_idx].fit[model_idx].chi2exp)
    println("   widx: ", widx, " model_idx: ", model_idx, " catidx: ", cat_idx)
    ## Best Res
    best_mod = f_tot_dltiso[model_idx]
    xdata = fitcat_c8_tot[cat_idx].xdata
    param = fitcat_c8_tot[cat_idx].fit[model_idx].param

    ph_res_best = 2 ./ (3 .* sqrt(3)) .* 1e5 * best_mod([0.0 phi2_ph phi4_ph], param)[1]; uwerr(ph_res_best)
    println("   best res: ", ph_res_best )
    ## histogram

    all_res = Vector{uwreal}()
    for (k, cat) in enumerate(fitcat_c8_tot)
        for (j, mod) in enumerate(f_tot_dltiso)
            push!(all_res, mod([0.0 phi2_ph phi4_ph], cat.fit[j].param)[1])
            # push!(all_res, mod([0.0 value(phi2_ph) value(phi4_ph)], cat.fit[j].param)[1])
        end
    end

    final_res, syst = 2 ./ (3 .* sqrt(3)) .* 1e5 .* model_average(all_res, ww_tot); uwerr(final_res)
    push!(RES, final_res)
    push!(SYST, syst)
    println("   Model ave:  ", final_res)
    println("   systematic: ", syst)


    #hist(value.(all_res), bins=800, histtype="stepfilled", alpha=0.5, ec="k", color="navy", weights=ww_tot)
    #display(gcf())
    #close()

end