using Revise, HVPobs
using PyPlot, LaTeXStrings
using BDIO, ADerrors, ALPHAio
import ADerrors: err


#================ SET UP PLOT PARAMETERS ==================#

rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
rcParams["text.usetex"] =  true
rcParams["mathtext.fontset"]  = "cm"
rcParams["font.size"] =20
rcParams["axes.labelsize"] =26
rcParams["axes.titlesize"] = 22
plt.rc("text", usetex=true) # set to true if a LaTeX installation is present


include("../../utils/const.jl")
include("../../utils/types.jl")
include("../../utils/plot_utils.jl")
include("../../utils/IO_BDIO.jl")
include("../../utils/tools.jl")
include("./func_comb_charm_conn.jl")

path_bdio_obs = "/Users/alessandroconigli/MyDrive/postdoc-mainz/projects/deltalpha/data"
path_store_pi = "/Users/alessandroconigli/MyDrive/postdoc-mainz/projects/deltalpha/PIdata/impr_deriv/new_strategy/Q_LD/"
path_plot = "/Users/alessandroconigli/MyDrive/postdoc-mainz/projects/deltalpha/plots/new_strategy/picc/Q_LD/"
path_phys_res = "/Users/alessandroconigli/MyDrive/postdoc-mainz/projects/deltalpha/physical_results/new_strategy/Q_LD/"


#======= PHYSICAL CONSTANTS ====================#
const MPI_ph = uwreal([134.9768, 0.0005], "mpi phys")
const MK_ph  = uwreal([495.011, 0.01], "mk phys")

const phi2_ph = (sqrt(8)*t0sqrt_ph * MPI_ph / hc)^2
const phi4_ph = (sqrt(8)*t0sqrt_ph)^2 * ((MK_ph/hc)^2 + 0.5*(MPI_ph/hc)^2)

# const Qgev = [3., 5., 9.] # Q^2
const Qgev = [4, 5, 6, 7, 8, 9, 12] ./ 16 # Q^2/4 # additional very high values

const Qmgev = 9.0 # Qm^2

enslist = sort([ "H101", "H102", "N101", "C101",
        "D450", "D452", # B450  removed
         "N202", "N203", "N200", "D200",  "E250",
        "N300", "J303", "E300", #  removed
         "J500"
])

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

    # println("using t0/a^2 for each ens")
    # phi2[k]  = 8 * t0ens_aux * m_pi^2
    # phi4[k]  = 8 * t0ens_aux * (m_k^2 + 0.5*m_pi^2 )
    println("using t0/a^2 from t0_sym Regensburg")
    phi2[k]  = 8 * t0(ens.beta) * m_pi^2
    phi4[k]  = 8 * t0(ens.beta) * (m_k^2 + 0.5*m_pi^2 )

    a28t0[k] =  a(ens.beta)^2 #1 / (8*t0ens_aux)
    t0ens[k] = t0ens_aux
end

#============= READ PI CHARM CHARM FROM BDIO =================#
# set 1 improvement coefficients
pi_cc_ll_s1 = Vector{Vector{uwreal}}(undef, 0)
pi_cc_lc_s1 = Vector{Vector{uwreal}}(undef, 0)

# set 2 improvement coefficients
pi_cc_ll_s2 = Vector{Vector{uwreal}}(undef, 0)
pi_cc_lc_s2 = Vector{Vector{uwreal}}(undef, 0)

# @warn "you are reading data not interpolated to target k_c! Should read PIcc_conn_interp.bdio! "
fb = BDIO_open(joinpath(path_store_pi, "PIcc_conn_interp.bdio"), "r")
tmp_res = Dict{Any, Any}()
res = Dict()
count=0
while ALPHAdobs_next_p(fb)
    count+=1
    d = ALPHAdobs_read_parameters(fb)
    sz = tuple(d["size"]...)
    extra = d["extra"]
    ks = collect(d["keys"])
    if extra["Ens"] ∉ enslist # != enslist[count]
        @info("Mismatch with EnsID in current BDIO uinfo2 ")
        println(extra["Ens"], " ", enslist)
        continue
    end
    res = ALPHAdobs_read_next(fb, size=sz, keys=ks)
    tmp_res[extra["Ens"]] = Dict{String, Array{uwreal}}(
        "picc_ll_s1" => res["picc_ll_s1"], 
        "picc_lc_s1" => res["picc_lc_s1"], 
        "picc_ll_s2" => res["picc_ll_s2"], 
        "picc_lc_s2" => res["picc_lc_s2"] 
    )

end
BDIO_close!(fb)

# rearrange the order 

for (k, ens) in enumerate(enslist)
    push!(pi_cc_ll_s1, tmp_res[ens]["picc_ll_s1"])
    push!(pi_cc_lc_s1, tmp_res[ens]["picc_lc_s1"])
    push!(pi_cc_ll_s2, tmp_res[ens]["picc_ll_s2"])
    push!(pi_cc_lc_s2, tmp_res[ens]["picc_lc_s2"])

end

##
##############################
## CREATE FIT CATEGORIES
##############################
NMOM = length(Qgev)
# set 1 picc
fitcat_cc_ll_s1 = [Vector{FitCat}(undef,0) for i=1:NMOM]
fitcat_cc_lc_s1 = [Vector{FitCat}(undef,0) for i=1:NMOM]
# set 2 picc
fitcat_cc_ll_s2 = [Vector{FitCat}(undef,0) for i=1:NMOM]
fitcat_cc_lc_s2 = [Vector{FitCat}(undef,0) for i=1:NMOM]

## no cuts in data
for s in 1:2
    xdata = [a28t0 phi2 phi4]
    if s == 1
        for q in 1:NMOM
            str = "all_data_set$(s)_q$(q)"
            # pi 33
            push!(fitcat_cc_ll_s1[q], FitCat(xdata, getindex.(pi_cc_ll_s1, q), str))
            push!(fitcat_cc_lc_s1[q], FitCat(xdata, getindex.(pi_cc_lc_s1, q), str))
        end
    elseif s == 2
        for q in 1:NMOM
            str = "all_data_set$(s)_q$(q)"
            # pi 33
            push!(fitcat_cc_ll_s2[q], FitCat(xdata, getindex.(pi_cc_ll_s2, q), str))
            push!(fitcat_cc_lc_s2[q], FitCat(xdata, getindex.(pi_cc_lc_s2, q), str))
        end
    end
end
## cuts in beta > 3.4
i_cutbeta = findall(x->x>3.4, getfield.(ensinfo, :beta))

for s in 1:2
    xdata = [a28t0[i_cutbeta] phi2[i_cutbeta] phi4[i_cutbeta]]
    if s == 1
        for q in 1:NMOM
            str = "betacut_set$(s)_q$(q)"
            # pi 33
            push!(fitcat_cc_ll_s1[q], FitCat(xdata, getindex.(pi_cc_ll_s1, q)[i_cutbeta], str))
            push!(fitcat_cc_lc_s1[q], FitCat(xdata, getindex.(pi_cc_lc_s1, q)[i_cutbeta], str))
        end
    elseif s == 2
        for q in 1:NMOM
            str = "betacut_set$(s)_q$(q)"
            # pi 33
            push!(fitcat_cc_ll_s2[q], FitCat(xdata, getindex.(pi_cc_ll_s2, q)[i_cutbeta], str))
            push!(fitcat_cc_lc_s2[q], FitCat(xdata, getindex.(pi_cc_lc_s2, q)[i_cutbeta], str))
        end

    end
end
##
#================= FITTING ====================#
# pi cc connected
for q in 1:NMOM

    for k_cat in eachindex(fitcat_cc_ll_s1[q])
        xdata = fitcat_cc_ll_s1[q][k_cat].xdata
        ydata_ll_s1 = fitcat_cc_ll_s1[q][k_cat].ydata
        ydata_ll_s2 = fitcat_cc_ll_s2[q][k_cat].ydata
        ydata_lc_s1 = fitcat_cc_lc_s1[q][k_cat].ydata
        ydata_lc_s2 = fitcat_cc_lc_s2[q][k_cat].ydata

        for (k_mod, model) in enumerate(f_tot_charm)
            println(k_mod)
            try
                fit_ll_s1 = fit_routine(model, value.(xdata), ydata_ll_s1, n_par_tot_charm[k_mod], pval=false)
                fit_ll_s2 = fit_routine(model, value.(xdata), ydata_ll_s2, n_par_tot_charm[k_mod], pval=false)
                fit_lc_s1 = fit_routine(model, value.(xdata), ydata_lc_s1, n_par_tot_charm[k_mod], pval=false)
                fit_lc_s2 = fit_routine(model, value.(xdata), ydata_lc_s2, n_par_tot_charm[k_mod], pval=false)
                
                # fit_ll_s1 = fit_data(model, value.(xdata), ydata_ll_s1, n_par_tot_charm[k_mod], correlated=false)
                # fit_ll_s2 = fit_data(model, value.(xdata), ydata_ll_s2, n_par_tot_charm[k_mod], correlated=false)
                # fit_lc_s1 = fit_data(model, value.(xdata), ydata_lc_s1, n_par_tot_charm[k_mod], correlated=false)
                # fit_lc_s2 = fit_data(model, value.(xdata), ydata_lc_s2, n_par_tot_charm[k_mod], correlated=false)
                
                push!(fitcat_cc_ll_s1[q][k_cat].fit, fit_ll_s1)
                push!(fitcat_cc_ll_s2[q][k_cat].fit, fit_ll_s2)
                push!(fitcat_cc_lc_s1[q][k_cat].fit, fit_lc_s1)
                push!(fitcat_cc_lc_s2[q][k_cat].fit, fit_lc_s2)
            catch
                println("failed")
            end
        end
    end
end

##########################
## PLOTS
#########################
using Statistics
ll = L"$\bar{\Pi}^{(c,c)}(Q^2/16)$"
plot_cl_all_set(fitcat_cc_ll_s1, fitcat_cc_ll_s2, fitcat_cc_lc_s1, fitcat_cc_lc_s2, nmom=3, path_plot=path_plot, ylab=ll, f_tot_isov=f_tot_charm)
plot_chiral_best_fit(fitcat_cc_lc_s1, path_plot=path_plot, tt=["Set", "1", "LC"], f_tot_isov=f_tot_charm, ylab=ll)
plot_cl_best_fit(fitcat_cc_lc_s2, path_plot=path_plot, tt=["Set", "2", "LC"], f_tot_isov=f_tot_charm, ylab=L"$(\Delta\alpha^{c,c})_{\mathrm{sub}}$")

cattot = [vcat(fitcat_cc_lc_s1[k]...) for k in eachindex(fitcat_cc_lc_s1)]
plot_mAve_summary(cattot, xlab=vcat(label_tot_charm,label_tot_charm), charge_factor=4/9, ylab=L"$\frac{4}{9}(\Delta\alpha^{c,c})_{\mathrm{sub}}$", tt=["Set", "1-2", "LC"], path_plot=nothing)
plot_mAve_summary(fitcat_cc_lc_s1, xlab=vcat(label_tot_charm,label_tot_charm), charge_factor=4/9, ylab=L"$\frac{4}{9}(\Delta\alpha^{c,c})_{\mathrm{sub}}$")

###################################
## RESULTS
##################################

RES = []
SYST = []
for q in 1:NMOM
    @info "Momentum no. $(q): $(Qgev[q]) GeV^2"
    fitcat_cc_mean = vcat(vcat(
                fitcat_cc_lc_s1[q],
                fitcat_cc_lc_s2[q])...
    )
    fitcat_cc_syst = vcat(vcat(fitcat_cc_lc_s1[q],
                    fitcat_cc_lc_s2[q])...
    )
    # ww_tot = get_w_from_fitcat(fitcat_cc_tot)

    ww_lc_s1 = get_w_from_fitcat(fitcat_cc_lc_s1[q])
    ww_lc_s2 = get_w_from_fitcat(fitcat_cc_lc_s2[q])

    ww_tot_mean = vcat(ww_lc_s1, ww_lc_s2)
    ww_tot_syst = vcat(ww_lc_s1, ww_lc_s2)

    w, widx  =  findmax(ww_tot_mean)
  
    model_idx = mod(widx, length(f_tot_charm))
    if model_idx == 0 
        model_idx = length(f_tot_charm)
    end
    println("   wmax: ", w, " model_idx: ", model_idx)
    model = f_tot_charm[model_idx]
    
    
    cat_idx = Int((widx - model_idx ) / length(f_tot_charm))+1
    if cat_idx <= 0
        cat_idx +=1
    end
    println("   best χ2/χ2exp: ", fitcat_cc_mean[cat_idx].fit[model_idx].chi2 / fitcat_cc_mean[cat_idx].fit[model_idx].chi2exp)
    println("   widx: ", widx, " model_idx: ", model_idx, " catidx: ", cat_idx)
    ## Best Res
    best_mod = f_tot_charm[model_idx]
    xdata = fitcat_cc_mean[cat_idx].xdata
    param = fitcat_cc_mean[cat_idx].fit[model_idx].param

    ph_res_best = 4/9 * best_mod([0.0 phi2_ph phi4_ph], param)[1]; uwerr(ph_res_best)
    println("   best res: ", ph_res_best )
    ## histogram

    all_res_mean = Vector{uwreal}()
    for (k, cat) in enumerate(fitcat_cc_mean)
        for (j, mod) in enumerate(f_tot_charm)
            push!(all_res_mean, mod([0.0 phi2_ph phi4_ph], cat.fit[j].param)[1])
            # push!(all_res, cat.fit[j].param[1] )
        end
    end
    all_res_syst = Vector{uwreal}()
    for (k, cat) in enumerate(fitcat_cc_syst)
        for (j, mod) in enumerate(f_tot_charm)
            push!(all_res_syst, mod([0.0 phi2_ph phi4_ph], cat.fit[j].param)[1])
            # push!(all_res, cat.fit[j].param[1] )
        end
    end
    final_res, _ = 4 ./ 9 .* model_average(all_res_mean, ww_tot_mean); uwerr(final_res)
    _, syst = 4 ./ 9 .* model_average(all_res_syst, ww_tot_syst); uwerr(final_res)
    push!(RES, final_res)
    push!(SYST, syst)
    println("   Model ave:  ", final_res)
    println("   systematic: ", syst)


    hist(value.(all_res_mean) .* 4 ./ 9 , bins=80, histtype="stepfilled", alpha=0.5, ec="k", color="navy", weights=ww_tot_mean, zorder=3)
    fill_betweenx([0,0.6], value(final_res).+err(final_res), value(final_res).-err(final_res), alpha=0.4, color="gold", zorder=2)
    errtot = sqrt(err(final_res)^2 + syst^2)
    fill_betweenx([0,0.6], value(final_res).+errtot, value(final_res).-errtot, alpha=0.4, color="tomato", zorder=1)
    xlim(value(final_res)-10*err(final_res), value(final_res)+10*err(final_res))
    ylabel(L"$\mathrm{Frequency}$")
    xlabel(L"$\mathit{\bar\Pi}^{cc}_{\mathrm{sub}}(Q^2)$")
    tight_layout()
    display(gcf())
    savefig(joinpath(path_plot, "hist", "hist_q$(q).pdf"))
    close()

end

## check error with systematics
aux = similar(RES)
for k in eachindex(RES)
    aux[k] = RES[k] + uwreal([0.0, SYST[k]], "Syst Picc conn")
end
uwerr.(aux)
## saving physical results in BDIO

io = IOBuffer()
write(io, "PICC connected physical results")
fb = ALPHAdobs_create(joinpath(path_phys_res, "PIcc_conn_QLD_physRes.bdio"), io)
for k in eachindex(RES)
    aux = RES[k] + uwreal([0.0, SYST[k]], "Syst Picc  QLD")
    ALPHAdobs_write(fb, aux)
end
ALPHAdobs_close(fb)

## saving systematics in txt file
using DelimitedFiles
open(joinpath(path_phys_res, "systematics.txt"), "a") do io
    writedlm(io, ["# pi cc conn "])
    writedlm(io, [Qgev SYST])
end

## test reading
fb = BDIO_open(joinpath(path_phys_res, "PIcc_conn_physRes.bdio"), "r")
res = []
while ALPHAdobs_next_p(fb)
    d = ALPHAdobs_read_parameters(fb)
    push!(res, ALPHAdobs_read_next(fb))
end
BDIO_close!(fb)


## test cancelling fluctuations
details(pi_cc_ll_s1[1][10])
mchist(pi_cc_ll_s1[1][10], "sqrtt0 [fm] err")
set_fluc_to_zero!(pi_cc_ll_s1[1][10], "sqrtt0 [fm] err")
mchist(pi_cc_ll_s1[1][10], "sqrtt0 [fm] err")
pi_cc_ll_s1[1][10] *=1

uwerr(pi_cc_ll_s1[1][10])
details(pi_cc_ll_s1[1][10])

set_fluc_to_zero(pi_cc_ll_s1[1][10], "C101")
