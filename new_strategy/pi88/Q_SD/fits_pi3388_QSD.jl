using Revise, HVPobs
using PyPlot, LaTeXStrings
using BDIO, ADerrors, ALPHAio
import ADerrors: err
using Statistics


#================ SET UP PLOT PARAMETERS ==================#

rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
rcParams["text.usetex"] =  true
rcParams["mathtext.fontset"]  = "cm"
rcParams["font.size"] = 20
rcParams["axes.labelsize"] = 26
rcParams["axes.titlesize"] = 22
plt.rc("text", usetex=true) # set to true if a LaTeX installation is present


include("../../utils/const.jl")
include("../../utils/types.jl")
include("../../utils/plot_utils.jl")
include("../../utils/IO_BDIO.jl")
include("../../utils/tools.jl")
include("./func_comb_pi3388_QSD.jl")

path_bdio_obs = "/Users/alessandroconigli/MyDrive/postdoc-mainz/projects/deltalpha/data"
path_store_pi = "/Users/alessandroconigli/MyDrive/postdoc-mainz/projects/deltalpha/PIdata/impr_deriv/new_strategy/Q_SD/"
path_plot = "/Users/alessandroconigli/MyDrive/postdoc-mainz/projects/deltalpha/plots/new_strategy/pi88/Q_SD/"
path_phys_res = "/Users/alessandroconigli/MyDrive/postdoc-mainz/projects/deltalpha/physical_results/new_strategy/Q_SD/"

#======= PHYSICAL CONSTANTS ====================#
const MPI_ph = uwreal([134.9768, 0.0005], "mpi phys")
const MK_ph  = uwreal([495.011, 0.01], "mk phys")

const phi2_ph = (sqrt(8)*t0sqrt_ph * MPI_ph / hc)^2
const phi4_ph = (sqrt(8)*t0sqrt_ph)^2 * ((MK_ph/hc)^2 + 0.5*(MPI_ph/hc)^2)

# const Qgev = [3., 5., 9.] # Q^2
const Qgev = [4, 5, 6, 7, 8, 9, 12] # Q^2 # additional very high values

const Qmgev = 9.0 # Qm^2
##

enslist = sort([ "H102", "N101", "C101", "C102", "D150",
         "N452", "N451", "D450", "D451", "D452", 
         "N200", "N203", "D251", "D200", "D201", "E250", 
         "J303", "J304", "E300", 
         "J501"]
)
ensinfo = EnsInfo.(enslist)

#============== READ PS MASSES AND t0 FROM BDIO FILES =================#
dir_path = filter(x-> basename(x) in enslist, readdir(path_bdio_obs, join=true))
spectrum_path = vcat(filter(!isempty, [filter(x-> occursin("spectrum.bdio", x)  , readdir(dir_path[k], join=true)) for k in eachindex(dir_path)])...)

phi2 = Vector{uwreal}(undef, length(enslist))
phi4 = Vector{uwreal}(undef, length(enslist))
a28t0 = Vector{uwreal}(undef, length(enslist))
t0ens = Vector{uwreal}(undef, length(enslist))
MPI_TOT = []

for (k, ens) in enumerate(ensinfo)
    m_pi  = read_BDIO(spectrum_path[k], "spectrum", "mpi")[1]
    push!(MPI_TOT, m_pi)
    m_k   = read_BDIO(spectrum_path[k], "spectrum", "mk")[1]
    t0ens_aux = read_BDIO(spectrum_path[k], "spectrum", "t0")[1]

    phi2[k]  = 8 * t0ens_aux * m_pi^2
    phi4[k]  = 8 * t0ens_aux * (m_k^2 + 0.5*m_pi^2 )
    a28t0[k] =  a(ens.beta)^2 #1 / (8*t0ens_aux)
    t0ens[k] = t0ens_aux
end

#============= READ PI (33-88) FROM BDIO =================#
# set 1 improvement coefficients
pi_3388_ll_s1 = Vector{Vector{uwreal}}(undef, 0)
pi_3388_lc_s1 = Vector{Vector{uwreal}}(undef, 0)

# set 2 improvement coefficients
pi_3388_ll_s2 = Vector{Vector{uwreal}}(undef, 0)
pi_3388_lc_s2 = Vector{Vector{uwreal}}(undef, 0)

fb = BDIO_open(joinpath(path_store_pi, "PI_3388_Q_SD.bdio"), "r")
tmp_res = Dict()
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
    tmp_res[extra["Ens"]] = Dict{String, Array{uwreal}}(

        "pi3388_ll_s1" => res["pi3388_ll_s1"],
        "pi3388_lc_s1" => res["pi3388_lc_s1"],
        "pi3388_ll_s2" => res["pi3388_ll_s2"],
        "pi3388_lc_s2" => res["pi3388_lc_s2"]
    )

end
BDIO_close!(fb)

# rearrange order
for (k,ens) in enumerate(enslist)
    push!(pi_3388_ll_s1, tmp_res[ens]["pi3388_ll_s1"])
    push!(pi_3388_lc_s1, tmp_res[ens]["pi3388_lc_s1"])
    push!(pi_3388_ll_s2, tmp_res[ens]["pi3388_ll_s2"])
    push!(pi_3388_lc_s2, tmp_res[ens]["pi3388_lc_s2"])
end


##
#========= ADD FVE CORRECTIONS ==========#
fb = BDIO_open(joinpath(path_store_pi, "FVE_HP_Q_SD_non_sub_kernel.bdio"), "r")
fvc = Dict()
while ALPHAdobs_next_p(fb)
    d = ALPHAdobs_read_parameters(fb)
    nobs = d["nobs"]
    dims = d["dimensions"]
    sz = tuple(d["size"]...)
    extra = d["extra"]
    ks = collect(d["keys"])
    fvc[extra["Ens"]] = ALPHAdobs_read_next(fb, size=sz, keys=ks)
end
BDIO_close!(fb)
@info("Adding Finite-Volume Corrections")
for (k,ens) in enumerate(ensinfo)
    corr_pi = fvc[ens.id]["fvc_pi"] 
    corr_k = fvc[ens.id]["fvc_k"] 

    pi_3388_ll_s1[k] .+= (corr_pi .+ corr_k .- 2/3*corr_k)
    pi_3388_ll_s2[k] .+= (corr_pi .+ corr_k .- 2/3*corr_k)
    pi_3388_lc_s1[k] .+= (corr_pi .+ corr_k .- 2/3*corr_k)
    pi_3388_lc_s2[k] .+= (corr_pi .+ corr_k .- 2/3*corr_k)
end


##############################
## CREATE FIT CATEGORIES
##############################
NMOM = length(Qgev)
# set 1 picc
fitcat_3388_ll_s1 = [Vector{FitCat}(undef,0) for i=1:NMOM]
fitcat_3388_lc_s1 = [Vector{FitCat}(undef,0) for i=1:NMOM]
# set 2 picc
fitcat_3388_ll_s2 = [Vector{FitCat}(undef,0) for i=1:NMOM]
fitcat_3388_lc_s2 = [Vector{FitCat}(undef,0) for i=1:NMOM]

## no cuts in data
for s in 1:2
    xdata = [a28t0 phi2 phi4]
    if s == 1
        for q in 1:NMOM
            str = "all_data_set$(s)_q$(q)"
            # pi 33 - 88
            push!(fitcat_3388_ll_s1[q], FitCat(xdata, getindex.(pi_3388_ll_s1, q), str))
            push!(fitcat_3388_lc_s1[q], FitCat(xdata, getindex.(pi_3388_lc_s1, q), str))
        end
    elseif s == 2
        for q in 1:NMOM
            str = "all_data_set$(s)_q$(q)"
            # pi 33 - 88
            push!(fitcat_3388_ll_s2[q], FitCat(xdata, getindex.(pi_3388_ll_s2, q), str))
            push!(fitcat_3388_lc_s2[q], FitCat(xdata, getindex.(pi_3388_lc_s2, q), str))
        end
    end
end

## cuts in phi2<0.45
i_cutphi2 = findall(x->x<0.45, value.(phi2))
for s in 1:2
    xdata = [a28t0[i_cutphi2] phi2[i_cutphi2] phi4[i_cutphi2]]
    if s == 1
        for q in 1:NMOM
            str = "mpicut_set$(s)_q$(q)"
            # pi 33 - 88
            push!(fitcat_3388_ll_s1[q], FitCat(xdata, getindex.(pi_3388_ll_s1, q)[i_cutphi2], str))
            push!(fitcat_3388_lc_s1[q], FitCat(xdata, getindex.(pi_3388_lc_s1, q)[i_cutphi2], str))
        end
    elseif s == 2
        for q in 1:NMOM
            str = "mpicut_set$(s)_q$(q)"
            # pi 33 - 88
            push!(fitcat_3388_ll_s2[q], FitCat(xdata, getindex.(pi_3388_ll_s2, q)[i_cutphi2], str))
            push!(fitcat_3388_lc_s2[q], FitCat(xdata, getindex.(pi_3388_lc_s2, q)[i_cutphi2], str))
        end
    end
end


##
#================= FITTING ====================#
for q in 1:NMOM

    for k_cat in eachindex(fitcat_3388_ll_s1[q])
        xdata = fitcat_3388_ll_s1[q][k_cat].xdata
        ydata_ll_s1 = fitcat_3388_ll_s1[q][k_cat].ydata
        ydata_ll_s2 = fitcat_3388_ll_s2[q][k_cat].ydata
        ydata_lc_s1 = fitcat_3388_lc_s1[q][k_cat].ydata
        ydata_lc_s2 = fitcat_3388_lc_s2[q][k_cat].ydata

        for (k_mod, model) in enumerate(f_tot_dltiso)
            println(q, " ", label_tot_dltiso[k_mod])
        
            println(k_mod)
            fit_ll_s1 = fit_routine(model, value.(xdata), ydata_ll_s1, n_par_tot_dltiso[k_mod], pval=false)
            fit_ll_s2 = fit_routine(model, value.(xdata), ydata_ll_s2, n_par_tot_dltiso[k_mod], pval=false)
            fit_lc_s1 = fit_routine(model, value.(xdata), ydata_lc_s1, n_par_tot_dltiso[k_mod], pval=false)
            fit_lc_s2 = fit_routine(model, value.(xdata), ydata_lc_s2, n_par_tot_dltiso[k_mod], pval=false)
            push!(fitcat_3388_ll_s1[q][k_cat].fit, fit_ll_s1)
            push!(fitcat_3388_ll_s2[q][k_cat].fit, fit_ll_s2)
            push!(fitcat_3388_lc_s1[q][k_cat].fit, fit_lc_s1)
            push!(fitcat_3388_lc_s2[q][k_cat].fit, fit_lc_s2)
        end
    end
end

##########################
## PLOTS
#########################

ll = L"$-\widehat{\Delta}_{ls}(Q^2)$"

plot_cl_all_set(fitcat_3388_ll_s1, fitcat_3388_ll_s2, fitcat_3388_lc_s1, fitcat_3388_lc_s2, nmom=NMOM, path_plot=path_plot, ylab=ll, f_tot_isov=f_tot_dltiso)
plot_chiral_best_fit(fitcat_3388_ll_s2, path_plot=path_plot, tt=["Set", "2", "LL"], nfit=0, f_tot_isov=f_tot_dltiso, nmom=NMOM, ylab=ll) # nfit=28 used for 2024 proceedings
plot_cl_best_fit(fitcat_pi33_ll_s2, nmom=NMOM, path_plot=nothing, tt=["Set", "2", "LL"], f_tot_isov=f_tot_isov, ylab=ll)

cattot = [vcat(fitcat_3388_ll_s2[k],fitcat_3388_lc_s2[k]...) for k in eachindex(fitcat_3388_lc_s1)]
plot_mAve_summary(cattot, nmom=NMOM, xlab=vcat(label_tot_dltiso,label_tot_dltiso,label_tot_dltiso,label_tot_dltiso), charge_factor=1., ylab=ll, tt=["Set", "2", "LL-LC"], path_plot=path_plot)
plot_mAve_summary(fitcat_pi33_ll_s2, nmom=NMOM, xlab=vcat(label_tot_isov,label_tot_isov,label_tot_isov), charge_factor=1., ylab=ll)


###################################
## RESULTS
##################################

RES = []
RES_all = []
SYST = []
for q in 1:NMOM
    @info "Momentum no. $(q): $(Qgev[q]) GeV^2"
    fitcat_pi3388_mean = vcat(vcat(fitcat_3388_ll_s2[q],
                fitcat_3388_lc_s2[q])...
    )
    fitcat_pi3388_syst = vcat(vcat(fitcat_3388_ll_s1[q],
    fitcat_3388_ll_s2[q],
    fitcat_3388_lc_s1[q],
    fitcat_3388_lc_s2[q])...
    )
    # ww_tot = get_w_from_fitcat(fitcat_pi3388_tot)

    ww_ll_s1 = get_w_from_fitcat(fitcat_3388_ll_s1[q])
    ww_ll_s2 = get_w_from_fitcat(fitcat_3388_ll_s2[q])
    ww_lc_s1 = get_w_from_fitcat(fitcat_3388_lc_s1[q])
    ww_lc_s2 = get_w_from_fitcat(fitcat_3388_lc_s2[q])

    ww_tot_mean = vcat(ww_ll_s2, ww_lc_s2)
    ww_tot_syst = vcat(ww_ll_s1, ww_ll_s2, ww_lc_s1, ww_lc_s2)

    w, widx  =  findmax(ww_tot_mean)
  
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
    println("   best χ2/χ2exp: ", fitcat_pi3388_mean[cat_idx].fit[model_idx].chi2 / fitcat_pi3388_mean[cat_idx].fit[model_idx].chi2exp)
    println("   widx: ", widx, " model_idx: ", model_idx, " catidx: ", cat_idx)
    
    ## Best Res
    best_mod = f_tot_dltiso[model_idx]
    xdata = fitcat_pi3388_mean[cat_idx].xdata
    param = fitcat_pi3388_mean[cat_idx].fit[model_idx].param

    ph_res_best =  1 / 3 * best_mod([0.0 phi2_ph phi4_ph], param)[1]; uwerr(ph_res_best)
    println("   best res: ", ph_res_best )
    ## histogram

    all_res_mean = Vector{uwreal}()
    for (k, cat) in enumerate(fitcat_pi3388_mean)
        for (j, mod) in enumerate(f_tot_dltiso)
            push!(all_res_mean, mod([0.0 phi2_ph phi4_ph], cat.fit[j].param)[1])
            #push!(all_res, mod([0.0 value(phi2_ph) value(phi4_ph)], cat.fit[j].param)[1])
        end
    end
    all_res_syst = Vector{uwreal}()
    for (k, cat) in enumerate(fitcat_pi3388_syst)
        for (j, mod) in enumerate(f_tot_dltiso)
            push!(all_res_syst, mod([0.0 phi2_ph phi4_ph], cat.fit[j].param)[1])
            #push!(all_res, mod([0.0 value(phi2_ph) value(phi4_ph)], cat.fit[j].param)[1])
        end
    end

    final_res, _ =  1 ./ 3 .* model_average(all_res_mean, ww_tot_mean); uwerr(final_res)
    final_res_all, syst =  1 ./ 3 .* model_average(all_res_syst, ww_tot_syst); uwerr(final_res)

    push!(RES, final_res)
    push!(RES_all, final_res_all)
    push!(SYST, syst)
    println("   Model ave:  ", final_res)
    println("   systematic: ", syst)
    println("\n")


    hist(value.(all_res_mean) ./ 3, bins=40, histtype="stepfilled", alpha=0.5, ec="k", color="navy", weights=ww_tot_mean, zorder=3)
    fill_betweenx([0,0.6], value(final_res).+err(final_res), value(final_res).-err(final_res), alpha=0.4, color="gold", zorder=2)
    errtot = sqrt(err(final_res)^2 + syst^2)
    fill_betweenx([0,0.6], value(final_res).+errtot, value(final_res).-errtot, alpha=0.4, color="tomato", zorder=1)
    xlim(value(final_res)-6*err(final_res), value(final_res)+6*err(final_res))
    ylabel(L"$\mathrm{Frequency}$")
    xlabel(L"$-\Delta_{ls}(\Delta\alpha)$")
    tight_layout()
    display(gcf())
    savefig(joinpath(path_plot, "hist", "hist_q$(q).pdf"))
    close()

end

## saving physical results in BDIO
io = IOBuffer()
write(io, "PI3388 SD physical results")
fb = ALPHAdobs_create(joinpath(path_phys_res, "PI3388_QSD_physRes.bdio"), io)
for k in eachindex(RES)
    aux = RES[k] + uwreal([0.0, SYST[k]], "Syst Pi3388 SD")
    ALPHAdobs_write(fb, aux)
end
ALPHAdobs_close(fb)

## test reading
fb = BDIO_open(joinpath(path_phys_res, "PI3388_QSD_physRes.bdio"), "r")
res = []
while ALPHAdobs_next_p(fb)
    d = ALPHAdobs_read_parameters(fb)
    push!(res, ALPHAdobs_read_next(fb))
end
BDIO_close!(fb)
