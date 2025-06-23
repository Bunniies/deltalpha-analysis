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


path_data_pi = "/Users/alessandroconigli/MyDrive/postdoc-mainz/projects/deltalpha/PIdata/impr_deriv/new_strategy/Q_SD/"
path_kappa_target = "/Users/alessandroconigli/MyDrive/postdoc-mainz/projects/deltalpha/physical_results/kappa_c_target.bdio"

## reading kappa charm target
kappa_target = Dict()
fb = BDIO_open(path_kappa_target, "r")
while ALPHAdobs_next_p(fb)
    d = ALPHAdobs_read_parameters(fb)
    sz = tuple(d["size"]...)
    ks = collect(d["keys"])
    kappa_target = ALPHAdobs_read_next(fb, size=sz, keys=ks)
end
BDIO_close!(fb)

# reading HVP data for kappaC_sim and kappaC_sim_plus

#kappaC_sim
fb = BDIO_open(joinpath(path_data_pi, "PIcc_conn.bdio"), "r")
hvp_charm_sim = Dict()
while ALPHAdobs_next_p(fb)
    d = ALPHAdobs_read_parameters(fb)
    nobs = d["nobs"]
    dims = d["dimensions"]
    sz = tuple(d["size"]...)
    extra = d["extra"]
    ks = collect(d["keys"])
    hvp_charm_sim[extra["Ens"]] = ALPHAdobs_read_next(fb, size=sz, keys=ks)
end
BDIO_close!(fb)

#kappaC_sim_plus
fb = BDIO_open(joinpath(path_data_pi, "PIcc_conn_plus.bdio"), "r")
hvp_charm_sim_plus = Dict()
while ALPHAdobs_next_p(fb)
    d = ALPHAdobs_read_parameters(fb)
    nobs = d["nobs"]
    dims = d["dimensions"]
    sz = tuple(d["size"]...)
    extra = d["extra"]
    ks = collect(d["keys"])
    hvp_charm_sim_plus[extra["Ens"]] = ALPHAdobs_read_next(fb, size=sz, keys=ks)
end
BDIO_close!(fb)


enslist = collect(keys(kappa_target))


## interpolation
store_res = Dict()

for (k, ens) in enumerate(enslist) # loop over ens
    store_res[ens] = Dict{String, Array{uwreal}}()
    kappa_tar_val = kappa_target[ens][1]

    for kk in collect(keys(hvp_charm_sim[ens])) # loop over 4 data set
        data_c_sim = hvp_charm_sim[ens][kk]
        data_c_sim_plus = hvp_charm_sim_plus[ens][kk]

        tmp_res_store = []
        for qq in eachindex(data_c_sim) # loop over momenta
            kappa_fit = [kcd_in[ens]["kappaC_sim"], kcd_in[ens]["kappaC_sim_plus"]]
            data_fit = [data_c_sim[qq], data_c_sim_plus[qq]]

            fitp, chi2 = lin_fit(kappa_fit, data_fit)

            hvp_target = y_lin_fit(fitp, value(kappa_tar_val))

            push!(tmp_res_store, hvp_target)

            fig = figure(figsize=(10,7.))
            title(ens)
            uwerr.(data_fit)
            uwerr(kappa_tar_val)
            uwerr(hvp_target)
            errorbar(kappa_fit, value.(data_fit), yerr=err.(data_fit), fmt="s", color="#0072B2", capsize=2)
            errorbar(value(kappa_tar_val), xerr=err(kappa_tar_val), value(hvp_target), err(hvp_target), fmt="d", color="#D55E00", capsize=2)

            xlabel(L"$\kappa_c$")
            ylabel(L"$\mathit{\bar\Pi}_c$")
            display(fig)
            close("all")

        end
        store_res[ens][kk] = tmp_res_store 
    end
end

## save data
io = IOBuffer()
write(io, "PIcc lowq charm-conn interp ")

fb = ALPHAdobs_create(joinpath(path_data_pi, "PIcc_conn_interp.bdio"), io)
for (k,ens) in enumerate(enslist)
    extra = Dict{String,Any}("Ens" => ens)
    data = Dict{String, Array{uwreal}}(
        "picc_ll_s1" => store_res[ens]["picc_ll_s1"],
        "picc_lc_s1" => store_res[ens]["picc_lc_s1"],
        "picc_ll_s2" => store_res[ens]["picc_ll_s2"],
        "picc_lc_s2" => store_res[ens]["picc_lc_s2"],
    )
    ALPHAdobs_write(fb, data, extra=extra)
end
ALPHAdobs_close(fb)

## test reading

fb = BDIO_open(joinpath(path_data_pi, "PI_blc_2Qm_interp.bdio"), "r")
res = Dict()
while ALPHAdobs_next_p(fb)
    d = ALPHAdobs_read_parameters(fb)
    nobs = d["nobs"]
    dims = d["dimensions"]
    sz = tuple(d["size"]...)
    extra = d["extra"]
    ks = collect(d["keys"])
    res[extra["Ens"]] = ALPHAdobs_read_next(fb, size=sz, keys=ks)
end
BDIO_close!(fb)




