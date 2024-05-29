using HVPobs, ADerrors

const path_data = "/Users/alessandroconigli/Lattice/data/HVP/2ptdata"
const path_rw   = "/Users/alessandroconigli/Lattice/data/HVP/rwf_deflated"

Gamma_l = ["V1V1", "V2V2", "V3V3", "V1T10", "V2T20", "V3T30"]
Gamma_c = ["V1V1c", "V2V2c", "V3V3c", "V1cT10", "V2cT20", "V3cT30"]
ensinfo = EnsInfo("H101")
sector = "light"

v1v1 = get_corr(path_data, ensinfo, sector, Gamma_l[1], path_rw=path_rw, frw_bcwd=false, L=1)
v2v2 = get_corr(path_data, ensinfo, sector, Gamma_l[2], path_rw=path_rw, frw_bcwd=false, L=1)
v3v3 = get_corr(path_data, ensinfo, sector, Gamma_l[3], path_rw=path_rw, frw_bcwd=false, L=1)

v1t10 = get_corr(path_data, ensinfo, sector, Gamma_l[4], path_rw=path_rw, frw_bcwd=false, L=1)
v2t20 = get_corr(path_data, ensinfo, sector, Gamma_l[5], path_rw=path_rw, frw_bcwd=false, L=1)
v3t30 = get_corr(path_data, ensinfo, sector, Gamma_l[6], path_rw=path_rw, frw_bcwd=false, L=1)

function compute_derivative(c::Corr)
    tvals = vcat(collect(0:length(c.obs)/2), -reverse(collect(1:length(c.obs)/2-1))...)
    # return tvals
    corr_times_t = c.obs .* tvals.^2
    corr_times_t_deriv = (corr_times_t[3:end] .- corr_times_t[1:end-2]) / 2
    corr_times_t_deriv[1] = corr_times_t[3] - corr_times_t[2]
    push!(corr_times_t_deriv, corr_times_t[end] - corr_times_t[end-1])
    println("length(t):      ", length(tvals))
    println("length(derive): ", length(corr_times_t_deriv))
    
    dcorr = (1 ./ tvals[2:end].^2) .* (corr_times_t_deriv .- 2 .* tvals[2:end] .* c.obs[2:end])

    return dcorr
end

tvals = compute_derivative(v1t10)

HVPobs.Obs.improve_derivative(v1t10.obs, std=false)