charge_factor = Dict(
    "88" => 1/3,
    "cc" => 4/9,
    "08" => 1/(6*sqrt(3))
)

function rational_func(; M=3, N=3)
    num = (x,p) -> sum(p[k] .* x.^k for k in 1:M)
    den = (x,p) -> 1. .+ sum(p[k+M] .* x.^k for k in 1:N)

    num1 = (x,p) -> sum(p[k+M+N] .* (x ./4).^k for k in 1:M)
    den1 = (x,p) -> 1. .+ sum(p[k+M+M+N] .* (x ./4).^k for k in 1:N)

    func(x,p) = num(x,p) ./ den(x,p) .- num1(x,p)./ den1(x,p)
    return func
end

function read_systematics(path)
    f = readdlm(joinpath(path, "systematics.txt"))
    delim = findall(x-> x=="#", f)
    dict_syst = Dict()
    for k in delim
        idx = k.I[1]
        name = join(f[idx,2:end])
        dict_syst[name] = Float64.(f[idx+1:idx+12, 2])
    end

    return dict_syst
end

function add_t0_err!(obs, t0phys)
    uwerr(t0phys)
    uwerr.(obs)
    for k in eachindex(obs)
        err_t0 = abs.(mchist(obs[k], "sqrtt0 [fm]") * 1e8 * err(t0phys))
        obs[k] = obs[k] + uwreal([0.0, err_t0[1]], "sqrtt0 [fm]")
    end
    return nothing
end

function add_t0_phi2_phi4_err!(obs, t0phys)
    uwerr(t0phys)
    uwerr.(obs)
    for k in eachindex(obs)
        err_t0_aux = mchist(obs[k], "sqrtt0 [fm]") * 1e8
        err_t0 = sqrt(err_t0_aux[1]^2 * err(t0phys)^2)
        obs[k] = obs[k] + uwreal([0.0, err_t0], "err t0")
    end
    return nothing
end