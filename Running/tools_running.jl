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