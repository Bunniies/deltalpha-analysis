using HVPobs, PyPlot, ADerrors

pp  = "/Users/alessandroconigli/Lattice/data/HVP/2ptdata/disc/C101/disc/2pt/88.npz"

res = read_disconnected_from_npz(pp, "H102")
vvdata = res["VV"].re_data
Tfalse = size(vvdata,2) 
Tt = Int64(96 / 2)
idx = vcat(collect(1:Tt), collect(Tfalse-Tt+1:Tfalse)...)
vvtest = vvdata[:, idx]
mean(vvtest, dims=1)#[45:51]

corr_vv = corr_obs(res["VV"])
##
yy = corr_vv.obs; uwerr.(yy)
xx = collect(1:length(yy))
errorbar(xx, value.(yy), err.(yy), fmt="s")
display(gcf())
close("all")


## OLD PRELIMINARY TESTS

using PyCall

np = pyimport("numpy")
pp  = "/Users/alessandroconigli/Lattice/data/HVP/2ptdata/disc/E250/disc/2pt/88.npz"
data = np.load(pp, allow_pickle=true)
kkk = collect(data.keys())


dataraw = get(data, "VcVc")
corrv1v1 = [] 
corrv2v2 = [] 
corrv3v3 = [] 
corrvv_av_first = []

using Statistics
vv_av  =dropdims(mean(dataraw[:,2:4,:], dims=2), dims=2)

for t in 1:192
    push!(corrv1v1, uwreal(real(dataraw[:,2,t]), "E2501"))
    push!(corrv2v2, uwreal(real(dataraw[:,3,t]), "E2501"))
    push!(corrv3v3, uwreal(real(dataraw[:,4,t]), "E2501"))
    push!(corrvv_av_first, uwreal(real(vv_av[:,t]), "E2501"))
end

corrvv_av_after = (corrv1v1 + corrv2v2 + corrv3v3)/3

cdatavv[:,2:4,:]



ncfg, nsrc= size(cdatavv)
T = size(cdatavv[1,2])[1]
icfg = collect(get(data, "icfg"))
rep =  String.((getindex.(split.(icfg, "n"),1)))
rep_len = [count(x->x==i, rep) for i in unique(rep) ]
idm_aux = parse.(Int64, getindex.(split.(icfg, "n") ,2))

re_data = Array{Float64}(undef, ncfg, T, nsrc-1)
im_data = similar(re_data)
for n in 1:ncfg
    for s in 1:nsrc-1
        re_data[n, :, s  ] = real(collect(cdatavv[n,s+1]))
        # im_data[n, :, s  ] = im(collect(cdatavv[n,s+1]))
    end
end
using Statistics
dropdims(mean(re_data, dims=3), dims=3)

function read_disconnected_from_npz(path::String, id::String)

    np = pyimport("numpy")
    data_raw = np.load(path, allow_pickle=true)
    
    icfg = collect(get(data, "icfg"))
    rep =  String.(getindex.(split.(icfg, "n"),1))
    rep_len = [count(x->x==i, rep) for i in unique(rep) ]
    idm_aux = parse.(Int64, getindex.(split.(icfg, "n") ,2))

    
    KEYS = ["VV", "VVc", "VT", "VcT"] 
    dict_res = Dict()

    for kk in KEYS
        data = get(data_raw, kk)

        ncfg, nsrc = size(data)
        T = size(data[1,2])[1]
        re_data = Array{Float64}(undef, ncfg, T, nsrc-1)

        for n in 1:ncfg
            for s in 1:nsrc-1
                re_data[n, :, s  ] = real(collect(cdatavv[n,s+1]))
            end
        end

        re_data = dropdims(mean(re_data, dims=3), dims=3)

        dict_res[kk] = CData(id, rep_len, re_data, re_data, idm_aux, kk)
    end
    return dict_res
end