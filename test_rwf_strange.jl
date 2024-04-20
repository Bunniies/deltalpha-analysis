path_rwf_strange = "/Users/alessandroconigli/MyDrive/postdoc-mainz/projects/deltalpha/resources/strange_sign.txt"
using DelimitedFiles
using Revise, HVPobs

function read_rwf_strange(path::String, id::String)
    f = readdlm(path)
    
    repLen = HVPobs.Data.CLS_CNFG[id]["repLen"]
    key = collect(keys(repLen))
    delim_ens = findall(x-> typeof(x)<:AbstractString && occursin(id, x) , f )

    strange_rwf = OrderedDict([k => ones(repLen[k]) for k in key])

    for (k,d) in enumerate(delim_ens)
        m = match(r"(?![0])\d{1,3}", split(f[d], "r")[2])
        m = isnothing(m) ? "r0" : "r"*m.match

        if m == key[k]
            idx_x = d.I[1]
            n_of_flagged_cnfg = f[idx_x+5,5][1]
            flagged_cnfg = filter(x-> typeof(x) == Int64, f[idx_x+4, 3:end])
            # println(flagged_cnfg)
            if n_of_flagged_cnfg != length(flagged_cnfg)
                error("Number of flagged configs does not match with the found flagged configs.")
            end
            strange_rwf[m][flagged_cnfg] .*= -1
        end
    end

    return strange_rwf
end

f = read_rwf_strange(path_rwf_strange, "H101")
delim_ens = findall(x-> typeof(x)<:AbstractString && occursin("H101", x) , f )
split(f[delim_ens[1]], "r")