using Revise
using HVPobs
using ADerrors
import ADerrors: err

include("../utils/const.jl")


const path_data = "/Users/alessandroconigli/Lattice/data/HVP/2ptdata"
const path_rw   = "/Users/alessandroconigli/Lattice/data/HVP/rwf_deflated"


ensinfo = EnsInfo("J304")

rwf_r0 = read_ms1("/Users/alessandroconigli/Lattice/data/HVP/rwf_deflated/J304r000.ms1.dat", v="2.0")

rwf_r1 = read_ms1("/Users/alessandroconigli/Lattice/data/HVP/rwf_deflated/J304r001.ms1.dat", v="1.6")

rwf = get_rw(path_rw, "J304")# read_ms1(path_rw, "J304")
v1v1 = get_corr(path_data, ensinfo, "light", "V1V1", path_rw=path_rw)

corrConnected(path_data, ensinfo, "light", path_rw=path_rw)
corrDisconnected(path_data, ensinfo, "08", path_rw=path_rw)