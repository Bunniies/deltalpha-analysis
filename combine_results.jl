using Revise, HVPobs
using PyPlot, LaTeXStrings
using BDIO, ADerrors, ALPHAio
using QuadGK
import ADerrors: err

path_phys_res = "/Users/alessandroconigli/MyDrive/postdoc-mainz/projects/deltalpha/physical_results"


## reading isovector SD
fb = BDIO_open(joinpath(path_phys_res, "PI33_SD_physRes.bdio"), "r")
res_33_SD = []
while ALPHAdobs_next_p(fb)
    d = ALPHAdobs_read_parameters(fb)
    push!(res_33_SD, ALPHAdobs_read_next(fb))
end
BDIO_close!(fb)
uwerr.(res_33_SD)

## reading isovector ILD

fb = BDIO_open(joinpath(path_phys_res, "PI33_ILD_physRes.bdio"), "r")
res_33_ILD = []
while ALPHAdobs_next_p(fb)
    d = ALPHAdobs_read_parameters(fb)
    push!(res_33_ILD, ALPHAdobs_read_next(fb))
end
BDIO_close!(fb)
uwerr.(res_33_ILD)

##
sum33 = res_33_SD .+ res_33_ILD; uwerr.(sum33)

details(sum33[3])
