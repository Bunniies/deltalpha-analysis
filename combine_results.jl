using Revise, HVPobs
using PyPlot, LaTeXStrings
using BDIO, ADerrors, ALPHAio
using QuadGK
import ADerrors: err

path_phys_res = "/Users/alessandroconigli/MyDrive/postdoc-mainz/projects/deltalpha/physical_results/scale_error_multimom/"


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


## reading charm connected

fb = BDIO_open(joinpath(path_phys_res, "PIcc_conn_physRes.bdio"), "r")
res_cc_conn = []
while ALPHAdobs_next_p(fb)
    d = ALPHAdobs_read_parameters(fb)
    push!(res_cc_conn, ALPHAdobs_read_next(fb))
end
BDIO_close!(fb)
uwerr.(res_cc_conn)

## reading systematics
using DelimitedFiles
fname = readdlm(joinpath(path_phys_res, "systematics.txt"))
syst_SD = fname[2:13, 2 ]
syst_ILD = fname[15:26, 2 ]
syst_cc_conn = fname[28:end, 2 ]


##
res_33_SD = [res_33_SD[k] + + uwreal([0.0, syst_SD[k]], "syst 33 SD") for k in eachindex(res_33_SD)]
res_33_ILD = [res_33_ILD[k] + + uwreal([0.0, syst_ILD[k]], "syst 33 ILD") for k in eachindex(res_33_ILD)]
sum33 = res_33_SD .+ res_33_ILD ; uwerr.(sum33)
sum33_syst = [sum33[k] + uwreal([0.0, syst_SD[k]], "syst 33 SD") +  uwreal([0.0, syst_ILD[k]], "syst 33 ILD") for k in eachindex(sum33)]
uwerr.(sum33_syst)


res_cc_conn_syst = [res_cc_conn[k] + + uwreal([0.0, syst_cc_conn[k]], "syst cc conn") for k in eachindex(res_cc_conn)]
uwerr.(res_cc_conn_syst)

