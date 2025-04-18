using ADerrors
using ALPHAio, BDIO
import ADerrors:err

rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
rcParams["text.usetex"] =  true
rcParams["mathtext.fontset"]  = "cm"
rcParams["font.size"] =20
rcParams["axes.labelsize"] =26
rcParams["axes.titlesize"] = 22
plt.rc("text", usetex=true) # set to true if a LaTeX installation is present

include("../utils/IO_BDIO.jl")
include("./tools_running.jl")
include("../utils/const.jl")

path_phys_res_highq = "/Users/alessandroconigli/MyDrive/postdoc-mainz/projects/deltalpha/physical_results/scale_error_artificial/piQ_Q4/"


const Qgev = [0.05, 0.1, 0.4, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0] # Q^2
NMOM = length(Qgev)
Qmgev =  9.0  # GeV^2


b33qm_pt_5loops = uwreal([1896.11*1e-5, 2.67*1e-5], "b33_qm_pt") # from wolfram code
b33qqm_highq = [(q/(4*Qmgev)) * b33qm_pt_5loops for q in Qgev] # using PT result


deltalc_highq = read_phys_res(path_phys_res_highq, "PI_blc_2Qm_physRes.bdio")

bcc_highq = Qgev ./ (4*Qmgev) .* deltalc_highq .+ 2 .* charge_factor["cc"] .* b33qqm_highq 