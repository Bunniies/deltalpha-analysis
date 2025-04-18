using HVPobs
using ALPHAio, BDIO
using PyPlot, LaTeXStrings
using ADerrors
import ADerrors:err
using OrderedCollections
using Interpolations
using QuadGK
using DelimitedFiles
#================ SET UP PLOT PARAMETERS ==================#

rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
rcParams["text.usetex"] =  true
rcParams["mathtext.fontset"]  = "cm"
rcParams["font.size"] =20
rcParams["axes.labelsize"] =26
rcParams["axes.titlesize"] = 22
plt.rc("text", usetex=true) # set to true if a LaTeX installation is present

path_phys_res_highq = "/Users/alessandroconigli/MyDrive/postdoc-mainz/projects/deltalpha/physical_results/scale_error_artificial/piQ_Q4/"
path_phys_res_lowq = "/Users/alessandroconigli/MyDrive/postdoc-mainz/projects/deltalpha/physical_results/scale_error_artificial/piQ4_piQ0/"
path_phys_res_fullq = "/Users/alessandroconigli/MyDrive/postdoc-mainz/projects/deltalpha/physical_results/scale_error_artificial/full_q/"
path_plot = "/Users/alessandroconigli/MyDrive/postdoc-mainz/projects/deltalpha/plots/running"

path_adpy_data = "/Users/alessandroconigli/MyDrive/postdoc-mainz/projects/deltalpha/running_jegerlhener/running_AdlerPy/charm_adler"

include("../utils/IO_BDIO.jl")
include("./tools_running.jl")
include("../utils/const.jl")

alpha0=1/137.035999180;
Mz = 91.1876

function read_adlerpy_datfile(path)
    fname = readdlm(path, comments=true)

    qval = fname[:,2]
    d_val = fname[:,3]
    d_err = fname[:,4]

    ad_res = []
    for k in eachindex(qval)
        push!(ad_res, uwreal([d_val[k], d_err[k]], "Dfunc Adpy"))
    end
    return qval, ad_res
end
function get_interp(qval, adler)
    itp = interpolate((qval,), (2) / (12*pi^2)  * adler , Gridded(Linear()))
    return itp
end

function compute_integral(itp, a, b)
    
    int_range = collect(range(a,b,length=10000))
    int_res = sum(itp(int_range)[1:end] .* push!(diff(int_range), 0.0) ./ int_range[1:end])
    return int_res 
end

const Qgev = [0.05, 0.1, 0.4, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0] # Q^2
NMOM = length(Qgev)
Qmgev =  9.0  # GeV^2


b33qm_pt_5loops = uwreal([1899.07*1e-5, 1.88*1e-5], "b33_qm_pt")
b33qqm_highq = [(q/(4*Qmgev)) * b33qm_pt_5loops for q in Qgev] # using PT result

# cc connected 
picc_highq_sub = read_phys_res(path_phys_res_highq, "PIcc_conn_physRes.bdio")
bcc_highq = read_phys_res(path_phys_res_highq, "PI_blc_2Qm_physRes.bdio")

picc_highq = picc_highq_sub .+ Qgev ./ (4*Qmgev) .* bcc_highq .+ 2 * charge_factor["cc"] * b33qqm_highq
#picc_highq ./= charge_factor["cc"]

##
qval, dfunc_adpy_msbar = read_adlerpy_datfile(joinpath(path_adpy_data, "adler_charm_msbar.dat"))
qval, dfunc_adpy_mpole = read_adlerpy_datfile(joinpath(path_adpy_data, "adler_charm_mpole.dat"))

itp_adpy_msbar = get_interp(qval, dfunc_adpy_msbar)
itp_adpy_mpole = get_interp(qval, dfunc_adpy_mpole)

store_res_adpy_msbar = []
store_res_adpy_mpole = []

for mom in collect(2:9)
    # AdPy MSBAR
    res_aux = compute_integral(itp_adpy_msbar, sqrt(mom/4), sqrt(mom))
    push!(store_res_adpy_msbar, res_aux)
    # AdPy mpole
    res_aux = compute_integral(itp_adpy_mpole, sqrt(mom/4), sqrt(mom))
    push!(store_res_adpy_mpole, res_aux)
end

store_res_adpy_mpole = store_res_adpy_mpole #* 1 / (12 *pi^2)
store_res_adpy_msbar = store_res_adpy_msbar #* 1 / (12 *pi^2)

uwerr.(store_res_adpy_msbar)
uwerr.(store_res_adpy_mpole)

##
uwerr.(picc_highq)
fig = figure(figsize=(10,8))

errorbar(Qgev, value.(picc_highq), err.(picc_highq),ms=10,mfc="none", fmt="d", label="Lattice", color="navy")

errorbar(collect(2:9).+0.1, value.(store_res_adpy_msbar), err.(store_res_adpy_msbar), fmt="^", mfc="none", color="tomato", ms=10, label=L"$\mathcal{O}(\alpha_s^3) \ \overline{\mathrm{MS}}$")
plot(collect(2:9).+0.1, value.(store_res_adpy_msbar), ls="dashed",color="tomato")

errorbar(collect(2:9), value.(store_res_adpy_mpole), err.(store_res_adpy_mpole), fmt="s", mfc="none", color="forestgreen", ms=10, label=L"$\mathcal{O}(\alpha_s^3)\ \mathrm{OS}$")
plot(collect(2:9), value.(store_res_adpy_mpole), ls="dashed", color="forestgreen")

xlim(1.5,9.5)
xlabel(L"$Q^2 \ [\mathrm{GeV}^2]$")
ylabel(L"$\widehat{\mathit{\Pi}}^{(c,c)}(Q^2)$")
legend()
tight_layout()
savefig(joinpath(path_plot, "charm_lattice_vs_pt.pdf"))
display(fig)