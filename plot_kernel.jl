##################################################################
# This file was created by Alessandro Conigli - 03/2024
# Here we compute all the correlator for a given ensemble and we then 
# compute the kernel for plotting all the contributions together
##################################################################

using Revise
using HVPobs
using PyPlot, LaTeXStrings
using OrderedCollections
using TimerOutputs
using BDIO, ADerrors
import ADerrors: err
using ALPHAio

include("utils/const.jl")
include("utils/types.jl")
include("utils/plot_utils.jl")
include("utils/IO_BDIO.jl")
include("utils/tools.jl")

#================ SET UP PLOT PARAMETERS ==================#

rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
rcParams["text.usetex"] =  true
rcParams["mathtext.fontset"]  = "cm"
rcParams["font.size"] =13
rcParams["axes.labelsize"] =22
rcParams["axes.titlesize"] = 18
plt.rc("text", usetex=true) # set to true if a LaTeX installation is present



const path_data = "/Users/alessandroconigli/Lattice/data/HVP/2ptdata"
const path_rw   = "/Users/alessandroconigli/Lattice/data/HVP/rwf_deflated"
const path_store_bdio = "/Users/alessandroconigli/MyDrive/postdoc-mainz/projects/deltalpha/corr/impr_deriv"
const path_plot = "/Users/alessandroconigli/MyDrive/postdoc-mainz/projects/deltalpha/plots/ensembles"

const IMPR = true
const IMPR_SET = "1"
const STD_DERIV = false
const RENORM = true

const KRNL = krnl_dα_qhalf      # standard kernel 
Qgev = [3., 5., 9.] # Q^2
Qmgev = 9.0 # Qm^2

ensinfo = EnsInfo("E250")

gll_ll, gll_lc = corrConnected(path_data, ensinfo, "light", path_rw=path_rw, impr=IMPR, impr_set=IMPR_SET, std=STD_DERIV)
gss_ll, gss_lc = corrConnected(path_data, ensinfo, "strange", path_rw=path_rw, impr=IMPR, impr_set=IMPR_SET, std=STD_DERIV)

g88_ll_conn = Corr(1/6 .*( gll_ll.obs + 2*gss_ll.obs), ensinfo.id, "G88_ll_conn" )
g88_lc_conn = Corr(1/6 .*( gll_lc.obs + 2*gss_lc.obs), ensinfo.id, "G88_lc_conn" )

gcc_ll_conn, gcc_lc_conn = corrConnected(path_data, ensinfo, "charm_plus", path_rw=path_rw, impr=IMPR, impr_set=IMPR_SET, std=STD_DERIV) 

g88_ll_disc, g88_lc_disc, _ = corrDisconnected(path_data, ensinfo, "88", path_rw=path_rw, impr=IMPR, impr_set=IMPR_SET, std=STD_DERIV)
gcc_ll_disc, gcc_lc_disc, gcc_cc_disc = corrDisconnected(path_data, ensinfo, "cc", path_rw=path_rw, impr=false, impr_set=IMPR_SET, std=STD_DERIV)
gc8_ll_disc, gc8_lc_disc, gc8_cc_disc = corrDisconnected(path_data, ensinfo, "c8", path_rw=path_rw, impr=false, impr_set=IMPR_SET, std=STD_DERIV)

g08_ll_disc, g08_lc_disc, _ = corrDisconnected(path_data, ensinfo, "08", path_rw=path_rw, impr=IMPR, impr_set=IMPR_SET, std=STD_DERIV)

g08_lc = Corr((1/2sqrt(3) .* (gll_lc.obs .- gss_lc.obs) .+ g08_lc_disc.obs ), ensinfo.id, "G08_lc")

if RENORM
    Z3 = get_Z3(ensinfo, impr_set=IMPR_SET) 
    renormalize!(gll_ll, Z3^2)  
    renormalize!(gll_lc, Z3)
    
    Z8 = get_Z8(ensinfo, impr_set=IMPR_SET)
    renormalize!(g88_ll_conn, Z8^2)
    renormalize!(g88_lc_conn, Z8)
    renormalize!(g88_ll_disc, Z8^2)
    renormalize!(g88_lc_disc, Z8)
    renormalize!(g08_lc, Z8)

    Zcc = Zvc_l[ensinfo.id]
    renormalize!(gcc_ll_conn, Zcc*Zcc)
    renormalize!(gcc_lc_conn, Zcc)
end

gdelta_iso_ll = Corr(-g88_ll_conn.obs - g88_ll_disc.obs + 0.5*gll_ll.obs, gll_ll.id, "deltaiso_ll")
gdelta_iso_lc = Corr(-g88_lc_conn.obs - g88_lc_disc.obs + 0.5*gll_lc.obs, gll_lc.id, "deltaiso_lc")


## COMPUTE TMR

Qlat = Qgev .* value.(t0sqrt_ph.^2) ./ t0(ensinfo.beta) ./hc^2 *1e6
qmlat = Qmgev * value(t0sqrt_ph^2) / t0(ensinfo.beta) /hc^2 *1e6

for (j,q) in enumerate(Qlat)
    if j in [1,3]
        continue
    end
    _, data33 = tmr_integrand(gll_lc, q, KRNL, pl=false, t0ens=t0(ensinfo.beta), data=true)
    _, datadltiso = tmr_integrand(gdelta_iso_lc, q, KRNL, pl=false, t0ens=t0(ensinfo.beta), data=true)
    _, datacc_conn = tmr_integrand(gcc_lc_conn, q, KRNL, pl=false, t0ens=t0(ensinfo.beta), data=true)
    _, datacc_disc = tmr_integrand(gcc_cc_disc, q, KRNL, pl=false, t0ens=t0(ensinfo.beta), data=true, wind=Window("SD"))
    _, datac8_disc = tmr_integrand(gc8_cc_disc, q, KRNL, pl=false, t0ens=t0(ensinfo.beta), data=true)
    _, data08 = tmr_integrand(g08_lc, q, KRNL, pl=false, t0ens=t0(ensinfo.beta), data=true )

    data33 .*= 0.5 ; uwerr.(data33)
    datadltiso .*= 10; uwerr.(datadltiso)
    datacc_conn .*= 4/9; uwerr.(datacc_conn)
    datacc_disc .*= -1000*4/9; uwerr.(datacc_disc)
    datac8_disc .*= -1000 * 2 / (3*sqrt(3)); uwerr.(datac8_disc)
    data08 .*= 10 * 1/3 ; uwerr.(data08)

    xx = collect(0:length(data33)-1) .* value(a(ensinfo.beta))
    errorbar(xx, value.(data33), err.(data33), fmt="v", color="tomato", capsize=4, label=L"$G^{(3,3)}$")
    errorbar(xx, value.(datadltiso), err.(datadltiso), fmt="^", color="forestgreen", capsize=4, label=L"$10\cdot(G^{(3,3)} - G^{(8,8)})$")
    errorbar(xx, value.(datacc_conn), err.(datacc_conn), fmt="d", color="royalblue", capsize=4, label=L"$\frac{4}{9}G^{(c,c)}_{\mathrm{conn}}$")
    errorbar(xx, value.(data08), err.(data08), fmt="o", color="black", capsize=4, label=L"$10\cdot\frac{1}{3}G^{(0,8)}$")
    errorbar(xx, value.(datacc_disc), err.(datacc_disc), fmt="d", color="gold", capsize=4, label=L"$1000\cdot\frac{4}{9}G^{(c,c)}_{\mathrm{disc}}$")
    errorbar(xx, value.(datac8_disc), err.(datac8_disc), fmt="o", color="plum", capsize=4, label=L"$1000\cdot\frac{2}{3\sqrt{3}}G^{(c,8)}_{\mathrm{disc}}$")


    xlabel(L"$t \ [\mathrm{fm}]$")
    ylabel(L"$G^{(d,e)}(t) \cdot K(t, Q^2)$")
    xlim(-0.05,2.0) 
    legend(ncol=1)
    tight_layout()
    display(gcf())
    savefig(joinpath(path_plot, ensinfo.id, "$(ensinfo.id)_dalpha_kernel_SD_vs_ILD.pdf"))
    close()
end