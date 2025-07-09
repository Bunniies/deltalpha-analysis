# here we read the various contributions computed in SD ID and LD and combine them
# together with the inclusion of the subtracted functions for isovector and charm

using Revise, HVPobs
using PyPlot, LaTeXStrings
using BDIO, ADerrors, ALPHAio
using QuadGK
import ADerrors: err
using Statistics

#================ SET UP PLOT PARAMETERS ==================#

rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
rcParams["text.usetex"] =  true
rcParams["mathtext.fontset"]  = "cm"
rcParams["font.size"] = 20
rcParams["axes.labelsize"] =26
rcParams["axes.titlesize"] = 22
plt.rc("text", usetex=true) # set to true if a LaTeX installation is present


include("../../utils/const.jl")
include("../../utils/IO_BDIO.jl")
include("../../Running/tools_running.jl")


path_res_physPoint = "/Users/alessandroconigli/MyDrive/postdoc-mainz/projects/deltalpha/physical_results/new_strategy"
path_SD = joinpath(path_res_physPoint, "Q_SD")
path_ID = joinpath(path_res_physPoint, "Q_ID")
path_LD = joinpath(path_res_physPoint, "Q_LD")

Qgev_SD = [4, 5, 6, 7, 8, 9, 12] # Q^2 GeV^2
Qgev_ID = Qgev_SD ./ 4
Qgev_LD = Qgev_SD ./ 16
Qgev_tot = vcat(Qgev_LD, Qgev_ID, Qgev_SD)
Qmgev =  9.0  # Q^2 GeV^2


## READING ISOVECTOR
# subtracted function
b33qm_pt_5loops = uwreal([1896.11*1e-5, 2.67*1e-5], "b33_qm_pt") # from wolfram code
b33qqm_SD = [(q/(4*Qmgev)) * b33qm_pt_5loops for q in Qgev_SD] # b33 at SD
b33qqm_ID = [(q/(4*Qmgev)) * b33qm_pt_5loops for q in Qgev_ID] # b33 at ID
b33qqm_LD = [(q/(4*Qmgev)) * b33qm_pt_5loops for q in Qgev_LD] # b33 at LD

# lattice data SD
pi33_SD_SD = read_phys_res(path_SD,"PI33_SD_physRes.bdio")
pi33_SD_ILD = read_phys_res(path_SD,"PI33_ILD_physRes.bdio")
pi33_SD = pi33_SD_SD + pi33_SD_ILD + b33qqm_SD

# lattice data ID
pi33_ID_sub = read_phys_res(path_ID, "PI33_ID_physRes.bdio")
pi33_ID = pi33_ID_sub + b33qqm_ID

# lattice data LD 
pi33_LD_sub = read_phys_res(path_LD, "PI33_LD_physRes.bdio")
pi33_LD = pi33_LD_sub + b33qqm_LD

# combining isovector
pi33 = pi33_SD + pi33_ID + pi33_LD
pi33_tot = vcat(pi33_LD, pi33_LD + pi33_ID, pi33_LD + pi33_ID + pi33_SD)

## READING ISOSCALAR
# lattice data SD
pi3388_SD = read_phys_res(path_SD, "PI3388_QSD_physRes.bdio")
pi88_SD = pi33_SD * charge_factor["88"] - pi3388_SD

#lattice data ID
pi3388_ID = read_phys_res(path_ID, "PI3388_QID_physRes.bdio")
pi88_ID = pi33_ID * charge_factor["88"] - pi3388_ID

# lattice data LD
pi88_LD = read_phys_res(path_LD, "PI88_LD_physRes.bdio")

# combining isoscalar
pi88 = (pi88_SD + pi88_ID + pi88_LD) ./ charge_factor["88"] # remove charge factor
pi88_tot = vcat(pi88_LD, pi88_LD + pi88_ID, pi88_LD + pi88_ID + pi88_SD) ./ charge_factor["88"]

## READING CHARM CONNECTED
bcc_raw = read_phys_res(path_SD, "PI_blc_2Qm_physRes.bdio")[1]
bcc_SD = [(q/(4*Qmgev)) * bcc_raw for q in Qgev_SD] # bcc at SD
bcc_ID = [(q/(4*Qmgev)) * bcc_raw for q in Qgev_ID] # bcc at ID
bcc_LD = [4/3 * (q/(4*Qmgev)) * bcc_raw for q in Qgev_LD] # bcc at LD

# lattice data SD
picc_SD_sub = read_phys_res(path_SD, "PIcc_conn_QSD_physRes.bdio")
picc_SD = picc_SD_sub .+ bcc_SD .+ 2 .* charge_factor["cc"] .* b33qqm_SD

# lattice data ID
picc_ID_sub = read_phys_res(path_ID, "PIcc_conn_QID_physRes.bdio")
picc_ID = picc_ID_sub .+ bcc_ID .+ 2 .* charge_factor["cc"] .* b33qqm_ID

# lattice data LD
picc_LD_sub = read_phys_res(path_LD, "PIcc_conn_QLD_physRes.bdio")
picc_LD = picc_LD_sub .+ bcc_LD .+ 2 .* charge_factor["cc"] .* b33qqm_LD

# combining charm connected
picc = (picc_SD + picc_ID + picc_LD) ./ charge_factor["cc"]
picc_tot = vcat(picc_LD, picc_LD + picc_ID, picc_LD + picc_ID + picc_SD) ./ charge_factor["cc"]

## READING MIXIN ISOSCALAR
# lattice data SD
pi08_SD = read_phys_res(path_SD, "PI08_Q_SD_physRes.bdio")

# lattice data ID
pi08_ID = read_phys_res(path_ID, "PI08_Q_ID_physRes.bdio")

# lattice data LD
pi08_LD = read_phys_res(path_LD, "PI08_Q_LD_physRes.bdio")

# combining pi08
pi08 = (pi08_SD + pi08_ID + pi08_LD) ./ charge_factor["08"]
pi08_tot = vcat(pi08_LD, pi08_LD + pi08_ID, pi08_LD + pi08_ID + pi08_SD) ./ charge_factor["08"]

## add scale setting error
add_t0_err!(pi33_tot, t0sqrt_ph_err); uwerr.(pi33_tot)
add_t0_err!(pi88_tot, t0sqrt_ph_err); uwerr.(pi88_tot)
add_t0_err!(picc_tot, t0sqrt_ph_err); uwerr.(picc_tot)
add_t0_err!(pi08_tot, t0sqrt_ph_err); uwerr.(pi08_tot)


##
#####################
## FIT to PI HVP
#####################

func33 = rational_func(M=3, N=3)
func88 = rational_func(M=3, N=3)
funccc = rational_func(M=3, N=3)
func08 = rational_func(M=2, N=3)

fit33 = fit_routine(func33, Qgev_tot, pi33_tot, 12)
fit88 = fit_routine(func88, Qgev_tot, pi88_tot, 12)
fitcc = fit_routine(funccc, Qgev_tot, picc_tot, 12)
fit08 = fit_routine(func08, Qgev_tot, pi08_tot, 10)

## plot
fig = figure(figsize=(10,7.))
ax = gca()
ax.tick_params(right=true)

xarr = range(0,12,length=100)

# pi33
errorbar(Qgev_tot, value.(pi33_tot), err.(pi33_tot), fmt="d", color="#009E73")
yarr33 = func33(xarr, fit33.param); uwerr.(yarr33)
fill_between(xarr, value.(yarr33).-err.(yarr33), value.(yarr33).+err.(yarr33), alpha=0.5, color="#009E73", label=L"$I=1$")
plot(xarr, value.(yarr33), color="#009E73")

# pi88
errorbar(Qgev_tot, value.(pi88_tot), err.(pi88_tot), fmt="d", color="#0072B2")
yarr88 = func88(xarr, fit88.param); uwerr.(yarr88)
fill_between(xarr, value.(yarr88).-err.(yarr88), value.(yarr88).+err.(yarr88), alpha=0.5, color="#0072B2", label=L"$I=0$")
plot(xarr, value.(yarr88), color="#0072B2")

# picc
errorbar(Qgev_tot, value.(picc_tot), err.(picc_tot), fmt="d", color="#D55E00")
yarrcc = funccc(xarr, fitcc.param); uwerr.(yarrcc)
fill_between(xarr, value.(yarrcc).-err.(yarrcc), value.(yarrcc).+err.(yarrcc), alpha=0.5, label=L"$\mathrm{charm}$", color="#D55E00")
plot(xarr, value.(yarrcc), color="#D55E00")

# pi08
errorbar(Qgev_tot, value.(pi08_tot), err.(pi08_tot), fmt="d", color="#E69F00")
yarr08 = func08(xarr, fit08.param); uwerr.(yarr08)
fill_between(xarr, value.(yarr08).-err.(yarr08), value.(yarr08).+err.(yarr08), alpha=0.5, label=L"$08$", color="#E69F00")
plot(xarr, value.(yarr08), color="#E69F00")



xlabel(L"$Q^2\ [\mathrm{GeV}^2]$")
ylabel(L"$\bar{\Pi}(Q^2)$")
legend()
tight_layout()
display(fig)
# savefig(joinpath(path_plot, "hvp_running_highq.pdf"))
close("all")
##