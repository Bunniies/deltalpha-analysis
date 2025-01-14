using HVPobs
using ALPHAio, BDIO
using PyPlot, LaTeXStrings
using ADerrors
import ADerrors:err
using OrderedCollections
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

include("../utils/IO_BDIO.jl")
include("./tools_running.jl")
include("../utils/const.jl")

const Qgev = [0.05, 0.1, 0.4, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0] # Q^2
NMOM = length(Qgev)
Qmgev =  9.0  # GeV^2


## COMPARING LIGHT HVP BETWEEN LATTICE AND PT 
# Q^2 - Q^2/4 only
pi33_SD = read_phys_res(path_phys_res_highq, "PI33_SD_physRes.bdio")
pi33_ILD = read_phys_res(path_phys_res_highq, "PI33_ILD_physRes.bdio")

b33qm_pt_5loops = uwreal([1895.26*1e-5, 2.62*1e-5], "b33_qm_pt")
b33qqm_highq = [(q/(4*Qmgev)) * b33qm_pt_5loops for q in Qgev] # using PT result
# b33qqm_highq = [(q/(4*Qmgev)) * log(2) / (4pi^2) for q in Qgev] # at LO
pi33_highq = (pi33_ILD + pi33_SD + b33qqm_highq ) 

pi33_pt = OrderedDict(
    "2" => uwreal([0.03535463041256225,0.02739426584616154],"msbar_2" ),
    "3" => uwreal([0.0194233109053707,0.002378607983786594],"msbar_3" ),
    "4" => uwreal([0.01965668758406193,0.0014002562074296404],"msbar_4" ),
    "5" => uwreal([0.020218076452051727,0.0010287575253706352],"msbar_5" ),
    "6" => uwreal([0.019637641535117285,0.0006766871504073361],"msbar_6" ),
    "7" => uwreal([0.01951275508421394,0.0004256577456237366],"msbar_7" ),
    "8" => uwreal([0.01996836352106275,0.0002962789634870446],"msbar_8" ),
    "9" => uwreal([0.019653422492354532,0.00024453913473834613], "msbar_9")
)


## plotting

uwerr.(pi33_highq)
uwerr.(values(pi33_pt))

fig = figure(figsize=(10,8))

errorbar(Qgev, value.(pi33_highq), err.(pi33_highq),ms=10,mfc="none", fmt="d", label=L"$\mathrm{Lattice}$")
errorbar(collect(2:9).+0.1, value.(values(pi33_pt)), err.(values(pi33_pt)), fmt="^", mfc="none", ms=10, label=L"$\mathrm{AdlerPy}$")


xlim(2.5,9.5)
ylim(0.0150, 0.0250)
xlabel(L"$Q^2 \ [\mathrm{GeV}^2]$")
ylabel(L"$\widehat{\mathit{\Pi}}^{(3,3)}(Q^2)$")
legend()
tight_layout()
savefig(joinpath(path_plot, "isovector_lattice_vs_pt.pdf"))
display(fig)
## COMPARING CHARM HVP BETWEEN LATTICE AND PT USING BOTH MSbar AND POLE MASSES
# Q^2 - Q^2/4 only

b33qm_pt_5loops = uwreal([1899.07*1e-5, 1.88*1e-5], "b33_qm_pt")
b33qqm_highq = [(q/(4*Qmgev)) * b33qm_pt_5loops for q in Qgev] # using PT result

# cc connected 
picc_highq_sub = read_phys_res(path_phys_res_highq, "PIcc_conn_physRes.bdio")
bcc_highq = read_phys_res(path_phys_res_highq, "PI_blc_2Qm_physRes.bdio")

picc_highq = picc_highq_sub .+ Qgev ./ (4*Qmgev) .* bcc_highq .+ 2 * charge_factor["cc"] * b33qqm_highq
#picc_highq ./= charge_factor["cc"]


picc_pt_msbar = OrderedDict(
    "0.05" => uwreal([6.735441697088523e-05, 8.859290467293803e-07 ],"msbar_005" ),
    "0.1" => uwreal([0.00013872089734313654, 1.99220095232246e-06],"msbar_01" ),
    "0.4" => uwreal([0.0005341238354721566, 6.059084379351581e-06],"msbar_04" ),
    "1" => uwreal([0.0012592981987082207,1.589126679440487e-05],"msbar_1" ),
    "2" => uwreal([0.002315894341390204,3.0233340570720026e-05],"msbar_2" ),
    "3" => uwreal([0.0033405656219024386, 3.8557270538983e-05],"msbar_3" ),
    "4" => uwreal([0.004107619292948906, 4.545520251470209e-05],"msbar_4" ),
    "5" => uwreal([0.004948131331661961,5.562605250991618e-05],"msbar_5" ),
    "6" => uwreal([0.005580640181799866,5.23016308104534e-05],"msbar_6" ),
    "7" => uwreal([0.006014283692962188, 4.885852633871269e-05],"msbar_7" ),
    "8" => uwreal([0.006674005538554557,5.301710047844647e-05],"msbar_8" ),
    "9" => uwreal([0.0070832138989291036, 5.5273640853745014e-05], "msbar_9")
)
picc_pt_mpole = OrderedDict(
    "0.05" => uwreal([6.228681445751013e-05, 5.657351504776487e-06 ],"mpole_005" ),
    "0.1" => uwreal([0.00011367502970462019, 8.82665459476582e-06 ],"mpole_01" ),
    "0.4" => uwreal([0.00048361923643717823, 4.874535626191511e-05 ],"mpole_04" ),
    "1" => uwreal([0.0011858989352538218, 9.728398802148765e-05 ],"mpole_1" ),
    "2" => uwreal([0.0021487118760959334,0.00018089109949931983],"mpole_2" ),
    "3" => uwreal([0.0030842685214353027, 0.00024762136093263985],"mpole_3" ),
    "4" => uwreal([0.003985309273003494,0.0002377454001919895],"mpole_4" ),
    "5" => uwreal([0.004503952569116253, 0.00027141204420679866],"mpole_5" ),
    "6" => uwreal([0.005086974324212858, 0.0003563353844951602],"mpole_6" ),
    "7" => uwreal([0.005640338337626476, 0.000303234517332303],"mpole_7" ),
    "8" => uwreal([0.0060831947308964485,0.00031551267981489876],"mpole_8" ),
    "9" => uwreal([0.006674360139956521, 0.00027009063951378775], "mpole_9")
)

## plotting

uwerr.(picc_highq)
uwerr.(values(picc_pt_mpole))
uwerr.(values(picc_pt_msbar))

fig = figure(figsize=(10,8))

errorbar(Qgev, value.(picc_highq), err.(picc_highq),ms=10,mfc="none", fmt="d", label="Lattice")

errorbar(Qgev.+0.1, value.(values(picc_pt_msbar)), err.(values(picc_pt_msbar)), fmt="^", mfc="none", ms=10, label=L"$\mathcal{O}(\alpha_s^3) \ \bar{\mathrm{MS}}$")

errorbar(Qgev.+0.2, value.(values(picc_pt_mpole)), err.(values(picc_pt_mpole)), fmt="s", mfc="none", ms=10, label=L"$\mathcal{O}(\alpha_s^3)\ \mathrm{OS}$")

xlim(1.5,9.5)
xlabel(L"$Q^2 \ [\mathrm{GeV}^2]$")
ylabel(L"$\widehat{\mathit{\Pi}}^{(c,c)}(Q^2)$")
legend()
tight_layout()
savefig(joinpath(path_plot, "charm_lattice_vs_pt.pdf"))
display(fig)

## COMPARE SUBTRACTED ISOVECTOR PIECE b33

b33qm_pt_5loops = uwreal([1895.26*1e-5, 2.62*1e-5], "b33_qm_pt") *1e5
b33qm_pt_5loops_wolfram = uwreal([1896.48*1e-5, 0.0], "b33_qm_pt_wolfram")

b33qm_allq = [(q/(4*Qmgev)) * b33qm_pt_5loops for q in Qgev]
b33qm_allq_wolfram = [(q/(4*Qmgev)) * b33qm_pt_5loops_wolfram for q in Qgev]
uwerr.(b33qm_allq)

fig = figure(figsize=(10,8))

errorbar(Qgev, value.(b33qm_allq), err.(b33qm_allq), fmt="s", ms=10,mfc="none", label="AdlerPy")
plot(Qgev, value.(b33qm_allq_wolfram), mfc="none", ms=10, label="Wolfram")
xlim(1.5,9.5)
xlabel(L"$Q^2 \ [\mathrm{GeV}^2]$")
ylabel(L"$b^{(3,3)(Q^2, Q_m^2)}$")
legend()
tight_layout()
savefig(joinpath(path_plot, "b33_lattice_vs_pt.pdf"))
display(fig)

## COMPARE SUBTRACTRED CHARM PIECE ΔlcB
bcc_highq = read_phys_res(path_phys_res_highq, "PI_blc_2Qm_physRes.bdio")
