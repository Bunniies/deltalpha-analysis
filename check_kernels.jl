using Revise, HVPobs
using PyPlot, LaTeXStrings
using BDIO, ADerrors, ALPHAio
import ADerrors: err

#================ SET UP PLOT PARAMETERS ==================#

rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
rcParams["text.usetex"] =  true
rcParams["mathtext.fontset"]  = "cm"
rcParams["font.size"] =13
rcParams["axes.labelsize"] =22
rcParams["axes.titlesize"] = 18
plt.rc("text", usetex=true) # set to true if a LaTeX installation is present


include("./utils/const.jl")
include("./utils/types.jl")
include("./utils/plot_utils.jl")
include("./utils/IO_BDIO.jl")
include("./utils/tools.jl")

path_plot = "/Users/alessandroconigli/MyDrive/postdoc-mainz/projects/deltalpha/plots/"


##
T = range(0,4, length=1000)
Q = 9. / hc^2  *1e6
Qm = 36. / hc^2  *1e6
k1 = krnl_dα.(T, Q/4)./ T.^3
k2 = krnl_dα_sub.(T, Q/4, Qm)./ T.^3
k3 = krnl_dα_qhalf.(T, Q)./ T.^3
k4 = krnl_dα_qhalf_sub.(T, Q, Qm)./ T.^3
replace!(k1, NaN =>0)
replace!(k2, NaN =>0)
replace!(k3, NaN =>0)
replace!(k4, NaN =>0)

plot(T, k1, label=L"$K_1$", color="red")
plot(T, k2, label=L"$K_2^{\mathrm{sub}}$", color="orange")
plot(T, k3, label=L"$K_3$", color="navy")
plot(T, k4, label=L"$K_4^{\mathrm{sub}}$", color="lightskyblue")
legend(ncol=2)
xlim(-0.1,2.)
xlabel(L"$t \ [\mathrm{fm}]$")
ylabel(L"$K(t,Q^2)/t^3 \ [\mathrm{fm}^{-1}]$")
tight_layout()
display(gcf())
savefig(joinpath(path_plot, "kernel_comparisons.pdf"))
close()