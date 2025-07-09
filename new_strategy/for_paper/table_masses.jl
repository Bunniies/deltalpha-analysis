using Revise, HVPobs
using PyPlot, LaTeXStrings
using BDIO, ADerrors, ALPHAio
import ADerrors: err
using DelimitedFiles

#================ SET UP PLOT PARAMETERS ==================#

rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
rcParams["text.usetex"] =  true
rcParams["mathtext.fontset"]  = "cm"
rcParams["font.size"] =20
rcParams["axes.labelsize"] =26
rcParams["axes.titlesize"] = 22
plt.rc("text", usetex=true) # set to true if a LaTeX installation is present


include("../../utils/const.jl")
include("../../utils/types.jl")
include("../../utils/plot_utils.jl")
include("../../utils/IO_BDIO.jl")
include("../../utils/tools.jl")

const path_bdio_obs = "/Users/alessandroconigli/MyDrive/postdoc-mainz/projects/deltalpha/data"


enslist = [ "H101", "H102", "N101", "C101", "C102", "D150",
          "B450", "N451", "N452", "D450", "D451", "D452",
         "N202", "N203", "N200", "D251", "D200", "D201", "E250",
          "J307", "J306", "J303", "J304", "E300", "F300",
         "J500", "J501"
]
ensinfo = EnsInfo.(enslist)

#============== READ PS MASSES AND t0 FROM BDIO FILES =================#
dir_path = filter(x-> basename(x) in enslist, readdir(path_bdio_obs, join=true))
spectrum_path = vcat(filter(!isempty, [filter(x-> occursin("spectrum.bdio", x)  , readdir(dir_path[k], join=true)) for k in eachindex(dir_path)])...)

phi2 = Vector{uwreal}(undef, length(enslist))
phi4 = Vector{uwreal}(undef, length(enslist))
mpi_tot = Vector{uwreal}(undef, length(enslist))
mk_tot = Vector{uwreal}(undef, length(enslist))
a28t0 = Vector{uwreal}(undef, length(enslist))
t0ens = Vector{uwreal}(undef, length(enslist))

for (k, ens) in enumerate(ensinfo)

    path_ens = filter(x->occursin(ens.id, basename(x)), spectrum_path)[1]
    m_pi  = read_BDIO(path_ens, "spectrum", "mpi")[1]
    m_k   = read_BDIO(path_ens, "spectrum", "mk")[1]
    t0ens_aux = read_BDIO(path_ens, "spectrum", "t0")[1]

    phi2[k]  = 8 * t0ens_aux * m_pi^2
    phi4[k]  = 8 * t0ens_aux * (m_k^2 + 0.5*m_pi^2 )
    a28t0[k] =  a(ens.beta)^2 #1 / (8*t0ens_aux)
    t0ens[k] = t0ens_aux
    mpi_tot[k] = m_pi
    mk_tot[k] = m_k
end

## Prepare latex table
io = open("/Users/alessandroconigli/Desktop/table_masses.txt", "w")
for (k,ens) in enumerate(ensinfo)
    write(io, ens.id, " & ", print_uwreal(mpi_tot[k]), " & ", print_uwreal(mk_tot[k]), " & ", print_uwreal(t0ens[k]), " & ", print_uwreal(phi2[k]), " & ", print_uwreal(phi4[k]), " \\\\ \n")
end
close(io)

