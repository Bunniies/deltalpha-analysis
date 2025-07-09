# Here we compute the FVC for delta alpha at some reference value. 
# We use the HP FVC up to t_ref, while we switch to MLL FVC beyond this value.

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

const path_fvc = "/Users/alessandroconigli/Lattice/data/HVP/FSE"
const path_bdio_obs = "/Users/alessandroconigli/MyDrive/postdoc-mainz/projects/deltalpha/data"
path_plot = "/Users/alessandroconigli/MyDrive/postdoc-mainz/projects/deltalpha/plots/new_strategy/ensembles/"

const path_store_fvc = "/Users/alessandroconigli/MyDrive/postdoc-mainz/projects/deltalpha/PIdata/impr_deriv/new_strategy/Q_SD/fvc_Lref"

#======= PHYSICAL CONSTANTS ====================#
const Qgev = [4, 5, 6, 7, 8, 9, 12]  # Q^2 # additional very high values
# const Qgev = [9.0]
const Qmgev = 9.0 # Qm^2


enslist = sort([ "H101", "H102", "N101", "C101", "C102", "D150",
          "B450", "N451", "N452", "D450", "D451", "D452",
         "N202", "N203", "N200", "D251", "D200", "D201", "E250",
          "J307", "J306", "J303", "J304", "E300", "F300",
         "J500", "J501"
])


ensinfo = EnsInfo.(enslist)

KRNLsub = krnl_dα_qhalf_sub # non-subtracted kernel
const Wind = Window("ILD") # SD or ILD distance window 

@warn("Always check the kernel that you are using!! \\ It has to match with the one used for the contribution you are considering")
#============== READ t0 FROM BDIO FILES =================#
dir_path = filter(x-> basename(x) in enslist, readdir(path_bdio_obs, join=true))
spectrum_path = vcat(filter(!isempty, [filter(x-> occursin("spectrum.bdio", x)  , readdir(dir_path[k], join=true)) for k in eachindex(dir_path)])...)

t0ens = Vector{uwreal}(undef, length(enslist))

for (k, ens) in enumerate(ensinfo)
    t0ens_aux = read_BDIO(spectrum_path[k], "spectrum", "t0")[1]
    t0ens[k] = t0ens_aux
end


##
#========== COMPUTE HP FVC ==========#

fvc_pi_hp = [Vector{Vector{uwreal}}(undef, length(Qgev)) for k in eachindex(ensinfo)]
fvc_k_hp = [Vector{Vector{uwreal}}(undef, length(Qgev)) for k in eachindex(ensinfo)]


for (k, ens) in enumerate(ensinfo)
    println("Ensemble: $(ens.id)")


    Qlat = Qgev .* t0sqrt_ph.^2 ./ t0ens[k] ./hc^2 *1e6
    qmlat = Qmgev  * t0sqrt_ph^2 / t0ens[k] /hc^2 *1e6

    # Linf
    fvc_raw_pi_Linf = get_fvc(joinpath(path_fvc, "fv_hp_Linf", "JKMPI_Mvmd"), ens.id)
    fvc_raw_k_Linf = get_fvc(joinpath(path_fvc, "fv_hp_Linf", "JKMK"), ens.id)
    # Lref
    fvc_raw_pi_Lref = get_fvc(joinpath(path_fvc, "fv_hp_Lref", "JKMPI_Mvmd"), ens.id)
    fvc_raw_k_Lref = get_fvc(joinpath(path_fvc, "fv_hp_Lref", "JKMK"), ens.id)

    # sum over the pion wrapping
    fvc_raw_pi_Linf = vcat(sum(fvc_raw_pi_Linf,dims=1)...)
    fvc_raw_k_Linf = vcat(sum(fvc_raw_k_Linf,dims=1)...)

    fvc_raw_pi_Lref = vcat(sum(fvc_raw_pi_Lref,dims=1)...)
    fvc_raw_k_Lref = vcat(sum(fvc_raw_k_Lref,dims=1)...)
    
    T = Float64.(collect(1:length(fvc_raw_pi_Linf)))

    for (j,q) in enumerate(Qlat)
        krnl = try
            KRNLsub.(T,q, qmlat)
        catch
            KRNLsub.(T,q)
        end

        # tmp = abs.((fvc_raw_pi_Linf .- fvc_raw_pi_Lref) .* krnl)
        fvc_pi_hp[k][j] = abs.((fvc_raw_pi_Linf .* krnl)) .- abs.((fvc_raw_pi_Lref .* krnl))
        fvc_k_hp[k][j]  = abs.((fvc_raw_k_Linf  .* krnl)) .- abs.((fvc_raw_k_Lref .* krnl))
    end
end
##
#========== COMPUTE MLL FVC ==========#
path_mll = joinpath(path_fvc, "fv_MLL")
fname = readdir(path_mll, join=true)
fname = filter!(x->occursin("Vref", basename(x)), fname) # take reference volume only
# fname = filter!(x->occursin("fvc.txt", basename(x)), fname) # take reference volume only

fvc_pi_mll = [Vector{Vector{uwreal}}(undef, length(Qgev)) for k in eachindex(ensinfo)]
for (k,ens) in enumerate(ensinfo)
    Qlat = Qgev .* t0sqrt_ph.^2 ./ t0ens[k] ./hc^2 *1e6
    qmlat = Qmgev  * t0sqrt_ph^2 / t0ens[k] /hc^2 *1e6

    p_ens = filter(x-> occursin(ens.id, basename(x)), fname)[1]
    println(p_ens)
    fname_tmp = readdlm(p_ens, comments=true)
    T = fname_tmp[:,1]

    corr_mll = Vector{uwreal}(undef, length(T))
    for t in eachindex(T)
        corr_mll[t] = uwreal([fname_tmp[t,2], fname_tmp[t,3]], ens.id*"_MLL_fve")
    end

    for (j,q) in enumerate(Qlat)
        krnl = try
            KRNLsub.(T,q, qmlat)
        catch
            KRNLsub.(T,q)
        end
        fvc_mll_tmr = -1 * corr_mll .* krnl
        # fvc_mll_tmr .= fvc_mll_tmr ./ value.(t0sqrt_ph) .* sqrt.(t0ens[k]) # this to transform to physical unit
        fvc_pi_mll[k][j] = fvc_mll_tmr
    end
end


## combine HP and MLL
pi_fvc_hp_tot = [Vector{uwreal}(undef, length(Qgev)) for k in eachindex(ensinfo)]
k_fvc_hp_tot = [Vector{uwreal}(undef, length(Qgev)) for k in eachindex(ensinfo)]
pi_fvc_hp_tstar = [Vector{uwreal}(undef, length(Qgev)) for k in eachindex(ensinfo)]
pi_fvc_mll = [Vector{uwreal}(undef, length(Qgev)) for k in eachindex(ensinfo)]
pi_fvc_hp_mll_tot = [Vector{uwreal}(undef, length(Qgev)) for k in eachindex(ensinfo)]


for (k,ens) in enumerate(ensinfo)
    T = HVPobs.Data.get_T(ens.id)
    tvals = collect(1:T)
    tvals_mll = collect(T-length(fvc_pi_mll[k][1]):T-1)


    for (j,q) in enumerate(Qgev)
        pi_hp_tot = sum(fvc_pi_hp[k][j] .* Wind(tvals .* value(a(ens.beta))))
        k_hp_tot = sum(fvc_k_hp[k][j] .* Wind(tvals .* value(a(ens.beta))))
        pi_hp_tstar = sum(fvc_pi_hp[k][j][1:tvals_mll[1]-1] .* Wind(collect(1:(tvals_mll[1]-1)) .* value(a(ens.beta))))
        
        pi_mll = sum(fvc_pi_mll[k][j] .* Wind(tvals_mll .* value(a(ens.beta))))

        tmr_combined = vcat(fvc_pi_hp[k][j][1:tvals_mll[1]-1], fvc_pi_mll[k][j] )
        pi_combined_tot = sum(tmr_combined .* Wind(tvals[1:end-1] .* value(a(ens.beta))))
        
        # store res
        pi_fvc_hp_tot[k][j] = pi_hp_tot
        k_fvc_hp_tot[k][j] = k_hp_tot
        pi_fvc_hp_tstar[k][j] = pi_hp_tstar
        pi_fvc_mll[k][j] = pi_mll
        pi_fvc_hp_mll_tot[k][j] = pi_combined_tot

    end

end

## Save in BDIO

@info("Saving FVE resulst in BDIO (HP&MLL for pion, HP for kaon)")
io = IOBuffer()
write(io, "FVC at Lref (HP&MLL for pion, HP for kaon)")

fb = ALPHAdobs_create(joinpath(path_store_fvc, "fve_Lref_subKernel_QSD_ILD.bdio"), io)
for (k,ens) in enumerate(ensinfo)
    extra = Dict{String, Any}("Ens" => ens.id)
    data = Dict{String, Array{uwreal}}(
        "fvc_pi" => pi_fvc_hp_mll_tot[k],
        "fvc_k"  => k_fvc_hp_tot[k]
    )
    ALPHAdobs_write(fb, data, extra=extra)
end 
ALPHAdobs_close(fb)
println("# Saving complete")

##
#========= TEST READING =========#
fb = BDIO_open(joinpath(path_store_fvc, "fve_Lref_subKernel_QLD.bdio"), "r")
res = Dict()
while ALPHAdobs_next_p(fb)
    d = ALPHAdobs_read_parameters(fb)
    nobs = d["nobs"]
    dims = d["dimensions"]
    sz = tuple(d["size"]...)
    extra = d["extra"]
    ks = collect(d["keys"])
    res[extra["Ens"]] = ALPHAdobs_read_next(fb, size=sz, keys=ks)
end
BDIO_close!(fb)


## Prepare latex table

io = open("/Users/alessandroconigli/Desktop/table_fve.txt", "w")
for (k,ens) in enumerate(ensinfo)
    if ens.kappa_l == ens.kappa_s
        total_correction = pi_fvc_hp_mll_tot[k][1] 
        write(io, ens.id, " & ", print_uwreal(pi_fvc_hp_mll_tot[k][1]*1e5), " & ", print_uwreal(pi_fvc_hp_tot[k][1]*1e5), " & ", "-", " & ", print_uwreal(total_correction*1e5), " \\\\ \n" )        
    else
        total_correction = pi_fvc_hp_mll_tot[k][1] + k_fvc_hp_tot[k][1]
        write(io, ens.id, " & ", print_uwreal(pi_fvc_hp_mll_tot[k][1]*1e5), " & ", print_uwreal(pi_fvc_hp_tot[k][1]*1e5), " & ", print_uwreal(k_fvc_hp_tot[k][1]*1e5), " & ", print_uwreal(total_correction*1e5), " \\\\ \n" )
    end
end
close(io)


## PLOT HP and MLL 


for (k,ens) in enumerate(ensinfo)

    for (j,q) in enumerate(Qgev[1:1])
        fig = figure(figsize=(10,7.5))
        T  = HVPobs.Data.get_T(ens.id)
        tvals = collect(1:T)
        uwerr.(fvc_pi_hp[k][j])
        errorbar(tvals, value.(fvc_pi_hp[k][j]), err.(fvc_pi_hp[k][j]), fmt="s", label="HP")

        tvals_mll = collect(T-length(fvc_pi_mll[k][j]):T-1)
        uwerr.(fvc_pi_mll[k][j])
        errorbar(tvals_mll, value.(fvc_pi_mll[k][j]), err.(fvc_pi_mll[k][j]), fmt="s", label="MLL")
        axvline(tvals_mll[1], ls="-.", color="black")

        legend()
        xlabel(L"$t/a$")
        ylabel(L"$G(t)K(t,Q^2)$")
        xlim(0, Int64(T/2))
        tight_layout()
        display(fig)
        #pp = joinpath(path_plot, ens.id, "fvc_hp_mll_1Gev2.pdf")
        #savefig(pp)
        close("all")

        res_hp = sum(fvc_pi_hp[k][j])*1e5; uwerr.(res_hp)
        res_mll = sum(fvc_pi_mll[k][j])*1e5; uwerr.(res_mll)
        res_hp_t_in = sum(fvc_pi_hp[k][j][1:tvals_mll[1]-1])*1e5; uwerr.(res_hp_t_in)
        combined = res_mll + res_hp_t_in; uwerr(combined)
        println("$(ens.id) at Q^2 = $(Qgev[j]):")
        println("HP:       ", res_hp)
        println("MLL:      ", res_mll)
        println("HP t<t_i: ", res_hp_t_in)
        println("combined: ", combined)
        println("\n")
    end
end


