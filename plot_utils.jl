function plot_tmr(tmr::Vector{uwreal}, ens::EnsInfo, q2::Float64; p::Union{String,Nothing}=nothing, info::String="_")
    tvals = collect(0:Int64(length(tmr))-1)
    yy = tmr ./ sqrt(8*t0(ens.beta))
    uwerr.(yy)
    errorbar(tvals, abs.(value.(yy)), err.(yy), fmt="s", mfc="none", color="royalblue", label=L"$Q^2= $"*"$(q2) "*L"$\mathrm{GeV}^2$" )
    xlabel(L"$t/a$")
    ylabel(L"$G(t) \cdot K(t, Q^2)$")
    legend()
    tight_layout()
    display(gcf())
    if !isnothing(p)
        tt = string("tmr_$(ens.id)_Q", Int64(q2), "_", info, ".pdf")
        savefig(joinpath(p, ens.id, tt))
    end
    close("all")

end

function plot_tmr_multiq(tmr::Vector{Vector{uwreal}}, ens::EnsInfo, q2::Vector{Float64}; p::Union{String,Nothing}=nothing, info::String="_")
    tvals = collect(0:Int64(length(tmr[1]))-1)
    cc=["forestgreen", "royalblue", "red"]
    for k in eachindex(tmr)
        yy = tmr[k] ./ sqrt.(8 .*t0.(ens.beta)); uwerr.(yy)
        errorbar(tvals, abs.(value.(yy)), err.(yy), fmt="o", color=cc[k], label=L"$Q^2=\ $"*"$(q2[k]) "*L"$\mathrm{GeV}^2$" )
        fill_between(tvals, abs.(value.(yy))- err.(yy), abs.(value.(yy))+ err.(yy), color=cc[k], alpha=0.2 )
    end
    xlabel(L"$t/a$")
    ylabel(L"$G(t) \cdot K(t, Q^2) / \sqrt{8t_0}$")
    legend()
    twiny()
    xx_fm =  tvals #.* value.(t0sqrt_ph ./ sqrt.(t0.(ens.beta)))
    fm_time = round.((xx_fm[1:6:end] ) .* value.(t0sqrt_ph ./ sqrt.(t0.(ens.beta))) , digits=1)
    xticks(xx_fm[1:6:end], fm_time)
    xlabel(L"$t\ [\mathrm{fm}]$", size=14)
    tight_layout()
    display(gcf())
    if !isnothing(p)
        tt = string("tmr$(info)_$(ens.id).pdf")
        savefig(joinpath(p, ens.id, tt))
    end
    close("all")
end