
function krnl_dα(t::Union{Float64, uwreal}, Q2::Union{Float64,uwreal})

    K = t^2 - 4/Q2 * sin((sqrt(Q2)*t)/2)^2
    return K
end

function krnl_dα_qhalf(t::Union{Float64, uwreal}, Q2::Union{Float64,uwreal})

    K = 16. / Q2 * sin((sqrt(Q2)*t)/4)^4
    return K
end
function krnl_dα_qhalf_sub(t::Union{Float64, uwreal}, Q2::Union{Float64,uwreal}, Q2m::Union{Float64,uwreal})
    K = 16. / Q2 * sin((sqrt(Q2)*t)/4)^4 - (Q2/Q2m^2) * sin((sqrt(Q2m)*t)/2)^4
    return K
end

function tmr_integrand(cc::Vector{uwreal}, q2::Union{Float64, uwreal}, krnl::Function; pl::Bool=true, t0ens::Union{uwreal, Nothing}=nothing, wind::Union{Nothing, Window}=nothing, data::Bool=false)
    
    if isnothing(t0ens) && !isnothing(wind)
        error("You assigned a window without passign t0ens. \n Conversion to t[fm] failed.")
    end
    T = length(cc)
    x0 = Float64.(collect(0:T-1)) 
    k = krnl.(x0, q2) #.* value.(sqrt.(t0ens) ./ t0sqrt_ph ).^(-2)
    
    integrand = cc .* k  #abs.(cc .* k)
    x0 = x0 .* value(t0sqrt_ph / sqrt(t0ens))
    if !isnothing(wind)
        integrand = integrand .* wind(x0)
    end


    Thalf = try
        Int64(T/2)
   catch
        Int64((T+1)/2)
   end
    if pl        
        uwerr.(integrand)
        if !isnothing(t0ens)
            # x0 = x0 .* value(t0sqrt_ph / sqrt(t0ens))
            xlabel(L"$t \ [\mathrm{fm}]$")
        else
            xlabel(L"$t/a$")
        end
        errorbar(x0[1:Thalf], value.(integrand)[1:Thalf], err.(integrand)[1:Thalf], fmt="s", mfc="none", color="forestgreen", capsize=2)
        ylabel(L"$G(t) \cdot K(t, Q^2)$")
        xlim(-0.05, 4)
        display(gcf())
        close("all")
    end
    !data ? (return sum(integrand[1:Thalf])) : (return sum(integrand[1:Thalf]), integrand[1:Thalf])
    return   integrand[1:Thalf]
end
tmr_integrand(cc::Corr, q2::Union{Float64, uwreal}, krnl::Function; pl::Bool=true, t0ens::Union{uwreal, Nothing}=nothing, wind::Union{Nothing, Window}=nothing, data::Bool=false) = tmr_integrand(cc.obs, q2, krnl; pl=pl, t0ens=t0ens, wind=wind, data=data)

function tmr_integrand(cc::Vector{uwreal}, q2::Union{Float64, uwreal}, q2m::Union{Float64, uwreal}, krnl::Function; pl::Bool=true, t0ens::Union{uwreal, Nothing}=nothing, wind::Union{Nothing, Window}=nothing, data::Bool=false)
    if isnothing(t0ens) && !isnothing(wind)
        error("You assigned a window without passign t0ens. \n Conversion to t[fm] failed.")
    end
    T = length(cc)
    x0 = Float64.(collect(0:T-1)) 
    k = krnl.(x0, q2, q2m) #.* value.(sqrt.(t0ens) ./ t0sqrt_ph ).^(-2)
    
    integrand = cc .* k # abs.(cc .* k)
    x0 = x0 .* value(t0sqrt_ph / sqrt(t0ens))
    if !isnothing(wind)
        integrand = integrand .* wind(x0)
    end


    Thalf = try
        Int64(T/2)
   catch
        Int64((T+1)/2)
   end
    if pl        
        uwerr.(integrand)
        if !isnothing(t0ens)
            # x0 = x0 .* value(t0sqrt_ph / sqrt(t0ens))
            xlabel(L"$t \ [\mathrm{fm}]$")
        else
            xlabel(L"$t/a$")
        end
        errorbar(x0[1:Thalf], value.(integrand)[1:Thalf], err.(integrand)[1:Thalf], fmt="s", mfc="none", color="forestgreen", capsize=2)
        ylabel(L"$G(t) \cdot K(t, Q^2)$")
        xlim(-0.05, 4)
        display(gcf())
        close("all")
    end
    !data ? (return sum(integrand[1:Thalf])) : (return sum(integrand[1:Thalf]), integrand[1:Thalf])
    return   integrand[1:Thalf] 
end
tmr_integrand(cc::Corr, q2::Union{Float64, uwreal}, q2m::Union{Float64, uwreal}, krnl; pl::Bool=true, t0ens::Union{uwreal, Nothing}=nothing, wind::Union{Nothing, Window}=nothing, data::Bool=false) = tmr_integrand(cc.obs, q2, q2m, krnl; pl=pl, t0ens=t0ens, wind=wind, data=data)


function tmr_integrand_3l(cc::Vector{Tu}, q2::Union{Tu, Tf}, q2m::Union{Tu, Tf}, krnl::Function; pl::Bool=true, t0ens::Union{Tu, Nothing}=nothing, wind::Bool=false, d=0.4, delta=0.15 ) where {Tu, Tf}
    
    T = length(cc)

    x0 = Float64.(collect(0:T-1)) 
    k = krnl.(x0, q2, q2m) #.* value.(sqrt.(t0ens) ./ t0sqrt_ph ).^(-2)
    x0 = Float64.(collect(0:T-1)) .* value(t0sqrt_ph / sqrt(t0ens))
    
    integrand = cc .* k
    
    Thalf = try
        Int64(T/2)
   catch
        Int64((T+1)/2)
   end

    if pl        
        uwerr.(integrand)

        cumtmr = cumsum(integrand)
        plot(x0, value.(cumtmr))
        xlabel(L"$t \ [\mathrm{fm}]$")
        ylabel(L"$\sum G(t) \cdot K(t, Q^2)$")
        if wind
            axvline(d, ls="dashed", color="black", lw=0.8)
        end
        display(gcf())
        close()

        errorbar(x0[1:Thalf], abs.(value.(integrand))[1:Thalf], err.(integrand)[1:Thalf], fmt="s", mfc="none", color="forestgreen", capsize=2)
        xlabel(L"$t/a$")
        ylabel(L"$G(t) \cdot K(t, Q^2)$")
        if wind
            axvline(d, ls="dashed", color="black", lw=0.8)
        end
        display(gcf())
        close("all")
    end

    if wind
        w = window_krnl.(x0, d=d, delta=delta)
        integrand .*=  w
    end
    return   integrand[1:Thalf] * sign(value(integrand[5]))
end

function window_krnl(x0; d=0.4, delta=0.15)
    w = 1 - 0.5 * (1 + tanh((x0-d)/delta))
    return w
end

function treelevel_continuum_correlator(t)
    return 1. / (2* pi^2 * t^3)
end

function get_meff_BMA(corr::Vector{Corr}, ens::EnsInfo; pl::Bool=true, 
    wpm::Union{Dict{Int64,Vector{Float64}},Dict{String,Vector{Float64}}, Nothing}=nothing, path_plt::Union{String,Nothing}=nothing, ll::LaTeXString=L"$m_{\pi}$")
   
    M = Vector{uwreal}(undef, length(corr))
    if !isnothing(path_plt) 
        path_plt = joinpath(path_plt, ens.id)  
    end

    
 
    if ens.id in ["A653", "A654", "D150", "B450", "N451", "D450", "D451", "D452", "D251", "E250"]
        taux = collect(2:15) 
        Thalf = Int(length(corr[1].obs)/2)+1
        T = length(corr[1].obs)
        tmin = Thalf .- taux
        tmax = Thalf .+ taux
        tmax = tmax[1:2:end]
        tmin = tmin[1:2:end] 
        @. const_fit(x,p) = p[2] * (exp(-p[1]*(x)) + exp(-p[1]*(T-x)) )# + p[4] * (exp(-p[3]*(x)) + exp(-p[3]*(T-x)) )
        param=2
        pbc = true
    else
        # @. const_fit(x, p) = p[1] + 0 * x
        param = 1 

        # tmin and tmax tuned on H101 and then scaled 
        tminaux = collect(18:44) #* round(Int64, sqrt(value(t0(ens.beta) / t0(3.4))))
        T = length(corr[1].obs)
        tmaxaux = collect(20:40) #* round(Int64, sqrt(value(t0(ens.beta) / t0(3.4))))
        tmax = T .- collect(reverse(tmaxaux))
        tmax = tmax[1:2:end]
        tmin = tminaux[1:2:end] 
        pbc = false
    end

    println("tmax: ", tmax)
    println("tmin: ", tmin)

    for k in eachindex(corr)
        t = string(L"$\kappa_1 = $" ,L" $\kappa_2 = $" )   
        if ens.id in ["A653", "A654","D150", "B450", "N451", "D450", "D451", "D452", "D251", "E250"]
            data = -1. .* corr[k].obs
        else
            _, data = meff(corr[k], [tmin[end-4],tmax[1]], pl=false, wpm=wpm, data=true)
        end

        name_bma = ll*"_BMA.pdf"
        isnothing(path_plt) ? tt = nothing : tt = joinpath(path_plt, name_bma)

        M[k], syst, all_res, FitDict = bayesian_av(const_fit, data, tmin, tmax, param, pl=pl, wpm=wpm, data=true, path_plt=tt, plt_title=t , label=ll, pbc=pbc)
    
        if pl
            fig = figure(figsize=(10,7.5))
            modstot = FitDict["mods"]
            AIC  = FitDict["AIC"]
            Wtot = FitDict["W"]
            IDX  = partialsortperm(AIC, 1:3)
            mods = modstot[IDX]
            main_res = all_res[IDX]
            W = Wtot[IDX]
            isnothing(wpm) ? uwerr.(data) : [uwerr(data[i],wpm) for i in eachindex(data)]
            vv = value.(data)
            ee = err.(data)
            xx = collect(0:length(vv)-1) 
            errorbar(xx, vv, ee, fmt="d", mfc="none", color="black", capsize=2)
            CC = ["forestgreen", "royalblue", "tomato"]
            for (ii,trange) in enumerate(mods)
                isnothing(wpm) ? uwerr(main_res[ii]) : uwerr(main_res[ii], wpm)
                vaux = value(main_res[ii])
                eaux = err(main_res[ii])
                regex = r"([0-9]+),([0-9]+)"
                aux = match(regex, trange)
                fill_between(parse(Float64,aux[1]):parse(Float64,aux[2]), vaux-eaux, vaux+eaux, color=CC[ii], alpha=0.5, label=round(W[ii], digits=3))
            end
            
            axhline(y=value(main_res[1])-14*err(main_res[1]), xmin=(tmin[1])/length(xx), xmax=(tmin[end])/length(xx), ls="dashed", color="green")
            axhline(y=value(main_res[1])-14.3*err(main_res[1]), xmin=(tmax[1])/length(xx), xmax=(tmax[end])/length(xx), ls="dashed", color="red")
            text(tmin[1],value(main_res[1])-13*err(main_res[1]), L"$\{t_{\mathrm{min}}\}$" )
            text(tmax[1],value(main_res[1])-13*err(main_res[1]), L"$\{t_{\mathrm{max}}\}$" )

            legend(title="Weights")
            xlabel(L"$t/a$")
            ylabel(ll)
            xlim(xx[1], xx[end])
            ylim(value(main_res[1]) -20*err(main_res[1]), value(main_res[1]) +20*err(main_res[1]))
            twiny()
            xlim(xx[1], xx[end])
            fm_time = round.((xx)[1:10:end] .* value(a(ens.beta)), digits=2)
            xticks(xx[1:10:end], fm_time)
            # xticks()
            xlabel(L"$t-t_0\ [\mathrm{fm}]$", fontsize=14)
            display(fig)
            if !isnothing(path_plt)
                savefig(joinpath(path_plt, ll*"_plateau.pdf"))
            end
            close("all")
        end
    end
    return M
end


function get_dec_const_BMA(corr_pp::Corr, corr_a0p::Corr, mps::uwreal, ens::EnsInfo; pl::Bool=true, 
    wpm::Union{Dict{Int64,Vector{Float64}},Dict{String,Vector{Float64}}, Nothing}=nothing, path_plt::Union{String,Nothing}=nothing, ll::LaTeXString=L"$m_{\pi}$")
   
    if !isnothing(path_plt) 
        path_plt = joinpath(path_plt, ens.id)  
    end

    
 
    if ens.id in ["A653", "A654","B450", "E250"]
        taux = collect(15:40) 
        Thalf = Int(length(corr_pp.obs)/2)+1
        T = length(corr_pp.obs)
        tmin = Thalf .- taux
        tmax = Thalf .+ taux
        tmax = tmax[1:2:end]
        tmin = tmin[1:2:end] 
        @. fit_corr(x,p) = p[1] * (exp(-p[2]*(x)) + exp(-p[2]*(T-x)) )
        @. fit_corr_a0p(x,p) = p[1] * (exp(-p[2]*(x)) - exp(-p[2]*(T-x)) )

        param=2
        pbc = true
    else
        #@. const_fit(x, p) = p[1] + 0 * x
        param = 1 

        # tmin and tmax tuned on H101 and then scaled 
        tminaux = collect(28:55) #* round(Int64, sqrt(value(t0(ens.beta) / t0(3.4))))
        T = length(corr_pp.obs)
        tmaxaux = collect(30:60) #* round(Int64, sqrt(value(t0(ens.beta) / t0(3.4))))
        tmax = T .- collect(reverse(tmaxaux))
        tmax = tmax[1:2:end]
        tmin = tminaux[1:2:end] 
        pbc = false
    end

    println("tmax: ", tmax)
    println("tmin: ", tmin)

    t = string(L"$\kappa_1 = $" ,L" $\kappa_2 = $" )   
    name_bma = ll*"_BMA.pdf"
    isnothing(path_plt) ? tt = nothing : tt = joinpath(path_plt, name_bma)
    if ens.id in [ "A653", "A654","B450", "E250"]
        data_pp = -1. .* corr_pp.obs
        data_a0p =  corr_a0p.obs
        idx = findall(x->x==0.0, value.(data_a0p))
        data_a0p[idx[2]] = uwreal([0.0, 1e-10], "test")
        idx = findall(x->x==0.0, value.(data_a0p))
        uwerr.(data_a0p)
        errorbar(collect(1:T), value.(data_a0p), err.(data_a0p), fmt="s")
        display(gcf())
        close("all")
        cp, syst, all_res, FitDict = bayesian_av(fit_corr, data_pp, tmin, tmax, param, pl=pl, wpm=wpm, data=true, path_plt=tt, plt_title=t , label=ll, pbc=pbc)
        caa, syst, all_res, FitDict = bayesian_av(fit_corr_a0p, data_a0p, tmin, tmax, param, pl=pl, wpm=wpm, data=true, path_plt=tt, plt_title=t , label=ll, pbc=pbc)

        f =  sqrt(2) / sqrt(mps) * caa / sqrt(cp)
    else
        _, data = meff(corr[k], [tmin[end-4],tmax[1]], pl=false, wpm=wpm, data=true)
    end


    if pl
        data = data_pp
        fig = figure(figsize=(10,7.5))
        modstot = FitDict["mods"]
        AIC  = FitDict["AIC"]
        Wtot = FitDict["W"]
        IDX  = partialsortperm(AIC, 1:3)
        mods = modstot[IDX]
        main_res = all_res[IDX]
        W = Wtot[IDX]
        isnothing(wpm) ? uwerr.(data) : [uwerr(data[i],wpm) for i in eachindex(data)]
        vv = value.(data)
        ee = err.(data)
        xx = collect(0:length(vv)-1) 
        errorbar(xx, vv, ee, fmt="d", mfc="none", color="black", capsize=2)
        CC = ["forestgreen", "royalblue", "tomato"]
        for (ii,trange) in enumerate(mods)
            isnothing(wpm) ? uwerr(main_res[ii]) : uwerr(main_res[ii], wpm)
            vaux = value(main_res[ii])
            eaux = err(main_res[ii])
            regex = r"([0-9]+),([0-9]+)"
            aux = match(regex, trange)
            fill_between(parse(Float64,aux[1]):parse(Float64,aux[2]), vaux-eaux, vaux+eaux, color=CC[ii], alpha=0.5, label=round(W[ii], digits=3))
        end
        
        axhline(y=value(main_res[1])-14*err(main_res[1]), xmin=(tmin[1])/length(xx), xmax=(tmin[end])/length(xx), ls="dashed", color="green")
        axhline(y=value(main_res[1])-14.3*err(main_res[1]), xmin=(tmax[1])/length(xx), xmax=(tmax[end])/length(xx), ls="dashed", color="red")
        text(tmin[1],value(main_res[1])-13*err(main_res[1]), L"$\{t_{\mathrm{min}}\}$" )
        text(tmax[1],value(main_res[1])-13*err(main_res[1]), L"$\{t_{\mathrm{max}}\}$" )
        legend(title="Weights")
        xlabel(L"$t/a$")
        ylabel(ll)
        xlim(xx[1], xx[end])
        ylim(value(main_res[1]) -20*err(main_res[1]), value(main_res[1]) +20*err(main_res[1]))
        twiny()
        xlim(xx[1], xx[end])
        fm_time = round.((xx)[1:12:end] .* value(a(ens.beta)), digits=2)
        xticks(xx[1:12:end], fm_time)
        # xticks()
        xlabel(L"$t-t_0\ [\mathrm{fm}]$", fontsize=14)
        display(fig)
        if !isnothing(path_plt)
            savefig(joinpath(path_plt, ll*"_plateau.pdf"))
        end
        close("all")
    end
    return f
end


function get_TIC(fitcat::FitCat)
    chi2tot = getfield.(fitcat.fit, :chi2)
    chi2exptot = getfield.(fitcat.fit, :chi2exp)
    ic = chi2tot .- 2 .* chi2exptot
    return ic
end

function get_bma_weight(ic::Vector{Float64})
    w = exp.(-0.5 .* ic)
    w = w ./ sum(w)
    return w
end

function get_w_from_fitcat(catvec::Vector{FitCat}; norm::Bool=false)
    ic = get_TIC.(catvec)
    ic = vcat(ic...)
    w = get_bma_weight(ic)
    if norm
        w_min = minimum(w)
        w_max = maximum(w)
        w = (w .- w_min) ./ (w_max - w_min)
    end
    return w
end

function model_average(results::Vector{uwreal}, ww)
    fin_res = sum(ww .* results)
    aux1 = sum( ww .* results.^2)
    aux2 = sum( ww .* results)^2
    syst = sqrt(abs(value(aux1 - aux2)))
    #syst = sum( )
    return fin_res, syst
    
end