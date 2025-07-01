function read_AdlerPy_by_component(path, name; N=400)
    fname = readdlm(joinpath(path,name), comments=true)
    #return fname
    qval = fname[:,2]
    adl_val = fname[:,3]
    adl_err = fname[:,4]
    adc_val = fname[:,5]
    adc_err = fname[:,6]
    adb_val = fname[:,7]
    adb_err = fname[:,8]
    adozi_val = fname[:,9]
    adozi_err = fname[:,10]

    adl, adc, adb, adozi = [Vector{uwreal}(undef, N) for k in 1:4]

    for k in eachindex(qval)
        adl[k] = uwreal([adl_val[k], adl_err[k]], "adl")
        adc[k] = uwreal([adc_val[k], adc_err[k]], "adc")
        adb[k] = uwreal([adb_val[k], adb_err[k]], "adb")
        adozi[k] = uwreal([adozi_val[k], adozi_err[k]], "adozi")
    end
    return qval, adl, adc, adb,  adozi
end

function read_std_trunc_error_data(path, name; N=400)
    qval, adl, adc, adb, adozi = read_AdlerPy_by_component(path, name, N=N)
    name_4l = split(name, "_")[1] *"_3loops.dat"
    _, adl_4l, adc_4l, adb_4l, adozi_4l = read_AdlerPy_by_component(path, name_4l, N=N)

    uwerr.(adl); uwerr.(adc); uwerr.(adb); uwerr.(adozi)
    uwerr.(adl_4l); uwerr.(adc_4l); uwerr.(adb_4l); uwerr.(adozi_4l)

    err_l = value.(adl) .- value.(adl_4l)
    err_c = value.(adc) .- value.(adc_4l)
    err_b = value.(adb) .- value.(adb_4l)
    err_ozi = value.(adozi) .- value.(adozi_4l)

    for k in eachindex(adl)
        adl[k] += uwreal([0.0, err_l[k]], "err")
        adc[k] += uwreal([0.0, err_c[k]], "err")
        adb[k] += uwreal([0.0, err_b[k]], "err")
        adozi[k] += uwreal([0.0, err_ozi[k]], "err")
    end

    return qval, adl, adc, adb, adozi
end

function read_pQCD1_datfile(path; N=400)
    fname = readdlm(path, comments=true)
    println(path)
    qval= -fname[1:N,1]
    mean_val = fname[1:N,3]

    path_err = joinpath(splitpath(path)[1:end-1])
    fname_pl = readdlm(joinpath(path_err,"data_err_alpha/pQCDAdler.dap52"), comments=true)
    fname_min = readdlm(joinpath(path_err,"data_err_alpha/pQCDAdler.dam52") , comments=true)
    err_input = (fname_pl[1:N,2] .- fname_min[1:N,2]) ./ 2

    fname_pl = readdlm(joinpath(path_err,"data_err_ope/pQCDAdler.dat52p"), comments=true)
    fname_min = readdlm(joinpath(path_err,"data_err_ope/pQCDAdler.dat52m") , comments=true)
    err_ope = (fname_pl[1:N,3] .- fname_min[1:N,3]) ./ 2
    
    tot_err = sqrt.(err_input.^2 + err_ope.^2)
    
    ad_res = []
    for k in eachindex(mean_val)
        push!(ad_res, uwreal([mean_val[k], tot_err[k]], "Dfunc pQCD1"))
    end
    return reverse(qval), reverse(ad_res)
end


function adler_lattice(Q)
    mcov = [1*0.0023^2 0.455*0.0023*0.015 0.17*0.0023*0.0006 0.641*0.0023*0.22 0.351*0.0023*0.19 0.0489*0.0023*0.0012;
        0.455*0.0023*0.015 1*0.015^2 0.823*0.015*0.0006 0.946*0.015*0.22 0.977*0.015*0.19 -0.0934*0.015*0.0012;
        0.17*0.0023*0.0006 0.823*0.0006*0.015 1*0.0006^2 0.642*0.0006*0.22 0.915*0.0006*0.19 0.0667*0.0006*0.0012;
        0.641*0.0023*0.22 0.946*0.22*0.015 0.642*0.22*0.0006 1*0.22^2 0.869*0.22*0.19 -0.044*0.22*0.0012;
        0.351*0.0023*0.19 0.977*0.19*0.015 0.915*0.19*0.0006 0.869*0.19*0.22 1*0.19^2 -0.115*0.19*0.0012;
        0.0489*0.0023*0.0012 -0.0934*0.0012*0.015 0.0667*0.0012*0.0006 -0.044*0.0012*0.22 -0.115*0.0012*0.19 1*0.0012^2
    ]
    mcov = round.(mcov, digits=18)
    # return mcov
    ave = [0.1094, 0.093, 0.0039, 2.85, 1.03, 0.0166]
    param = cobs(ave, mcov, "2022_param")
    Q2 = Q*Q
    return (12*pi^2*Q2*(-(param[1]*(-1 + param[5]*Q2^2 + 2*param[6]*Q2^3)) + 
           Q2*(param[3]*Q2*(3 + 2*param[4]*Q2 + param[5]*Q2^2) + param[2]*(2 + param[4]*Q2 - param[6]*Q2^3))))/(1 + param[4]*Q2 + param[5]*Q2^2 + param[6]*Q2^3)^2

end




function get_integral(data, dl) # apply trapezodial rule assuming data to be equally dl distant 
    return sum((data[1:end-1] + data[2:end])*0.5*dl)
end