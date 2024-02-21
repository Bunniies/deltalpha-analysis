function read_BDIO(path::String, uinfo::Int64)
    r = Vector{uwreal}(undef, 0)
    fb = BDIO_open(path, "r")
    BDIO_seek!(fb)
    
    if BDIO_get_uinfo(fb) == uinfo
        push!(r, read_uwreal(fb))
    end

    while BDIO_seek!(fb, 2) 
        if BDIO_get_uinfo(fb) == uinfo
            push!(r, read_uwreal(fb))
        end
    end

    BDIO_close!(fb)
    return r
end

function read_BDIO(path::String, type::String, obs::String)
    dict_dalpha = Dict(
        "t0"        => 0,
        "fvc"       => 1,
        "pi_33_ll"  => 2,
        "pi_33_lc"  => 3,
        "pi_88_ll_conn" => 4,
        "pi_88_lc_conn" => 5,
        "pi_08_ll_conn" => 6,
        "pi_08_lc_conn" => 7
    )
    dict_spectrum = Dict(
        "t0"     => 0,
        "mpi"    => 1,
        "mk"     => 2
    )
    dict_tree_level = Dict(
        "3l_33_ll" => 0,
        "3l_33_lc" => 1
    )

    dict2dict = Dict("dalpha" => dict_dalpha, "spectrum" => dict_spectrum, "3level" => dict_tree_level)
    
    if !(type in keys(dict2dict))
        error("Incorrect type.\ntype = $(keys(dict2dict))")
    end
    if !(obs in keys(dict2dict[type]))
        error("Incorrect obs.\nobs = $(keys(dict2dict[type]))")
    end
    return read_BDIO(path, dict2dict[type][obs])
end
