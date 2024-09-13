const Zvc_l = Dict(
"A653" => uwreal([1.32281, 0.00072], "Zvc"),
"A654" => uwreal([1.30495, 0.00106], "Zvc"),

"H101" => uwreal([1.20324, 0.00071], "Zvc"),
"H102" => uwreal([1.19743, 0.00100], "Zvc"),
"H105" => uwreal([1.18964, 0.00075], "Zvc"),
"N101" => uwreal([1.18964, 0.00075], "Zvc"),
"C101" => uwreal([1.18500, 0.00044], "Zvc"),
"C102" => uwreal([1.18500, 0.00044], "Zvc"), # TO BE UPDATED!
"D150" => uwreal([1.18500, 0.00044], "Zvc"), # TO BE UPDATED!

"B450" => uwreal([1.12972, 0.00083], "Zvc"),
"S400" => uwreal([1.11159, 0.00089], "Zvc"),
"N451" => uwreal([1.11412, 0.00058], "Zvc"),
"D450" => uwreal([1.10790, 0.00026], "Zvc"),
"D451" => uwreal([1.10790, 0.00026], "Zvc"), # TO BE UPDATED!
"D452" => uwreal([1.10790, 0.00023], "Zvc"),

"N202" => uwreal([1.04843, 0.00085], "Zvc"), 
"N203" => uwreal([1.04534, 0.00039], "Zvc"),
"N200" => uwreal([1.04012, 0.00025], "Zvc"),
"D251" => uwreal([1.04012, 0.00045], "Zvc"), # TO BE UPDATED! 
"D200" => uwreal([1.03587, 0.00022], "Zvc"),
"D201" => uwreal([1.03587, 0.00022], "Zvc"), # TO BE UPDATED!
"E250" => uwreal([1.03310, 0.00011], "Zvc"),

"J307" => uwreal([0.96037, 0.00039], "Zvc"), # TO BE UPDATED!
"N300" => uwreal([0.97722, 0.00060], "Zvc"),
"N302" => uwreal([0.97241, 0.00030], "Zvc"),
"J306" => uwreal([0.96037, 0.00039], "Zvc"), # TO BE UPDATED!
"J303" => uwreal([0.96037, 0.00039], "Zvc"),
"J304" => uwreal([0.96037, 0.00039], "Zvc"), # TO BE UPDATED!
"E300" => uwreal([0.96639, 0.00026], "Zvc"),
"F300" => uwreal([0.96639, 0.00026], "Zvc"), # TO BE UPDATED!

"J500" => uwreal([0.93412, 0.00051], "Zvc"),
"J501" => uwreal([0.93412, 0.00051], "Zvc") # TO BE UPDATED!
)

m_ens = Dict(
    # pion and kaon masses extracted from 2206.06582v2
    # rho masses extracted from 2203.08676v2
    "A653" => Dict("m_pi" => uwreal([0.21193,0.00091], "mpi"), "m_K" => uwreal([0.21193,0.00091], "mK"), "m_rho" => uwreal([NaN, NaN],    "mrho")),
    "A654" => Dict("m_pi" => uwreal([0.16647,0.00121], "mpi"), "m_K" => uwreal([0.22712,0.00089], "mK"), "m_rho" => uwreal([NaN, NaN],    "mrho")),
    
    "H101" => Dict("m_pi" => uwreal([0.18217,0.00062], "mpi"), "m_K" => uwreal([0.18217,0.00062], "mK"), "m_rho" => uwreal([0.375,0.002], "mrho")),
    "H102" => Dict("m_pi" => uwreal([0.15395,0.00071], "mpi"), "m_K" => uwreal([0.19144,0.00057], "mK"), "m_rho" => uwreal([0.358,0.003], "mrho")),
    "H105" => Dict("m_pi" => uwreal([0.12136,0.00124], "mpi"), "m_K" => uwreal([0.20230,0.00061], "mK"), "m_rho" => uwreal([0.338,0.011], "mrho")),
    "N101" => Dict("m_pi" => uwreal([0.12150,0.00055], "mpi"), "m_K" => uwreal([0.20158,0.00031], "mK"), "m_rho" => uwreal([0.340,0.004], "mrho")),
    "C101" => Dict("m_pi" => uwreal([0.09569,0.00073], "mpi"), "m_K" => uwreal([0.20579,0.00034], "mK"), "m_rho" => uwreal([0.326,0.003], "mrho")),
    "C102" => Dict("m_pi" => uwreal([0.09671,0.00078], "mpi"), "m_K" => uwreal([0.21761,0.00047], "mK"), "m_rho" => uwreal([0.323,0.004], "mrho")), # rho updated
    "D150" => Dict("m_pi" => uwreal([0.05654,0.00094], "mpi"), "m_K" => uwreal([0.20835,0.00035], "mK"), "m_rho" => uwreal([0.326,0.003], "mrho")), # rho not upadted
    
    "B450" => Dict("m_pi" => uwreal([0.16063,0.00045], "mpi"), "m_K" => uwreal([0.16063,0.00045], "mK"), "m_rho" => uwreal([0.337,0.001], "mrho")),
    "S400" => Dict("m_pi" => uwreal([0.13506,0.00044], "mpi"), "m_K" => uwreal([0.17022,0.00039], "mK"), "m_rho" => uwreal([0.312,0.004], "mrho")),
    "N451" => Dict("m_pi" => uwreal([0.11072,0.00029], "mpi"), "m_K" => uwreal([0.17824,0.00018], "mK"), "m_rho" => uwreal([0.302,0.004], "mrho")),
    "D450" => Dict("m_pi" => uwreal([0.08329,0.00043], "mpi"), "m_K" => uwreal([0.18384,0.00018], "mK"), "m_rho" => uwreal([0.303,0.008], "mrho")),
    "D451" => Dict("m_pi" => uwreal([0.08359,0.00030], "mpi"), "m_K" => uwreal([0.19402,0.00014], "mK"), "m_rho" => uwreal([0.303,0.008], "mrho")), # rho not updated
    "D452" => Dict("m_pi" => uwreal([0.05941,0.00055], "mpi"), "m_K" => uwreal([0.18651,0.00015], "mK"), "m_rho" => uwreal([0.297,0.008], "mrho")), # rho not updated
    
    "H200" => Dict("m_pi" => uwreal([0.13535,0.00060], "mpi"), "m_K" => uwreal([0.13535,0.00060], "mK"), "m_rho" => uwreal([0.286,0.003], "mrho")),
    "N202" => Dict("m_pi" => uwreal([0.13424,0.00031], "mpi"), "m_K" => uwreal([0.13424,0.00031], "mK"), "m_rho" => uwreal([0.280,0.003], "mrho")),
    "N203" => Dict("m_pi" => uwreal([0.11254,0.00024], "mpi"), "m_K" => uwreal([0.14402,0.00020], "mK"), "m_rho" => uwreal([0.268,0.001], "mrho")),
    "N200" => Dict("m_pi" => uwreal([0.09234,0.00031], "mpi"), "m_K" => uwreal([0.15071,0.00023], "mK"), "m_rho" => uwreal([0.252,0.002], "mrho")),
    "D251" => Dict("m_pi" => uwreal([0.09293,0.00016], "mpi"), "m_K" => uwreal([0.15041,0.00012], "mK"), "m_rho" => uwreal([0.250,0.002], "mrho")), # rho not updated
    "D200" => Dict("m_pi" => uwreal([0.06507,0.00028], "mpi"), "m_K" => uwreal([0.15630,0.00015], "mK"), "m_rho" => uwreal([0.250,0.002], "mrho")),
    "D201" => Dict("m_pi" => uwreal([0.06499,0.00043], "mpi"), "m_K" => uwreal([0.16309,0.00024], "mK"), "m_rho" => uwreal([0.239,0.002], "mrho")), # rho updated
    "E250" => Dict("m_pi" => uwreal([0.04170,0.00041], "mpi"), "m_K" => uwreal([0.15924,0.00009], "mK"), "m_rho" => uwreal([0.251,0.004], "mrho")),
    
    "N300" => Dict("m_pi" => uwreal([0.10569,0.00023], "mpi"), "m_K" => uwreal([0.10569,0.00023], "mK"), "m_rho" => uwreal([0.222,0.003], "mrho")),
    "N302" => Dict("m_pi" => uwreal([0.08690,0.00034], "mpi"), "m_K" => uwreal([0.11358,0.00028], "mK"), "m_rho" => uwreal([0.216,0.003], "mrho")),
    "J303" => Dict("m_pi" => uwreal([0.06475,0.00018], "mpi"), "m_K" => uwreal([0.11963,0.00016], "mK"), "m_rho" => uwreal([0.200,0.002], "mrho")),
    "J304" => Dict("m_pi" => uwreal([0.06550,0.00020], "mpi"), "m_K" => uwreal([0.13187,0.00017], "mK"), "m_rho" => uwreal([0.200,0.002], "mrho")), # rho  updated
    "E300" => Dict("m_pi" => uwreal([0.04393,0.00016], "mpi"), "m_K" => uwreal([0.12372,0.00010], "mK"), "m_rho" => uwreal([0.198,0.002], "mrho")),
    
    "J500" => Dict("m_pi" => uwreal([0.08153,0.00019], "mpi"), "m_K" => uwreal([0.08153,0.00019], "mK"), "m_rho" => uwreal([0.173,0.001], "mrho")), # rho updated
    "J501" => Dict("m_pi" => uwreal([0.06582,0.00023], "mpi"), "m_K" => uwreal([0.08794,0.00022], "mK"), "m_rho" => uwreal([0.167,0.002], "mrho"))  # rho updated
)

const b_values = [3.34, 3.40, 3.46, 3.55, 3.70, 3.85]
const hc = 197.3269804 #MeV fm

# scale setting
const t0sqrt_ph = uwreal([0.1449, 1e-8], "sqrtt0 [fm]")  # Regensburg with artificial error
const t0sqrt_ph_err = uwreal([0.1449, 0.007], "sqrtt0 [fm] err")  # Regensburg with correct error

# const t0sqrt_ph_noerr = uwreal([0.1449, 0.0000], "sqrtt0 [fm]")  # Regensburg without error
# println("scale setting from BRUNO ET AL!!!")
const t0sqrt_ph_bruno = uwreal([0.1467, 0.0017], "sqrtt0 Bruno [fm]") 

#2211.03744
const t0_data = [2.204, 2.872, 3.682, 5.162, 8.613, 14.011]
const t0_error = [5, 10, 12, 16, 25, 39] .* 1e-3

const Ct0 = zeros(6, 6)
for i = 1:6
    Ct0[i,i] = t0_error[i] ^ 2    
end

const t0_ = cobs(t0_data, Ct0, "t0sym/a2")
const a_ = t0sqrt_ph ./ sqrt.( t0_)

t0(beta::Float64) = t0_[b_values .== beta][1]
a(beta::Float64)  = a_[b_values .== beta][1]

