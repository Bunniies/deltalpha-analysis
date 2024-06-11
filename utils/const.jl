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

"N300" => uwreal([0.97722, 0.00060], "Zvc"),
"N302" => uwreal([0.97241, 0.00030], "Zvc"),
"J303" => uwreal([0.96037, 0.00039], "Zvc"),
"J304" => uwreal([0.96037, 0.00039], "Zvc"), # TO BE UPDATED!
"E300" => uwreal([0.96639, 0.00026], "Zvc"),

"J500" => uwreal([0.93412, 0.00051], "Zvc"),
"J501" => uwreal([0.93412, 0.00051], "Zvc") # TO BE UPDATED!
)

const b_values = [3.34, 3.40, 3.46, 3.55, 3.70, 3.85]
const hc = 197.3269804 #MeV fm

# scale setting
# const t0sqrt_ph = uwreal([0.1443, 0.0007], "sqrtt0 [fm]")  # Regensburg
println("scale setting from BRUNO ET AL!!!")
const t0sqrt_ph = uwreal([0.1467, 0.0017], "sqrtt0 [fm]") 

#1608.08900
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

