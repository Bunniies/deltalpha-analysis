using ADerrors

mean_val = 1896.11
err_lambda = mean_val - 1897.88
err_trunc = mean_val - 1894.11
err_tot = sqrt(err_lambda^2 + err_trunc^2)
b33_from_wolfram = uwreal([mean_val, err_tot], "b33" ) 

Qm2 = 9 
Q2 = [0.05, 0.1, 0.4, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0]

for q2 in Q2
    tmp = q2 / (4* Qm2) * b33_from_wolfram
    e_lambda = q2 / (4* Qm2) * err_lambda
    e_trunc = q2 / (4* Qm2) * err_trunc
    uwerr(tmp)
    println("at Q2=$(q2): $(value(tmp)) +- $(e_lambda) +- $(e_trunc)")
end