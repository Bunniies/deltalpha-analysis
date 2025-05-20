################################
# ISOVECTOR CHANNEL
###############################
charm_basemodel(x,p) = p[1] .+ p[2] .* x[:,1]  .+ p[3] .* (x[:,2] .- value.(phi2_ph))  #.+ p[4] .* (x[:,3] .- value.(phi4_ph))

# phi2lin(x)   = (x[:,2] .- value.(phi2_ph))                              # Φ2 - Φ2^{ph}
# phi4lin(x)   = (x[:,3] .- value.(phi4_ph))                              # Φ4 - Φ4^{ph}
a3cutoff(x)  = x[:,1].^(3/2)                                            # (a^2/8t0)^{3/2}
a4cutoff(x)  = x[:,1].^(2)                                              # (a^2/8t0)^{2}
a2phi2(x)    = x[:,1] .* (x[:,2] .- value(phi2_ph))                     # (a^2/8t0)*(ϕ2 - ϕ2^{ph})     
a3phi2(x)    = x[:,1].^(3/2) .* (x[:,2] .- value.(phi2_ph))              # (a^2/8t0)^{3/2}*(ϕ2 - ϕ2^{ph})    
phi2sqr(x)   = (x[:,2].^2 .- value.(phi2_ph).^2)                          # (ϕ2^2 - (ϕ2^{ph})^2) 
phi2log(x)   = x[:,2] .* log.(x[:,2]) .- value.(phi2_ph) .* log.(value.(phi2_ph))  #(ϕ2log(ϕ2) - ϕ2^{ph}log(ϕ2^{ph}))

model_var_list  = [ a3cutoff, a4cutoff, a2phi2, a3phi2, phi2sqr, phi2log]
model_var_label = [ "a3", "a4", "a2phi2", "a3phi2", "phi2sqr", "phi2log"]
model_map = [Bool.([i,j,k,l,m,n]) for i=0:1 for j=0:1 for k=0:1 for l=0:1 for m=0:1 for n=0:1]

# model_var_list  = [a3cutoff, a2phi2, a2phi4, phi2sqr,  a3phi2, phi2log]
# model_var_label = ["a3", "a2phi2", "a2phi4", "phi2sqr",  "a3phi2", "phi2log"]
# model_map = [Bool.([i,j,k,l,m,n]) for i=0:1 for j=0:1 for k=0:1 for l=0:1 for m=0:1 for n=0:1] 

n_par_var = length(model_var_list) # number of extra parameters
n_par_tot_charm =  [3]
label_tot_charm = Vector{Vector{String}}(undef, 0)
push!(label_tot_charm, ["a2"])

f_tot_charm = Vector{Vector{Function}}(undef, n_par_var+1)
f_tot_charm[1] = Vector{Function}(undef, 1)
f_tot_charm[1][1] = (x,p) -> charm_basemodel(x,p)

for n = 2:n_par_var+1
    aux = filter(x->sum(x)==n-1, model_map)
    f_tot_charm[n] = Vector{Function}(undef, length(aux))
    f_aux = []
    for (k, a) in enumerate(aux)
        if "a4" ∈ model_var_label[a] && "a3" ∉ model_var_label[a]
            continue
        end
        if "a3phi2" ∈ model_var_label[a] #&& "a2phi2" ∉ model_var_label[a]
            continue
        end
        if "phi2log" ∈ model_var_label[a] #|| "phi2sqr" ∈ model_var_label[a]
            continue
        end
        if [ "a3", "a4", "a2phi2", "a3phi2", "phi2sqr", "phi2log"] == model_var_label[a]
            continue
        end
        push!(n_par_tot_charm, n_par_tot_charm[1]+n-1)
        push!(label_tot_charm, model_var_label[a])
        f_tot_charm[n][k] = (x,p) -> charm_basemodel(x,p) .+ sum([p[i+n_par_tot_charm[1]] for i=1:(n-1)] .* (fill(x, n-1) .|> model_var_list[a]))
        # push!(f_aux,  (x,p) -> isov_basemodel(x,p) .+ sum([p[i+n_par_tot_isov[1]] for i=1:(n-1)] .* (fill(x, n-1) .|> model_var_list[a])))
    end
    f_tot_charm[n] = [f_tot_charm[n][l] for l in eachindex(f_tot_charm[n]) if isassigned(f_tot_charm[n], l)]
end

f_tot_charm = vcat(f_tot_charm...)