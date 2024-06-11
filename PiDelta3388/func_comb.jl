################################
# ISOSCALAR CHANNEL
###############################

@. dltiso_basemodel(x,p) = (x[:,3] - 1.5*x[:,2]) * (p[1] + p[2]*x[:,1])

phi42(x) = (x[:,3] - 1.5*x[:,2]).^2  # phi4 -1.5phi2
a3cutoff(x)  =  (x[:,3] - 1.5*x[:,2]) .* x[:,1].^(3/2)      # (a^2/8t0)^{3/2}
phi4sq(x) =  (x[:,3] - 1.5*x[:,2]) .* x[:,3]                # phi4


model_var_list = [phi42, a3cutoff, phi4sq]
model_var_label = ["phi4-1.5phi2", "a3", "phi4"]
model_map = [Bool.([i,j,k]) for i=0:1 for j=0:1 for k=0:1]

n_par_var = length(model_var_list) # number of extra parameters
n_par_tot_dltiso =  [2]
label_tot_dltiso = Vector{Vector{String}}(undef, 0)
push!(label_tot_dltiso, ["a2"])

f_tot_dltiso = Vector{Vector{Function}}(undef, n_par_var+1)
f_tot_dltiso[1] = Vector{Function}(undef, 1)
f_tot_dltiso[1][1] = (x,p) -> dltiso_basemodel(x,p)

for n = 2:n_par_var+1
    aux = filter(x->sum(x)==n-1, model_map)
    f_tot_dltiso[n] = Vector{Function}(undef, length(aux))
    f_aux = []
    for (k, a) in enumerate(aux)
        
        push!(n_par_tot_dltiso, n_par_tot_dltiso[1]+n-1)
        push!(label_tot_dltiso, model_var_label[a])
        f_tot_dltiso[n][k] = (x,p) -> dltiso_basemodel(x,p) .+ sum([p[i+n_par_tot_dltiso[1]] for i=1:(n-1)] .* (fill(x, n-1) .|> model_var_list[a]))
        # push!(f_aux,  (x,p) -> isov_basemodel(x,p) .+ sum([p[i+n_par_tot_isov[1]] for i=1:(n-1)] .* (fill(x, n-1) .|> model_var_list[a])))
    end
    f_tot_dltiso[n] = [f_tot_dltiso[n][l] for l in eachindex(f_tot_dltiso[n]) if isassigned(f_tot_dltiso[n], l)]
end

f_tot_dltiso = vcat(f_tot_dltiso...)
