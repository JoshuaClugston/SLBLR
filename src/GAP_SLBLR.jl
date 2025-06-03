# using Pkg; Pkg.add("SCIP")
using ArgParse, JuMP, CPLEX, Distributions, LinearAlgebra

include("./get_data.jl") ## include function which reads in data file to matrices

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin 
        "--data_file"
            arg_type = String 
            default  = "./gap_d/d05100" ## ./gap_e/e201600, ./gap_e/e801600, ./gap_d/d05100
        "--max_iter"
            arg_type = Int64
            default  = 1000
        "--time_limit"
            arg_type = Float64
            default  = 600 ## seconds 
        "--nu"
            arg_type = Float64
            default  = 2.0
        "--zeta"
            arg_type = Float64
            default = 1/1.5
        "--epsilon_gap"
            arg_type = Float64
            default  = 0.02 ### 
        "--epsilon"
            arg_type = Float64
            default  = 0.01 
        "--step_init" 
            arg_type = Float64
            default = 0.02 
    end
    return parse_args(s)
end
args = parse_commandline() 

A,c,b = get_data(args["data_file"]) 

m = size(A,1)
n = size(A,2)

### set parameters
# λ0     = rand(Uniform(10,20),m) ## each λ0_i in [90,100] i = 1,...,m
γ      = 1/m
ν      = args["nu"]
ϵ_gap  = args["epsilon_gap"] 
ϵ      = args["epsilon"]
s0     = args["step_init"]
ζ      = args["zeta"]
time_limit = args["time_limit"]
max_iter = args["max_iter"]
qmax = -1e9
num_subproblems = 20 ## number of subproblems to consider for heuristics

###############################################################################################################################################################

###############################################################################################################################################################

function get_q0(A,b,c,s0)
    model = Model(CPLEX.Optimizer)
    set_silent(model)
    m,n = size(A)
    c1 = Dict()

    @variable(model, x[1:m, 1:n] >=0)

    @constraint(model, [i=1:m], sum( A[i,j]*x[i,j] for j in 1:n ) <= b[i]) 

    for j = 1:n 
        c1[j] = @constraint(model, sum( x[i,j] for i in 1:m) == 1)
    end

    @objective(model, Min, sum( sum(c[i,j]*x[i,j] for i in 1:m ) for j in 1:n))

    optimize!(model)

    x_f = value.(x)

    λ0 = dual.(values(c1))

    gx0 = Float64[]
    for j in 1:n
        push!(gx0, sum(x_f[i,j] for i in 1:m) - 1)
    end

    Lx0λ0 = sum(sum(  c[i,j] * x_f[i,j] for j in 1:n) for i in 1:m) + sum(λ0[j] * gx0[j] for j in 1:n)

    q0 = Lx0λ0 + s0 * norm(gx0)^2

    return q0, λ0 
end

q0, λ0 = get_q0(A,b,c, s0) 

######

### main algorithm 
function SLBLR(λ0, s0, γ, ν, ζ, q0, ϵ_gap, ϵ, max_iter, time_limit, qmax, A, b, c, num_subproblems, args) # ζ = 1/1.5, ν = 2, s0 = 0.02, γ = 1/m
    m,n = size(A)
    ck = Dict() ## for first iteration 
    value = true
    iter_counter = 0 ## initialize the iteration counter 
    k = 1
    i = 1
    j = 1
    λk = λ0
    sk = s0
    qj = q0 
    x̃k, obj = solve_relaxed_problem(A,b,c,λk,true) ## get initial x
    qλk = obj
    gap_progress = [] ## empty array for the gap progress
    λ = [] ## empty array for multiplier
    args["elapsed_time"] = 0
    timer = time() ## initialize the total time limit 
    timer1 = time()
    while true
        ### minimize the lagrangian subproblem ensure that the sufficient decrease condition is satisfied before accepting the solution
        if k != 1
            # println("--------------------- SOLVING SUBPROLEM AT ITERATION $(k) ------------------------------------------- \n\n")
            x̃k_new, value = solve_subproblem(A,b,c, λk, x̃k, i)
            x̃k = x̃k_new
        end

        gx̃ = [] # empty array for subgradient values at this iteration
        for j in 1:n
            push!(gx̃, sum(x̃k[i,j] for i in 1:m) - 1) 
        end 
        
        Lx̃λ̃ = sum(sum(  c[i,j] * x̃k[i,j] for j in 1:n) for i in 1:m) + sum(λk[j] * gx̃[j] for j in 1:n) ## get the langrangian
        println(Lx̃λ̃)
 
        sk = ζ * γ * (qj - Lx̃λ̃) / norm(gx̃)^2 # update the stepsize; set minimum stepsize to s0. 

        λk_new = λk + sk * gx̃ ## update the multipliers 

        qmax = max(qmax, sk * norm(gx̃)^2 / γ + Lx̃λ̃)

        # update indices
        if i == m # corresponding to the number of 1:m (rows)
            i = 1
        else 
            i += 1
        end
        k += 1 # number of iterations 

        ## solve the feasibility problem using CPLEX
        # println("\n\n --------------------------- SOLVING CONSTRAINT SATISFIABILITY PROBLEM AT ITERATION $(k) ---------------------------\n ")
        ret = constraint_satisfiability_problem(λk, λk_new, m, k, ck, ν, sk, gx̃)

        if ret == INFEASIBLE ## check if constraint problem is infeasible 
            qj = qmax
            println("\n\n Constraint satisfiability problem is infeasible at iteration $(k). Updating qj = $(qj), j = $(j+1) ----------------------------\n\n")
            qmax = -1e9 ## really small number
            j += 1
            empty!(ck) ## reset the constraint set for the feasibility problem 
        else
            ck = ret
        end

        λk = λk_new ## update λk 

        push!(λ, λk)
        push!(gap_progress, (qj - Lx̃λ̃)/qj)

        # println("---------------------- THIS IS THE qλk VALUE: ", qλk, ". WHILE THE qj value is $(qj) -------------------------------------------------- \n\n")

        # if (qj - qλk)/qj < ϵ || iter_counter >= max_iter || args["elapsed_time"] >= time_limit 
        if (qj - Lx̃λ̃)/qj < ϵ || iter_counter >= max_iter || args["elapsed_time"] >= time_limit 
            println("SOLVING FEASIBILITY MILP AT ITERATION $(k)")
            x̂, fx̂ = feasibility_milp(A,b, x̃k, num_subproblems) # feasibility problem for original MILP with heuristics -- num_subproblems set outside of function described the number of subproblem solutions to consider
            x̃k = x̂

            if (fx̂ - Lx̃λ̃)/fx̂ <= ϵ_gap
                println("feasible objective value found after $(time() - timer1)s: $(fx̂). Gap reported: $((fx̂ - Lx̃λ̃)/fx̂)") ### report the amount of time needed, number of iterations, and objective function value at optimality
                return gap_progress, λ
            else
                println("feasible objective value found: $(fx̂). Not optimal after $(args["elapsed_time"])s and $(k) iterations. Continuing to solve.")
            end 
            if args["elapsed_time"] >= time_limit || iter_counter >= max_iter ### terminate if the total alloted time or iterations have been exceeded
                # println("feasible objective value found: $(fx̂). Not optimal after $(args["elapsed_time"])s and $(k) iterations. Terminate with gap = $((fx̂ - qλk)/fx̂).")
                println("feasible objective value found: $(fx̂). Not optimal after $(args["elapsed_time"])s and $(k) iterations. Terminate with gap = $((fx̂ - Lx̃λ̃)/fx̂).")
                return gap_progress, λ
            end
        end

        # qλk = solve_relaxed_problem(A,b,c,λk, false) 
        
        args["elapsed_time"] = time() - timer ## update the main timer 
        println("\n\n ------------------- Current elasped time at iteration $(k) is: $(args["elapsed_time"])s -----------------------------")
        # println(" ------------------- Current qλk value at iteration $(k) is: $(qλk) --------------------------------------------------")
        # println(" ------------------- Current qj value at iteration $(k) is: $(qj) ---------------------------------------------------- ")
        # println(" ------------------- Progress of (qj - qλk)/qj at iteration $(k) is: $((qj - qλk)/qj) -------------------------------- \n\n")
        iter_counter += 1 ## update the iteration count 
    end
end

## Solve the subproblem for a particular ik using sufficient decrease condition
function solve_subproblem(A, b, c, λk, xk, ik)
    m, n = size(A)
    model = Model(CPLEX.Optimizer)
    # set_attribute(model, "CPXPARAM_MIP_Display", 5)
    # set_attribute(model, "CPXPARAM_Simplex_Display", 2)
    set_silent(model)
    @variable(model, x[1:m, 1:n], Bin)

    @constraint(model, sum(A[ik, j] * x[ik, j] for j in 1:n) <= b[ik])

    @objective(model, Min, sum((c[ik, j] + λk[j]) * x[ik, j] for j in 1:n))

    optimize!(model)
    xk_1 = value.(x)

    ### check the surrogate optimality condition 
    sum_lb = sum((c[i, j] + λk[j]) * xk_1[i, j] for i in 1:m, j in 1:n)
    sum_up = sum((c[i, j] + λk[j]) * xk[i, j] for i in 1:m, j in 1:n) 

    if sum_lb < sum_up
        x_out = xk_1
        return x_out, true
    else
        return xk, false
    end

end

## want to get record of constraints, and add new constraints on top of already used constraints ... 
function constraint_satisfiability_problem(λk, λk_new, n, k, ck, ν, sk, gx̃) # need to account for number of times this problem is infeasible (j), ... goal is to get an infeasible solution ...
# function constraint_satisfiability_problem(λk, m, k, ck, sk, gx̃)
    model = Model(CPLEX.Optimizer) 
    set_silent(model)

    @variable(model, λ[1:n])

    ck[k] = @constraint(model, sum( (λ[j] - λk_new[j])^2 for j in 1:n) <= sum( (λ[j] - λk[j])^2 for j in 1:n) )
    # ck[k] = @constraint(model, sum( (λ[j] - λk_new[j])^2 for j in 1:n) <= (1-2*ν*sk)*sum( (λ[j] - λk[j])^2 for j in 1:n) ) # constant multiplied to help speed up convergence 
    # ck[k] = @constraint(model, 2* sum((λ[i] - λk[i]) * gx̃[i] for i in 1:n) >= sk * norm(gx̃)^2) # linearized constraint
    optimize!(model)

    if termination_status(model) == INFEASIBLE
        return termination_status(model)
    elseif termination_status(model) == OPTIMAL
        return ck 
    end
end 


function feasibility_milp(A,b, xk, num_subproblems) 
    model = Model(CPLEX.Optimizer)
    set_silent(model)
    m, n = size(A)

    @variable(model, x[i=1:m,j=1:n], start = xk[i,j], Bin)
    # @variable(model, z[i=1:m], Bin) ## for heuristic 
    
    @constraint(model, [i=1:m], sum( A[i,j]*x[i,j] for j in 1:n) <= b[i])
    @constraint(model, [j=1:n], sum( x[i,j] for i in 1:m) == 1)

    ## heuristics to perturb solution obtained from previous iterations 
    # @constraint(model, sum(z[i] for i in 1:m) <= num_subproblems)
    # @constraint(model, [i=1:m, j=1:n], x[i,j] - xk[i,j] <= z[i])
    # @constraint(model, [i=1:m, j=1:n], x[i,j] - xk[i,j] >= -z[i])

    @objective(model, Min, sum( sum(c[i,j]*x[i,j] for i in 1:m ) for j in 1:n))

    optimize!(model)

    x_out = value.(x)
    obj_out = objective_value(model)

    return x_out, obj_out
end

function solve_relaxed_problem(A,b,c, λk, init)
    model = Model(CPLEX.Optimizer)
    set_silent(model) 
    m, n = size(A)

    @variable(model, x[i=1:m,j=1:n], Bin)
    # @variable(model, 0 <= x[i=1:m,j=1:n] <= 1)

    @constraint(model, [i=1:m], sum(A[i,j]*x[i,j] for j in 1:n) <= b[i])

    @objective(model, Min, sum(  (c[i,j] + λk[j])*x[i,j] for j in 1:n, i in 1:m) - sum(λk[j] for j in 1:n))

    optimize!(model)

    if init == true 
        x_out = value.(x)
        obj = objective_value(model)
        return x_out, obj
    else
        return objective_value(model)
    end
end

gap_progress, λ = SLBLR(λ0, s0, γ, ν, ζ, q0, ϵ_gap, ϵ, max_iter, time_limit, qmax, A, b, c, num_subproblems, args)

println(gap_progress)
# println(λ)

