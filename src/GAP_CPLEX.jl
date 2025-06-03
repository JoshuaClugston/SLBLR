using JuMP, CPLEX, DelimitedFiles
include("./get_data.jl")


### define matrices and vectors corresponding to data
# instance 1
A_1, c_1, b_1 = get_data("./data/gap_e/e201600")

# instance 2
A_2, c_2, b_2 = get_data("../data/gap_e/e801600")

# # instance 3 
A_3, c_3, b_3 = get_data("../data/gap_d/d05100")


### define a function which describes the model, taking as input the data 
### Function will return the optimal objective value. 
### Remaining output will be written to an external file, which is accessible with all submission documents.

function solve_model(A,c,b)
    model = Model(CPLEX.Optimizer) # construct and empty model using CPLEX as the solver
    set_attribute(model, "CPX_PARAM_EPAGAP", 1e-8)
    set_attribute(model, "CPX_PARAM_EPGAP", 1e-8)
    set_attribute(model, "CPX_PARAM_MIPDISPLAY", 5)
     
    # get dimensions
    J = size(A,2) # number of columns
    I = size(A,1) # number of rows 
    
    # define variables
    @variable(model, x[i=1:I,j=1:J], Bin)

    # define constraints
    @constraint(model, [i=1:I], sum(A[i,j]*x[i,j] for j in 1:J) <= b[i])
    @constraint(model, [j=1:J], sum(x[i,j] for i in 1:I) == 1)

    # define objective function 
    @objective(model, Min, sum( sum(c[i,j]*x[i,j] for j in 1:J) for i in 1:I))

    optimize!(model)

    return [objective_value(model), value.(x)]  # return the solution
end

# println("Instance one optimal objective value = ", solve_model(A_1,c_1,b_1))
# println("Instance two optimal objective value = ", solve_model(A_2,c_2,b_2))
println("Instance three optimal objective value = ", solve_model(A_3,c_3,b_3)[1])
