using DelimitedFiles

### load in data used for each instance 
function get_data(data_file)

    # define empty matrices and vectors for data storage
    A = []
    c = []
    b = []
    
    # initialize m and n outside of the loop 
    m = 0
    n = 0

    line_counter = 1 # keep track of the number of lines

    for row in eachrow(readdlm(data_file))
        if line_counter == 1
            m = row[1]
            n = row[2]
        else
            if line_counter <= m+1
                c = push!(c, row) ## define the c matrix
            elseif line_counter > m+1 && line_counter <= 2*m+1
                A = push!(A, row) ## define the A matrix
            else 
                b = row
            end
        end 

        line_counter += 1 # increase the line count
    end

    ## convert from vector of vectors to matrices
    c = mapreduce(permutedims, vcat, c)
    A = mapreduce(permutedims, vcat, A)

    ## remove unneeded elements from b
    b = filter(e->e∉[""], b)

    return A,c,b
end