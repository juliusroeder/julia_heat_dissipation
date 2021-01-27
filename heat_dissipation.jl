
function read_data(path_images)
    temp_matrix_string = readlines(path_images)[4:end]
    conduct_matrix_string = readlines(path_images)[4:end]

    temp_matrix = zeros(Float64, row_num * col_num)
    conduct_matrix = zeros(Float64, row_num * col_num)

    j = 1

    for (temp_string, conduct_string) in zip(temp_matrix_string, conduct_matrix_string)
        temp_line = split(temp_string)
        conduct_line = split(conduct_string)
        for (t_num, c_num) in zip(temp_line, conduct_line)
            t = parse(Float64, t_num)
            c = parse(Float64, c_num)
            temp_matrix[j] = io_tmin + t * (io_tmax - io_tmin) / maxv
            conduct_matrix[j] = c
            j += 1
        end
    end
    return (temp_matrix, conduct_matrix)
end

function make_halo_smear_cols(temp_matrix, conduct_matrix)
    # make halo rows
    temp_matrix = reshape(temp_matrix, row_num, col_num)
    temp_matrix = vcat(reshape(temp_matrix[1,:], 1, col_num), temp_matrix)
    temp_matrix = vcat(temp_matrix, reshape(temp_matrix[end,:], 1, col_num))

    conduct_matrix = reshape(conduct_matrix, row_num, col_num)
    conduct_matrix = vcat(reshape(conduct_matrix[1,:], 1, col_num), conduct_matrix)
    conduct_matrix = vcat(conduct_matrix, reshape(conduct_matrix[end,:], 1, col_num))

    # make side columns
    first_col = temp_matrix[:,1]
    last_col = temp_matrix[:,end]
    temp_matrix = hcat(temp_matrix, first_col)
    temp_matrix = hcat(last_col, temp_matrix)

    first_col = conduct_matrix[:,1]
    last_col = conduct_matrix[:,end]
    conduct_matrix = hcat(conduct_matrix, first_col)
    conduct_matrix = hcat(last_col, conduct_matrix)
    return (temp_matrix, conduct_matrix)
end

macro swap!(a::Symbol,b::Symbol)
       blk = quote
         if typeof($(esc(a))) != typeof($(esc(b)))
           throw(ArgumentError("Arrays of different type"))
         else
             c = $(esc(a))
             $(esc(a)) = $(esc(b))
             $(esc(b)) = c
           end
         end
         return blk
       end

io_tmin = -100.0
io_tmax = 100.0
maxv = 655535
row_num = 150
col_num = 100
max_iter = 10

path_images = "images/pat1_100x150.pgm"

temp_matrix = zeros(Float64, row_num * col_num)
conduct_matrix = zeros(Float64, row_num * col_num)
new_temp_matrix = zeros(Float64, row_num + 2,  col_num + 2)

temp_matrix, conduct_matrix = read_data(path_images)
temp_matrix, conduct_matrix = make_halo_smear_cols(temp_matrix, conduct_matrix)


# make halo cells
new_temp_matrix[1,:] = temp_matrix[1,:]
new_temp_matrix[end,:] = temp_matrix[end,:]

for t in 1:max_iter
    for i in 2:row_num
        for j in 2:col_num
            direct_neighbors = temp_matrix[i - 1, j] + temp_matrix[i + 1, j]
                             + temp_matrix[i, j - 1] + temp_matrix[i, j + 1]
            indirect_neighbours = temp_matrix[i - 1, j - 1] + temp_matrix[i - 1, j + 1]
                                + temp_matrix[i + 1, j - 1] + temp_matrix[i + 1, j + 1]
            new_temp_matrix[i, j] = temp_matrix[i, j] + direct_neighbors + indirect_neighbours
        end
    end

    #check if done; if yes break

    #check if report

    #smear columns

    #swap new and old matrix
    swap!(new_temp_matrix, temp_matrix)
end

#report
