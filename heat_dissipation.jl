using Printf

function read_data(path_images)
    io_tmin = 0
    io_tmax = 100.0

    temp_matrix_string = readlines(path_images)
    conduct_matrix_string = readlines(path_images)

    col_num = parse(Int64, split(temp_matrix_string[2])[1])
    row_num = parse(Int64, split(temp_matrix_string[2])[2])
    maxv = parse(Int64, temp_matrix_string[3])
    maxv_conduct = parse(Int64, conduct_matrix_string[3])

    temp_matrix = zeros(Float64, row_num * col_num)
    conduct_matrix = zeros(Float64, row_num * col_num)

    j = 1
    for (temp_string, conduct_string) in zip(temp_matrix_string[4:end], conduct_matrix_string[4:end])
        temp_line = split(temp_string)
        conduct_line = split(conduct_string)
        for (t_num, c_num) in zip(temp_line, conduct_line)
            t = parse(UInt32, t_num)
            temp_matrix[j] = io_tmin + (t * (io_tmax - io_tmin) / maxv)

            c = parse(UInt32, c_num)
            conduct_matrix[j] = c / maxv_conduct
            j += 1
        end
    end
    return (row_num, col_num, temp_matrix, conduct_matrix)
end

function make_halo_smear_cols(temp_matrix_in, conduct_matrix_in, row_num, col_num)
    # make halo rows
    temp_matrix = reshape(temp_matrix_in, col_num, row_num)' # got to transpose!!
    temp_matrix = vcat(reshape(temp_matrix[1,:], 1, col_num), temp_matrix)
    temp_matrix = vcat(temp_matrix, reshape(temp_matrix[end,:], 1, col_num))

    conduct_matrix = reshape(conduct_matrix_in, col_num, row_num)' # got to transpose!!
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

function max_not_nan(x::Float64,y::Float64)
    return ifelse(x > y, x, y)
end

function do_compute!(row_num, col_num, new_temp_matrix, temp_matrix, conduct_matrix)
    max_abs_diff = 0.0
    @inbounds @simd for j in 2:(col_num+1)
        @inbounds @simd for i in 2:(row_num+1)
            @views local_temp_matrix = temp_matrix[i-1:i+1 , j-1:j+1]

            direct_neighbors = (local_temp_matrix[1,2] + local_temp_matrix[3,2]
                             + local_temp_matrix[2,1] + local_temp_matrix[2,3])

            indirect_neighbours = (local_temp_matrix[1,1] + local_temp_matrix[1,3]
                                + local_temp_matrix[3,1] + local_temp_matrix[3,3])

            weight = conduct_matrix[i, j]
            rest_weight = 1 - weight

            new_temp_matrix[i, j] = (local_temp_matrix[2, 2] *weight
                                  + direct_neighbors * (rest_weight * direct_fixed)
                                  + indirect_neighbours * (rest_weight * indirect_fixed))
            # direct_neighbors = (temp_matrix[i - 1, j] + temp_matrix[i + 1, j]
            #                  + temp_matrix[i, j - 1] + temp_matrix[i, j + 1])
            #
            # indirect_neighbours = (temp_matrix[i - 1, j - 1] + temp_matrix[i - 1, j + 1]
            #                     + temp_matrix[i + 1, j - 1] + temp_matrix[i + 1, j + 1])
            #
            # weight = conduct_matrix[i, j]
            # rest_weight = 1 - weight
            #
            # new_temp_matrix[i, j] = (temp_matrix[i, j] * weight
            #                       + direct_neighbors * (rest_weight * direct_fixed)
            #                       + indirect_neighbours * (rest_weight * indirect_fixed))
            temp_val = abs(new_temp_matrix[i, j] - temp_matrix[i, j])
            max_abs_diff = max_not_nan(max_abs_diff, temp_val)
        end
    end
    return (max_abs_diff)
end

function report(i, new_temp_matrix, row_num, col_num, max_abs_diff, time, reporting_frequency)
    min_val = typemax(Float64)
    max_val = typemin(Float64)
    average = 0
    for k in 2:(row_num+1)
        for j in 2:(col_num+1)
            min_val = min(min_val, new_temp_matrix[k, j])
            max_val = max(max_val, new_temp_matrix[k, j])
            average += new_temp_matrix[k, j]
        end
    end

    time = time / 1e9
    flops = (row_num*col_num * (FPOPS_PER_POINT_PER_ITERATION * i + (i/reporting_frequency)))/time
    average /= (row_num * col_num)
    @printf("%-13zu % .6e % .6e % .6e % .6e % .6e % .6e \n", i, min_val, max_val, max_abs_diff, average, time, flops)
end

function main(args)
    #should be command line
    max_iter = 100
    max_abs_diff = 1000.0
    error= 0.0001
    reporting_freq = 10

    i = 0
    path_images = "images/pat1_5000x5000.pgm"
    row_num, col_num, temp_matrix1, conduct_matrix1 = read_data(path_images)
    temp_matrix, conduct_matrix = make_halo_smear_cols(temp_matrix1, conduct_matrix1, row_num, col_num)

    # make halo cells
    new_temp_matrix = zeros(Float64, row_num + 2,  col_num + 2)
    new_temp_matrix = deepcopy(temp_matrix)

    println("   Iterations        T(min)        T(max)       T(diff)        T(avg)          Time        FLOP/s")

    # time loop
    before = time_ns()
    while i < max_iter && max_abs_diff > error
        #swap new and old matrix
        new_temp_matrix, temp_matrix = temp_matrix, new_temp_matrix
        # do the actual simulation
        max_abs_diff = do_compute!(row_num, col_num, new_temp_matrix, temp_matrix, conduct_matrix)
        i += 1

        #smear columns
        new_temp_matrix[:, end] = new_temp_matrix[:, 2]
        new_temp_matrix[:, 1] = new_temp_matrix[:, end - 1]

        #check if report
        if ( i % reporting_freq == 0)
            elapsed_time = time_ns() - before
            report(i, new_temp_matrix, row_num, col_num, max_abs_diff, elapsed_time, reporting_freq)
        end
    end
    elapsed_time = time_ns() - before
    report(i, new_temp_matrix, row_num, col_num, max_abs_diff, elapsed_time, reporting_freq)
end

const direct_fixed =  0.25 * sqrt(2) / (sqrt(2)+1)
const indirect_fixed =  0.25 / (sqrt(2)+1)
const FPOPS_PER_POINT_PER_ITERATION = 12

main(ARGS)
