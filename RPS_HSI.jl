using Random, HDF5, Printf, Statistics, Base.Threads
using Plots

function RPS_HSI(Lsize, reproduction_rate, selection_rate, mobility, para, rn)
    start_time = time()
    Random.seed!(rn)

    # Initialize lattice
    Lattice = rand(0:3, Lsize, Lsize)

    # Compute mobility and rates
    M = para * 10^(-mobility * 0.1)
    eps = M * (Lsize^2) / 2
    r1 = reproduction_rate / (reproduction_rate + selection_rate + eps)
    r2 = selection_rate / (reproduction_rate + selection_rate + eps)
    r3 = eps / (reproduction_rate + selection_rate + eps)

    # Define neighbor offsets
    neighbor_shifts = [1 0; -1 0; 0 1; 0 -1]

    # Setup HDF5 file
    # file_dir = "/Volumes/yoonD/RPS/RPS_inter_test_$(rn).h5"
    file_dir = "/home/ty/Desktop/yoonD/RPS/inter/RPS_inter_test_$(rn).h5"
    dataset1 = "/HSI/$(rn)"
    dataset2 = "/NumS/$(rn)"
    
    h5open(file_dir, "cw") do f
        # 데이터셋이 존재하면 삭제
        haskey(f, dataset1) && delete_object(f, dataset1)
        haskey(f, dataset2) && delete_object(f, dataset2)
        
        # 새 데이터셋 생성
        create_dataset(f, dataset1, Float64, ((1,4,3), (-1, 4, 3)); chunk=(1,4,3))
        create_dataset(f, dataset2, Float64, ((1,3), (-1, 3)); chunk=(1,3))
    end

    generation = 0
    Flag = true

    while Flag
        Point = true
        G = 1
        generation += 1

        # Pre-generate random numbers and coordinates for all cells
        p = rand(Lsize^2)
        R = rand(1:Lsize, (Lsize^2, 2))
        rr = rand(1:4, Lsize^2)

        # Compute the neighbor indices with boundary condition handling
        Cpre = mod1.(R .+ neighbor_shifts[rr, :], Lsize)
        C = sub2ind(Lsize, Cpre[:, 1], Cpre[:, 2])
        R_idx = sub2ind(Lsize, R[:, 1], R[:, 2])

        # Precompute common conditions outside the loop
        p_r1 = p .< r1
        p_r1_r2 = p .< (r1 + r2)
        p_r1_r2_r3 = p .< (r1 + r2 + r3)

        while Point
            for i in 1:Lsize^2
                main = Lattice[R_idx[i]]
                neighbor = Lattice[C[i]]

                if main == 0
                    continue
                end

                # Reproduction
                if p_r1[i]
                    if neighbor == 0
                        G += 1
                        Lattice[C[i]] = main
                    end

                # Selection
                elseif p_r1_r2[i]
                    if neighbor != 0 && (neighbor - main == 1 || (main == 3 && neighbor == 1))
                        G += 1
                        Lattice[C[i]] = 0
                    end

                # Mobility (Swap)
                elseif p_r1_r2_r3[i]
                    G += 1
                    Lattice[C[i]] = main
                    Lattice[R_idx[i]] = neighbor
                end

                if G == Lsize^2
                    Point = false
                    break
                end
            end
        end

        # Count species
        nA = count(==(1), Lattice)
        nB = count(==(2), Lattice)
        nC = count(==(3), Lattice)
        nExt = sum([nA, nB, nC] .== 0)

        # Compute HSI
        A_indices = findall(x -> x == 1, Lattice)
        B_indices = findall(x -> x == 2, Lattice)
        C_indices = findall(x -> x == 3, Lattice)

        HSI_A = compute_HSI(Lattice, A_indices, A_indices, neighbor_shifts, Lsize, 2, 3, 1)
        HSI_B = compute_HSI(Lattice, B_indices, B_indices, neighbor_shifts, Lsize, 3, 1, 2)
        HSI_C = compute_HSI(Lattice, C_indices, C_indices, neighbor_shifts, Lsize, 1, 2, 3)

        # Write data to HDF5 file
        h5open(file_dir, "r+") do f
            dset1 = f[dataset1]
            dset2 = f[dataset2]
            curr_size = size(dset1, 1)
            new_size = curr_size + 1
            
            # 데이터셋 크기 확장
            HDF5.set_extent_dims(dset1, (new_size, 4, 3))
            HDF5.set_extent_dims(dset2, (new_size, 3))
            
            # 데이터 쓰기
            dset1[new_size, :, 1] = HSI_A
            dset1[new_size, :, 2] = HSI_B
            dset1[new_size, :, 3] = HSI_C  
            dset2[new_size, :] = [nA, nB, nC]
        end

        @printf("rn=%d, species=%d, %d, %d, nExt=%d, generation=%d\n", rn, nA, nB, nC, nExt, generation)

        # Stop conditions
        if nExt == 2
            Flag = false
        end
    end
    
    println("Elapsed time: ", time() - start_time, " seconds")
end

# Helper function for sub2ind
function sub2ind(Lsize, row, col)
    return (col .- 1) .* Lsize .+ row
end

# Helper function to compute HSI
function compute_HSI(Lattice, rows, cols, neighbor_shifts, Lsize, prey_species, predator_species, same_species)
    num_cells = length(rows)

    if num_cells == 0
        return [0.0, 0.0, 0.0, 0.0]  # 4개의 float 값을 반환
    end

    neighbors = zeros(Int, num_cells, size(neighbor_shifts, 1))

    for k in axes(neighbor_shifts, 1)
        neighbor_rows = mod1.(getindex.(rows, 1) .+ neighbor_shifts[k, 1], Lsize)
        neighbor_cols = mod1.(getindex.(cols, 2) .+ neighbor_shifts[k, 2], Lsize)
        neighbors[:, k] = Lattice[CartesianIndex.(neighbor_rows, neighbor_cols)]
    end

    prey_count = sum(neighbors .== prey_species, dims=2)
    predator_count = sum(neighbors .== predator_species, dims=2)
    empty_count = sum(neighbors .== 0, dims=2)
    same_count = sum(neighbors .== same_species, dims=2)

    return [
        mean((1.0*prey_count .+ 1.0*empty_count .-1.0*predator_count) ./ 4),
        mean((0.0*prey_count .+ 1.0*empty_count .-1.0*predator_count) ./ 4),
        mean((0.0*prey_count .+ 1.0*empty_count .-2.0*predator_count) ./ 4),
        mean((0.0*prey_count .+ 2.0*empty_count .-1.0*predator_count) ./ 4)
    ]
end

# Example usage:
function get_parameters()
    println("Enter parameters in the following order:")
    println("Lsize, reproduction_rate, selection_rate, mobility, para, rn_start, rn_end")
    println("Example: 200 1.0 1.0 30 3.0 1 1000")
    
    while true
        try
            input = readline()
            params = parse.(Float64, split(input))
            if length(params) != 7
                throw(ArgumentError("Incorrect number of parameters"))
            end
            return (
                Int(params[1]),  # Lsize
                params[2],       # reproduction_rate2222
                params[3],       # selection_rate
                Int(params[4]),  # mobility
                params[5],       # para
                Int(params[6]),  # rn_start
                Int(params[7])   # rn_end
            )
        catch e
            println("Invalid input. Please try again.")
            println("Error: ", e)
        end
    end
end

# Get parameters
Lsize, reproduction_rate, selection_rate, mobility, para, rn_start, rn_end = get_parameters()

#Main loop
Threads.@threads for rn in rn_start:rn_end
    RPS_HSI(Lsize, reproduction_rate, selection_rate, mobility, para, rn)
end

# rn=1
# h5open("/Volumes/yoonJL/RPS_inter_test.h5", "r") do file
#     hsi_data = read(file["/HSI/$(rn)"])
#     nums_data = read(file["/NumS/$(rn)"])
    
#     println("HSI 데이터:")
#     println(hsi_data)
#     println("\n종 수 데이터:")
#     println(nums_data)
# end


# end

