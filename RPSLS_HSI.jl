using Random, HDF5, Printf, Plots, Statistics, Base.Threads

function RPSLS_HSI(Lsize, reproduction_rate, selection_rate, mobility, para, rn)
    start_time = time()
    Random.seed!(rn)

    # Initialize lattice with random species (0 to 5)
    Lattice = rand(0:5, Lsize, Lsize)

    # Compute mobility and interaction rates
    M = para * 10^(-mobility * 0.1)
    eps = M * (Lsize^2) / 2
    r1 = reproduction_rate / (reproduction_rate + selection_rate + eps)
    r2 = selection_rate / (reproduction_rate + selection_rate + eps)
    r3 = eps / (reproduction_rate + selection_rate + eps)

    # Define neighbor offsets
    neighbor_shifts = [1 0; -1 0; 0 1; 0 -1]

    # Setup HDF5 file and create datasets for Lattice and NumS
    file_dir = "/Volumes/yoonD/RPS/inter/RPSLS_inter_test_$(rn).h5"
    # file_dir = "/home/ty/Desktop/yoonD/RPS/inter/RPSLS_inter_test_$(rn).h5"
    dataset1 = "/HSI/$(rn)"
    dataset2 = "/NumS/$(rn)"

    # Create the datasets with -1 for extendable dimension
    h5open(file_dir, "cw") do f
        if !haskey(f, dataset1)
            create_dataset(f, dataset1, Float64, ((1,4,5), (-1, 4, 5)); chunk=(1,4,5))
        end
        if !haskey(f, dataset2)
            create_dataset(f, dataset2, Float64, ((1,5), (-1, 5)); chunk=(1,5))
        end
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
                    if neighbor != 0 && (neighbor - main == 1 || 
                        (main == 5 && neighbor == 1) || 
                        (main == 1 && neighbor == 4) || 
                        (main == 4 && neighbor == 2) || 
                        (main == 2 && neighbor == 5) || 
                        (main == 5 && neighbor == 3) || 
                        (main == 3 && neighbor == 1))
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
        nD = count(==(4), Lattice)
        nE = count(==(5), Lattice)
        nExt = sum([nA, nB, nC, nD, nE] .== 0)

         # Compute HSI
         A_indices = findall(x -> x == 1, Lattice)
         B_indices = findall(x -> x == 2, Lattice)
         C_indices = findall(x -> x == 3, Lattice)
         D_indices = findall(x -> x == 4, Lattice)
         E_indices = findall(x -> x == 5, Lattice)

         HSI_A = compute_HSI(Lattice, A_indices, neighbor_shifts, Lsize, [2,4], [5,3])
         HSI_B = compute_HSI(Lattice, B_indices, neighbor_shifts, Lsize, [3,5], [1,4])
         HSI_C = compute_HSI(Lattice, C_indices, neighbor_shifts, Lsize, [4,1], [2,5])
         HSI_D = compute_HSI(Lattice, D_indices, neighbor_shifts, Lsize, [2,5], [3,1])
         HSI_E = compute_HSI(Lattice, E_indices, neighbor_shifts, Lsize, [1,3], [4,2])

        # Write data to HDF5 file
        h5open(file_dir, "r+") do f
            dset1 = f[dataset1]
            dset2 = f[dataset2]
            curr_size = size(dset1, 1)
            new_size = curr_size + 1

            # Resize dataset before appending new data
            HDF5.set_extent_dims(dset1, (new_size, 4, 5))
            HDF5.set_extent_dims(dset2, (new_size, 5))

            # Write the new generation's data
            dset1[new_size, :, 1] = HSI_A
            dset1[new_size, :, 2] = HSI_B
            dset1[new_size, :, 3] = HSI_C
            dset1[new_size, :, 4] = HSI_D
            dset1[new_size, :, 5] = HSI_E
            dset2[new_size, :] = [nA, nB, nC, nD, nE]
        end

        @printf("rn=%d, species=%d, %d, %d, %d, %d, nExt=%d, generation=%d\n", rn, nA, nB, nC, nD, nE, nExt, generation)

        # Stop conditions
        if mobility <= 30 && nExt == 4
            Flag = false
        end
    end
    println("Elapsed time: ", time() - start_time, " seconds")
end

# Helper function for sub2ind
function sub2ind(Lsize, row, col)
    return (col .- 1) .* Lsize .+ row
end

function compute_HSI(Lattice, indices, neighbor_shifts, Lsize, prey_species, predator_species)
    num_cells = length(indices)

    if num_cells == 0
        return [0.0, 0.0, 0.0, 0.0]
    end

    neighbors = zeros(Int, num_cells, size(neighbor_shifts, 1))

    for k in axes(neighbor_shifts, 1)
        neighbor_rows = mod1.(getindex.(indices, 1) .+ neighbor_shifts[k, 1], Lsize)
        neighbor_cols = mod1.(getindex.(indices, 2) .+ neighbor_shifts[k, 2], Lsize)
        neighbors[:, k] = Lattice[CartesianIndex.(neighbor_rows, neighbor_cols)]
    end

    prey_count = sum(in.(neighbors, Ref(Set(prey_species))), dims=2)
    predator_count = sum(in.(neighbors, Ref(Set(predator_species))), dims=2)
    empty_count = sum(neighbors .== 0, dims=2)

    return [mean((1.0*prey_count .+ 1.0*empty_count .-1.0*predator_count) ./ 4),
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
                params[2],       # reproduction_rate
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

# Main loop
Threads.@threads for rn in rn_start:rn_end
    RPSLS_HSI(Lsize, reproduction_rate, selection_rate, mobility, para, rn)
end