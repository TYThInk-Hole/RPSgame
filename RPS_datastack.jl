using Random, HDF5, Printf, Statistics, SparseArrays, StatsBase, Base.Threads

function RPS_intra(Lsize, reproduction_rate, selection_rate, mobility, intra1, intra2, intra3, ext, para, rn)
    start_time = time()
    Random.seed!(rn)

    # Initialize random lattice (values 0 to 3)
    Lattice = rand(0:3, Lsize, Lsize)

    # Track individual age (using sparse matrices)
    Trace = spzeros(Float64, Lsize, Lsize) .+ (Lattice .> 0)

    # Mobility and exchange rate
    M = para * 10^(-mobility * 0.1)
    eps = M * (Lsize^2) * 1/2

    # Compute rates
    rate_sum = reproduction_rate + selection_rate + intra1 + intra2 + intra3 + eps
    r1 = reproduction_rate / rate_sum
    r2 = selection_rate / rate_sum
    r3 = intra1 / rate_sum
    r4 = intra2 / rate_sum
    r5 = intra3 / rate_sum
    r6 = eps / rate_sum

    # Neighbor shifts (4 directions)
    neighbor_shifts = [1 0; -1 0; 0 1; 0 -1]

    # HDF5 file setup
    file_dir = "/Volumes/figures/RPS_intra_$rn.h5"
    # file_dir = "/home/ty/Desktop/yoonD/RPS/intra/RPS_intra_$rn.h5"
    dataset1 = "/Lattice/$rn"
    dataset2 = "/Trace/$rn"
    dataset3 = "/NumS/$rn"

    # HDF5 파일 열기 및 데이터셋 확인/생성
    h5open(file_dir, "cw") do f
        # dataset1 처리
        if haskey(f, dataset1)
            delete_object(f, dataset1)
        end
        create_dataset(f, dataset1, Float64, ((1, 200, 200), (-1, 200, 200)); chunk=(1, 200, 200))

        # dataset2 처리
        if haskey(f, dataset2)
            delete_object(f, dataset2)
        end
        create_dataset(f, dataset2, Float64, ((1, 200, 200), (-1, 200, 200)); chunk=(1, 200, 200))

        # dataset3 처리
        if haskey(f, dataset3)
            delete_object(f, dataset3)
        end
        create_dataset(f, dataset3, Float64, ((1, 3), (-1, 3)); chunk=(1, 3))
    end

    generation = 0
    Flag = true

    while Flag
        Point = true
        G = 1
        generation += 1
     
        # Generate random values and indices
        p = rand(Lsize^2)
        R = rand(1:Lsize, (Lsize^2, 2))
        rr = rand(1:4, Lsize^2)

        # Compute neighbor indices with periodic boundary conditions
        Cpre = mod1.(R .+ neighbor_shifts[rr, :], Lsize)
        C = sub2ind(Lsize, Cpre[:, 1], Cpre[:, 2])
        R_idx = sub2ind(Lsize, R[:, 1], R[:, 2])
        # Precompute conditions for all interactions
        p_r1 = p .< r1
        p_r1_r2 = p .< (r1 + r2)
        p_r1_r2_r3 = p .< (r1 + r2 + r3)
        p_r1_r2_r3_r4 = p .< (r1 + r2 + r3 + r4)
        p_r1_r2_r3_r4_r5 = p .< (r1 + r2 + r3 + r4 + r5)
        p_r1_r2_r3_r4_r5_r6 = p .< (r1 + r2 + r3 + r4 + r5 + r6)

        while Point
            for i in 1:Lsize^2
                main = Lattice[R_idx[i]]
                neighbor = Lattice[C[i]]

                main_trace = Trace[R_idx[i]]
                neighbor_trace = Trace[C[i]]

                if main == 0
                    continue
                end

                # 재생산
                if p_r1[i] 
                    if neighbor == 0
                        G += 1
                        Lattice[C[i]] = main
                        Trace[C[i]] = 1
                    end
                # 선택
                elseif p_r1_r2[i] 
                    if neighbor != 0 && (neighbor - main == 1 || (main == 3 && neighbor == 1))
                        G += 1
                        Lattice[C[i]] = 0
                        Trace[C[i]] = 0
                    end
                # 종내 선택
                elseif p_r1_r2_r3[i]
                    if main == 1 && neighbor == 1
                        G += 1
                        Lattice[C[i]] = 0
                        Trace[C[i]] = 0
                    end
                elseif p_r1_r2_r3_r4[i]
                    if main == 2 && neighbor == 2 
                        G += 1
                        Lattice[C[i]] = 0
                        Trace[C[i]] = 0
                    end
                elseif p_r1_r2_r3_r4_r5[i] 
                    if main == 3 && neighbor == 3
                        G += 1
                        Lattice[C[i]] = 0
                        Trace[C[i]] = 0
                    end
                # 이동성
                elseif p_r1_r2_r3_r4_r5_r6[i]
                    G += 1
                    Lattice[C[i]], Lattice[R_idx[i]] = main, neighbor
                    Trace[C[i]], Trace[R_idx[i]] = main_trace, neighbor_trace
                end

                if G == Lsize^2
                    Point = false
                    break
                end
            end
        end

        # Trace 업데이트 부분 수정
        Trace = Trace .+ (Lattice .!= 0)

        # Extract species ages from the trace

            # Save histogram data (extend dataset1)
            h5open(file_dir, "r+") do f
                dset1 = f[dataset1]
                curr_size = size(dset1, 1)
                new_size = curr_size + 1
                HDF5.set_extent_dims(dset1, (new_size, 200, 200))
                HDF5.set_extent_dims(dset2, (new_size, 200, 200))
                
                dset1[new_size, :, :] = Lattice
                dset2[new_size, :, :] = Trace
            end

        # Count species numbers
        nA = count(==(1), Lattice)
        nB = count(==(2), Lattice)
        nC = count(==(3), Lattice)
        nExt = sum([nA, nB, nC] .== 0)

        # Save NumS data
        h5open(file_dir, "r+") do f
            dset3 = f[dataset3]
            curr_size = size(dset3, 1)
            new_size = curr_size + 1
            HDF5.set_extent_dims(dset3, (new_size, 3))
            dset3[new_size, :] = [nA, nB, nC]
        end

        @printf("rn=%d, species=%d, %d, %d, nExt=%d, generation=%d\n", rn, nA, nB, nC, nExt, generation)

        # Stop if extinction or generation limit reached
        if nExt == ext 
            Flag = false
        end
    end

    println("Elapsed time: ", time() - start_time, " seconds")
end

# Helper function to compute histograms for species (without normalization)
# function histogram_data(species_data, bin_edges; normalize=false)
#     h = fit(Histogram, species_data, bin_edges)
#     counts = h.weights
#     if normalize
#         counts = counts ./ sum(counts)
#     end    
#     x_values = range(0, stop=bin_edges[end], length=length(counts))
#     return x_values, counts
# end

# Helper function for sub2ind
function sub2ind(Lsize, row, col)
    return (col .- 1) .* Lsize .+ row
end

# Example usage:
function get_parameters()
    println("Enter parameters in the following order:")
    println("Lsize, reproduction_rate, selection_rate, mobility, intra1, intra2, intra3, ext, para, rn_start, rn_end")
    println("Example: 200 1.0 1.0 30 2.5 1.0 0.5 1 1.0 1 1000")
    
    while true
        try
            input = readline()
            params = parse.(Float64, split(input))
            if length(params) != 11
                throw(ArgumentError("Incorrect number of parameters"))
            end
            return (
                Int(params[1]),  # Lsize
                params[2],       # reproduction_rate
                params[3],       # selection_rate
                Int(params[4]),  # mobility
                params[5],       # intra1
                params[6],       # intra2
                params[7],       # intra3
                Int(params[8]),  # ext
                params[9],       # para
                Int(params[10]), # rn_start
                Int(params[11])  # rn_end
            )
        catch e
            println("Invalid input. Please try again.")
            println("Error: ", e)
        end
    end
end

# Get parameters
Lsize, reproduction_rate, selection_rate, mobility, intra1, intra2, intra3, ext, para, rn_start, rn_end = get_parameters()

# Main loop for simulations
# for rn in rn_start:rn_end
#     RPS_intra(Lsize, reproduction_rate, selection_rate, mobility, intra1, intra2, intra3, ext, para, rn)
# end
# Main loop parallelized using Threads.@threads
Threads.@threads for rn in rn_start:rn_end
    RPS_intra(Lsize, reproduction_rate, selection_rate, mobility, intra1, intra2, intra3, ext, para, rn)
end