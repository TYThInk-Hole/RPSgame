using Random, HDF5, Printf

function RPS_inter_individualtt(Lsize, reproduction_rate, selection_rate, mobility, para, rn)
    start_time = time()
    Random.seed!(rn)

    # Initialize lattice
    Lattice = rand(0:3, Lsize, Lsize)

    # Compute mobility and rates
    M = para * 10^(-mobility * 0.1)
    eps = M * (Lsize^2) / 2
    r1 = reproduction_rate / (reproduction_rate + selection_rate + eps)
    r2 = selection_rate / (reproduction_rate + eps)
    r3 = eps / (reproduction_rate + selection_rate + eps)

    # Define neighbor offset
    A = [1 0; -1 0; 0 1; 0 -1]

    # Setup HDF5 file
    file_dir = "/Volumes/yoonOWC/RPS_inter_test.h5"
    dataset1 = "/$(rn)/Lattice"
    dataset2 = "/$(rn)/NumS"

    # Create HDF5 datasets
    h5open(file_dir, "cw") do f
        create_dataset(f, dataset1, datatype(Float64), dataspace(Lsize, Lsize, 1000000); chunk=(Lsize, Lsize, 1))
        create_dataset(f, dataset2, datatype(Float64), dataspace(1000000, 3); chunk=(1, 3))
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
        Cpre = mod1.(R .+ A[rr, :] .- 1, Lsize)
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

        # Write data to HDF5 file
        h5open(file_dir, "r+") do f
            dset1 = f[dataset1]
            dset2 = f[dataset2]
            dset1[:, :, generation] = Lattice
            dset2[generation, :] = [nA, nB, nC]
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

# Example usage:
Lsize = 200
reproduction_rate = 1.0
selection_rate = 1.0
mobility = 30
para = 3.0
for rn in 1:100
    RPS_inter_individualtt(Lsize, reproduction_rate, selection_rate, mobility, para, rn)
end