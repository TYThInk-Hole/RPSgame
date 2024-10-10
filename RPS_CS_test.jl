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
    # file_dir = "/Volumes/yoonD/RPS/intra/RPS_intra.h5"
    file_dir = "/home/ty/Desktop/yoonD/RPS/intra/RPS_intra.h5"
    group_name = @sprintf("intra_%.3f_%.3f_%.3f", intra1, intra2, intra3)
    dataset1 = "$group_name/Histogram/$rn"
    dataset2 = "$group_name/NumS/$rn"

    # 데이터셋 생성 (각 스레드별로 수행)
    h5open(file_dir, "r+") do f
        g = f[group_name]
        if !haskey(g["Histogram"], "$rn")
            create_dataset(g["Histogram"], "$rn", Float64, ((1, 50, 6), (-1, 50, 6)), chunk=(1, 50, 6))
        end
        if !haskey(g["NumS"], "$rn")
            create_dataset(g["NumS"], "$rn", Float64, ((1, 3), (-1, 3)), chunk=(1, 3))
        end
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
        SA = Trace .* (Lattice .== 1)
        SB = Trace .* (Lattice .== 2)
        SC = Trace .* (Lattice .== 3)
        
        spe_a = nonzeros(SA)
        spe_b = nonzeros(SB)
        spe_c = nonzeros(SC)

        bin_max = maximum([isempty(spe_a) ? 0 : maximum(spe_a),
                           isempty(spe_b) ? 0 : maximum(spe_b),
                           isempty(spe_c) ? 0 : maximum(spe_c)])

        if bin_max > 0
            bin_edges = range(0, bin_max, length=51)  # bin 중앙값

            a_values, h_a = histogram_data(spe_a, bin_edges)
            b_values, h_b = histogram_data(spe_b, bin_edges)
            c_values, h_c = histogram_data(spe_c, bin_edges)

            # Save histogram data (extend dataset1)
            h5open(file_dir, "r+") do f
                dset1 = f[dataset1]
                curr_size = size(dset1, 1)
                new_size = curr_size + 1
                HDF5.set_extent_dims(dset1, (new_size, 50, 6))

                dset1[new_size, :, 1] = a_values
                dset1[new_size, :, 2] = b_values
                dset1[new_size, :, 3] = c_values
                dset1[new_size, :, 4] = h_a
                dset1[new_size, :, 5] = h_b
                dset1[new_size, :, 6] = h_c
            end
        end

        # Count species numbers
        nA = count(==(1), Lattice)
        nB = count(==(2), Lattice)
        nC = count(==(3), Lattice)
        nExt = sum([nA, nB, nC] .== 0)

        # Save NumS data
        h5open(file_dir, "r+") do f
            dset2 = f[dataset2]
            curr_size = size(dset2, 1)
            new_size = curr_size + 1
            HDF5.set_extent_dims(dset2, (new_size, 3))
            dset2[new_size, :] = [nA, nB, nC]
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
function histogram_data(species_data, bin_edges; normalize=false)
    h = fit(Histogram, species_data, bin_edges)
    counts = h.weights
    if normalize
        counts = counts ./ sum(counts)
    end    
    x_values = range(0, stop=bin_edges[end], length=length(counts))
    return x_values, counts
end

# Helper function for sub2ind
function sub2ind(Lsize, row, col)
    return (col .- 1) .* Lsize .+ row
end

# Example usage:
function parse_command_line_args()
    if length(ARGS) != 14
        println("오류: 인수 개수가 잘못되었습니다")
        println("사용법: julia ERPS_CS_test.jl Lsize reproduction_rate selection_rate mobility intra1_start intra1_end intra1_step intra2 intra3 ext para rn_start rn_end")
        exit(1)
    end

    return (
        parse(Int, ARGS[1]),    # Lsize
        parse(Float64, ARGS[2]),  # reproduction_rate
        parse(Float64, ARGS[3]),  # selection_rate
        parse(Int, ARGS[4]),    # mobility
        parse(Float64, ARGS[5]),  # intra1_start
        parse(Float64, ARGS[6]),  # intra1_end
        parse(Float64, ARGS[7]),  # intra1_step
        parse(Float64, ARGS[8]),  # intra2
        parse(Float64, ARGS[9]),  # intra3
        parse(Int, ARGS[10]),   # ext
        parse(Float64, ARGS[11]), # para
        parse(Int, ARGS[12]),   # rn_start
        parse(Int, ARGS[13])    # rn_end
    )
end

# 메인 함수에서 파일 및 그룹 구조 생성
function main()
    Lsize, reproduction_rate, selection_rate, mobility, intra1_start, intra1_end, intra1_step, intra2, intra3, ext, para, rn_start, rn_end = parse_command_line_args()

    file_dir = "/home/ty/Desktop/yoonD/RPS/intra/RPS_intra.h5"

    for intra1 in intra1_start:intra1_step:intra1_end
        group_name = @sprintf("intra_%.3f_%.3f_%.3f", intra1, intra2, intra3)

        h5open(file_dir, "cw") do f
            if !haskey(f, group_name)
                create_group(f, group_name)
            end
            g = f[group_name]
            if !haskey(g, "Histogram")
                create_group(g, "Histogram")
            end
            if !haskey(g, "NumS")
                create_group(g, "NumS")
            end
        end

        Threads.@threads for rn in rn_start:rn_end
            RPS_intra(Lsize, reproduction_rate, selection_rate, mobility, intra1, intra2, intra3, ext, para, rn)
        end
    end
end

# 메인 함수 실행
main()