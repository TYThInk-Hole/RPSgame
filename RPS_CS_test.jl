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

    # HDF5 파일 설정
    # file_dir = "/home/ty/Desktop/yoonD/RPS/intra/RPS_intra.h5"
    file_dir = "/Volumes/yoonD/RPS/intra/RPS_intra.h5"
    group_name = @sprintf("intra_%.1f_%.1f_%.1f", intra1, intra2, intra3)
    dataset1 = "$group_name/Histogram/$rn"
    dataset2 = "$group_name/NumS/$rn"

    # 데이터를 메모리에 저장
    histogram_data = zeros(Float64, 1, 50, 6)
    nums_data = zeros(Float64, 1, 3)

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

            a_values, h_a = compute_histogram(spe_a, bin_edges)
            b_values, h_b = compute_histogram(spe_b, bin_edges)
            c_values, h_c = compute_histogram(spe_c, bin_edges)

            # Save histogram data (extend dataset1)
            histogram_data[1, :, 1] = a_values
            histogram_data[1, :, 2] = b_values
            histogram_data[1, :, 3] = c_values
            histogram_data[1, :, 4] = h_a
            histogram_data[1, :, 5] = h_b
            histogram_data[1, :, 6] = h_c
        end

        # Count species numbers
        nA = count(==(1), Lattice)
        nB = count(==(2), Lattice)
        nC = count(==(3), Lattice)
        nExt = sum([nA, nB, nC] .== 0)

        # Save NumS data
        nums_data[1, :] = [nA, nB, nC]

        @printf("rn=%d, species=%d, %d, %d, nExt=%d, generation=%d\n", rn, nA, nB, nC, nExt, generation)

        # Stop if extinction or generation limit reached
        if nExt == ext 
            Flag = false
        end
    end

    # 시뮬레이션이 끝난 후 데이터를 HDF5 파일에 저장
    lock = ReentrantLock()
    lock_and_write(file_dir, dataset1, histogram_data, dataset2, nums_data, lock)

    println("Elapsed time: ", time() - start_time, " seconds")
end

# Helper function to compute histograms for species (without normalization)
function compute_histogram(species_data, bin_edges; normalize=false)
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

function lock_and_write(file_dir, dataset1, histogram_data, dataset2, nums_data, lock)
    lock() do
        h5open(file_dir, "cw") do f
            if !haskey(f, split(dataset1, "/")[1])
                create_group(f, split(dataset1, "/")[1])
            end
            g = f[split(dataset1, "/")[1]]

            # dataset1 처리
            if haskey(g, join(split(dataset1, "/")[2:end], "/"))
                delete_object(g, join(split(dataset1, "/")[2:end], "/"))
            end
            create_dataset(g, join(split(dataset1, "/")[2:end], "/"), histogram_data)

            # dataset2 처리
            if haskey(g, join(split(dataset2, "/")[2:end], "/"))
                delete_object(g, join(split(dataset2, "/")[2:end], "/"))
            end
            create_dataset(g, join(split(dataset2, "/")[2:end], "/"), nums_data)
        end
    end
end

# Example usage:
function parse_command_line_args()
    if length(ARGS) != 14
        println("오류: 인수 개수가 올바르지 않습니다")
        println("사용법: julia ERPS_CS_test.jl Lsize reproduction_rate selection_rate mobility intra1_start intra1_end intra2_start intra2_end intra3_start intra3_end ext para start_rn end_rn")
        exit(1)
    end

    return (
        parse(Int, ARGS[1]),    # Lsize
        parse(Float64, ARGS[2]),  # reproduction_rate
        parse(Float64, ARGS[3]),  # selection_rate
        parse(Int, ARGS[4]),    # mobility
        parse(Float64, ARGS[5]),  # intra1_start
        parse(Float64, ARGS[6]),  # intra1_end
        parse(Float64, ARGS[7]),  # intra2_start
        parse(Float64, ARGS[8]),  # intra2_end
        parse(Float64, ARGS[9]),  # intra3_start
        parse(Float64, ARGS[10]), # intra3_end
        parse(Int, ARGS[11]),   # ext
        parse(Float64, ARGS[12]), # para
        parse(Int, ARGS[13]),   # start_rn
        parse(Int, ARGS[14])    # end_rn
    )
end

function run_simulation(Lsize, reproduction_rate, selection_rate, mobility, intra1_start, intra1_end, intra2_start, intra2_end, intra3_start, intra3_end, ext, para, rn)
    for intra1 in intra1_start:0.1:intra1_end
        for intra2 in intra2_start:0.1:intra2_end
            for intra3 in intra3_start:0.1:intra3_end
                RPS_intra(Lsize, reproduction_rate, selection_rate, mobility, intra1, intra2, intra3, ext, para, rn)
            end
        end
    end
end

# 메인 스크립트 부분
Lsize, reproduction_rate, selection_rate, mobility, intra1_start, intra1_end, intra2_start, intra2_end, intra3_start, intra3_end, ext, para, start_rn, end_rn = parse_command_line_args()

lock = ReentrantLock()

Threads.@threads for rn in start_rn:end_rn
    run_simulation(Lsize, reproduction_rate, selection_rate, mobility, intra1_start, intra1_end, intra2_start, intra2_end, intra3_start, intra3_end, ext, para, rn)
end

println("모든 시뮬레이션이 완료되었습니다!")