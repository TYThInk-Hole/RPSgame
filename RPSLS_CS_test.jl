using Random, HDF5, Printf, Statistics, SparseArrays, StatsBase, Base.Threads

function RPSLS_intra(Lsize, reproduction_rate, selection_rate, mobility, intra1, intra2, intra3, intra4, intra5, ext, para, rn)
    start_time = time()
    Random.seed!(rn)

    # Initialize random lattice (values 0 to 3)
    Lattice = rand(0:5, Lsize, Lsize)

    # Track individual age (using sparse matrices)
    Trace = spzeros(Float64, Lsize, Lsize) .+ (Lattice .> 0)

    # Mobility and exchange rate
    M = para * 10^(-mobility * 0.1)
    eps = M * (Lsize^2) * 1/2

    # Compute rates
    rate_sum = reproduction_rate + selection_rate + intra1 + intra2 + intra3 + intra4 + intra5 + eps
    r1 = reproduction_rate / rate_sum
    r2 = selection_rate / rate_sum
    r3 = intra1 / rate_sum
    r4 = intra2 / rate_sum
    r5 = intra3 / rate_sum
    r6 = intra4 / rate_sum
    r7 = intra5 / rate_sum
    r8 = eps / rate_sum

    # Neighbor shifts (4 directions)
    neighbor_shifts = [1 0; -1 0; 0 1; 0 -1]

    # HDF5 file setup
    file_dir = "/home/ty/Desktop/yoonD/RPS/intra/RPSLS_intra.h5"
    group_name = @sprintf("intra_%.3f_%.3f_%.3f_%.3f_%.3f", intra1, intra2, intra3, intra4, intra5)
    dataset1 = "$group_name/Histogram/$rn"
    dataset2 = "$group_name/NumS/$rn"

    # HDF5 파일 열기 및 데이터셋 확인/생성
    h5open(file_dir, "cw") do f
        # dataset1 처리
        if haskey(f, dataset1)
            delete_object(f, dataset1)
        end
        create_dataset(f, dataset1, Float64, ((1, 50, 6), (-1, 50, 6)); chunk=(1, 50, 6))

        # dataset2 처리
        if haskey(f, dataset2)
            delete_object(f, dataset2)
        end
        create_dataset(f, dataset2, Float64, ((1, 5), (-1, 5)); chunk=(1, 5))

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
        p_r1_r2_r3_r4_r5_r6_r7 = p .< (r1 + r2 + r3 + r4 + r5 + r6 + r7)
        p_r1_r2_r3_r4_r5_r6_r7_r8 = p .< (r1 + r2 + r3 + r4 + r5 + r6 + r7 + r8)
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
                elseif p_r1_r2_r3_r4_r5_r6[i] 
                    if main == 4 && neighbor == 4
                        G += 1
                        Lattice[C[i]] = 0
                        Trace[C[i]] = 0
                    end
                elseif p_r1_r2_r3_r4_r5_r6_r7[i] 
                    if main == 5 && neighbor == 5
                        G += 1
                        Lattice[C[i]] = 0
                        Trace[C[i]] = 0
                    end
                # 이동성
                elseif p_r1_r2_r3_r4_r5_r6_r7_r8[i]
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
        SD = Trace .* (Lattice .== 4)
        SE = Trace .* (Lattice .== 5)
        
        spe_a = nonzeros(SA)
        spe_b = nonzeros(SB)
        spe_c = nonzeros(SC)
        spe_d = nonzeros(SD)
        spe_e = nonzeros(SE)

        bin_max = maximum([isempty(spe_a) ? 0 : maximum(spe_a),
                           isempty(spe_b) ? 0 : maximum(spe_b),
                           isempty(spe_c) ? 0 : maximum(spe_c),
                           isempty(spe_d) ? 0 : maximum(spe_d),
                           isempty(spe_e) ? 0 : maximum(spe_e)])

        if bin_max > 0            
            bin_edges = range(0, bin_max, length=51)  # bin 중앙값

            x_values, h_a = histogram_data(spe_a, bin_edges)
            _, h_b = histogram_data(spe_b, bin_edges)
            _, h_c = histogram_data(spe_c, bin_edges)
            _, h_d = histogram_data(spe_d, bin_edges)
            _, h_e = histogram_data(spe_e, bin_edges)

            # Save histogram data (extend dataset1)
            h5open(file_dir, "r+") do f
                dset1 = f[dataset1]
                curr_size = size(dset1, 1)
                new_size = curr_size + 1
                HDF5.set_extent_dims(dset1, (new_size, 50, 6))

                dset1[new_size, :, 1] = x_values
                dset1[new_size, :, 2] = h_a
                dset1[new_size, :, 3] = h_b
                dset1[new_size, :, 4] = h_c
                dset1[new_size, :, 5] = h_d
                dset1[new_size, :, 6] = h_e
            end
        end

        # Count species numbers
        nA = count(==(1), Lattice)
        nB = count(==(2), Lattice)
        nC = count(==(3), Lattice)
        nD = count(==(4), Lattice)
        nE = count(==(5), Lattice)
        nExt = sum([nA, nB, nC, nD, nE] .== 0)

        # Save NumS data
        h5open(file_dir, "r+") do f
            dset2 = f[dataset2]
            curr_size = size(dset2, 1)
            new_size = curr_size + 1
            HDF5.set_extent_dims(dset2, (new_size, 5))
            dset2[new_size, :] = [nA, nB, nC, nD, nE]
        end

        @printf("rn=%d, species=%d, %d, %d, %d, %d, nExt=%d, generation=%d\n", rn, nA, nB, nC, nD, nE, nExt, generation)

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
    x_values = (bin_edges[1:end-1] .+ bin_edges[2:end]) ./ 2  # bin 중앙값
    return x_values, counts
end

# Helper function for sub2ind
function sub2ind(Lsize, row, col)
    return (col .- 1) .* Lsize .+ row
end

# 명령줄 인수를 파싱하는 함수로 get_parameters() 함수를 대체
function parse_command_line_args()
    if length(ARGS) != 13
        println("오류: 인수 개수가 잘못되었습니다")
        println("사용법: julia RPSLS_CS_test.jl Lsize reproduction_rate selection_rate mobility intra1 intra2 intra3 intra4 intra5 ext para rn_start rn_end")
        exit(1)
    end

    return (
        parse(Int, ARGS[1]),    # Lsize
        parse(Float64, ARGS[2]),  # reproduction_rate
        parse(Float64, ARGS[3]),  # selection_rate
        parse(Int, ARGS[4]),    # mobility
        parse(Float64, ARGS[5]),  # intra1
        parse(Float64, ARGS[6]),  # intra2
        parse(Float64, ARGS[7]),  # intra3
        parse(Float64, ARGS[8]),  # intra4
        parse(Float64, ARGS[9]),  # intra5
        parse(Int, ARGS[10]),   # ext
        parse(Float64, ARGS[11]), # para
        parse(Int, ARGS[12]),   # rn_start
        parse(Int, ARGS[13])    # rn_end
    )
end

# 메인 함수에서 파일 및 그룹 구조 생성
function main()
    Lsize, reproduction_rate, selection_rate, mobility, intra1, intra2, intra3, intra4, intra5, ext, para, rn_start, rn_end = parse_command_line_args()

    file_dir = "/home/ty/Desktop/yoonD/RPS/intra/RPSLS_intra.h5"
    group_name = @sprintf("intra_%.3f_%.3f_%.3f_%.3f_%.3f", intra1, intra2, intra3, intra4, intra5)

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

    # 스레드 실행
    Threads.@threads for rn in rn_start:rn_end
        RPSLS_intra(Lsize, reproduction_rate, selection_rate, mobility, intra1, intra2, intra3, intra4, intra5, ext, para, rn)
    end
end

# 메인 함수 실행
main()