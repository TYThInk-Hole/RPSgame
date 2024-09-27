using Random, HDF5, Printf, Statistics, Plots

rn = 1;
Lsize = 200;
neighbor_shifts = [1 0; -1 0; 0 1; 0 -1]

file_dir = "/home/ty/Desktop/yoonD/RPS/RPS_intra_$(rn).h5"

data_Lattice = h5read(file_dir, "/Lattice/$(rn)")
data_Trace = h5read(file_dir, "/Trace/$(rn)")
data_NumS = h5read(file_dir, "/NumS/$(rn)")
L = size(data_NumS, 1)

# 필요한 데이터 미리 계산
A_MM_list = []
B_MM_list = []
C_MM_list = []

# Helper function to compute HSI
function compute_HSI(Lattice, indices, neighbor_shifts, Lsize, prey_species, predator_species)
    num_cells = length(indices)

    if num_cells == 0
        return 0.0, zeros(Float64, size(Lattice))
    end

    result_matrix = zeros(Float64, size(Lattice))

    for idx in indices
        row, col = Tuple(idx)
        total_count = 0.0

        for shift in eachrow(neighbor_shifts)
            neighbor_row = mod1(row + shift[1], Lsize)
            neighbor_col = mod1(col + shift[2], Lsize)
            neighbor_value = Lattice[neighbor_row, neighbor_col]

            if neighbor_value == predator_species
                total_count += 10
            elseif neighbor_value == 0
                total_count += 1
            end
        end

        result_matrix[row, col] = total_count
    end

    hsi_value = result_matrix[indices]
    return hsi_value, result_matrix
end

for i in 1:L
    Lattice = data_Lattice[i, :, :]
    Trace = data_Trace[i, :, :]

    # 종별 인덱스 찾기
    A_indices = findall(Lattice .== 1)
    B_indices = findall(Lattice .== 2)
    C_indices = findall(Lattice .== 3)

    # HSI 계산
    _, count_matrix_A = compute_HSI(Lattice, A_indices, neighbor_shifts, Lsize, 2, 3)
    _, count_matrix_B = compute_HSI(Lattice, B_indices, neighbor_shifts, Lsize, 3, 1)
    _, count_matrix_C = compute_HSI(Lattice, C_indices, neighbor_shifts, Lsize, 1, 2)

    A_T = (Lattice .== 1) .* Trace
    B_T = (Lattice .== 2) .* Trace
    C_T = (Lattice .== 3) .* Trace

    # 각 종의 최대 나이 계산
    max_A_T = Int(maximum(A_T))
    max_B_T = Int(maximum(B_T))
    max_C_T = Int(maximum(C_T))

    # 평균값 저장을 위한 배열 초기화
    A_MM = Float64[]
    B_MM = Float64[]
    C_MM = Float64[]

    # 평균값 계산
    for j in 1:max_A_T
        indices = findall(A_T .== j)
        if !isempty(indices)
            values = count_matrix_A[indices]
            push!(A_MM, mean(values))
        end
    end

    for j in 1:max_B_T
        indices = findall(B_T .== j)
        if !isempty(indices)
            values = count_matrix_B[indices]
            push!(B_MM, mean(values))
        end
    end

    for j in 1:max_C_T
        indices = findall(C_T .== j)
        if !isempty(indices)
            values = count_matrix_C[indices]
            push!(C_MM, mean(values))
        end
    end

    push!(A_MM_list, A_MM)
    push!(B_MM_list, B_MM)
    push!(C_MM_list, C_MM)

    @printf("rn=%d, process = %d/%d\n", rn, i, L)
end

# 애니메이션 생성
if !isempty(A_MM_list) && !isempty(B_MM_list) && !isempty(C_MM_list)
    anim = @animate for idx in 1:L
        A_MM = A_MM_list[idx]
        B_MM = B_MM_list[idx]
        C_MM = C_MM_list[idx]
        NumS = data_NumS[idx, :]

        p1 = scatter(1:length(A_MM), A_MM, label = "species A", color = :red, markershape = :circle)
        p2 = scatter(1:length(B_MM), B_MM, label = "species B", color = :green, markershape = :circle)
        p3 = scatter(1:length(C_MM), C_MM, label = "species C", color = :blue, markershape = :circle)
        p4 = scatter(1:length(NumS), NumS, label = "NumS", color = :purple, markershape = :circle)

        plot(p1, p2, p3, p4, layout = (2, 2), size = (800, 600))
    end

    mp4(anim, "species_animation.mp4", fps = 30)  # gif 대신 mp4로 저장
    println("애니메이션이 'species_animation.gif'로 저장되었습니다.")
else
    println("애니메이션을 생성할 데이터가 충분하지 않습니다.")
end