using Random, HDF5, Printf, Statistics, Plots
using Base.Threads

rn = 1;
Lsize = 200;
neighbor_shifts = [1 0; -1 0; 0 1; 0 -1]

# file_dir = "/home/ty/Desktop/yoonD/RPS/RPS_intra_$(rn).h5"
file_dir = "/Volumes/yoonD/RPS/RPS_intra_$(rn).h5"
data_Lattice = h5read(file_dir, "/Lattice/$(rn)")
data_Trace = h5read(file_dir, "/Trace/$(rn)")
data_NumS = h5read(file_dir, "/NumS/$(rn)")
L = size(data_NumS, 1)

# 필요한 데이터 미리 계산
A_MM_death = Vector{Vector{Float64}}(undef, L)
A_MM_birth = Vector{Vector{Float64}}(undef, L)
B_MM_death = Vector{Vector{Float64}}(undef, L)
B_MM_birth = Vector{Vector{Float64}}(undef, L)
C_MM_death = Vector{Vector{Float64}}(undef, L)
C_MM_birth = Vector{Vector{Float64}}(undef, L)

# Helper function to compute HSI
function compute_HSI!(birth_rates, death_rates, Lattice, indices, neighbor_shifts, Lsize, prey_species, predator_species)
    fill!(birth_rates, 0)
    fill!(death_rates, 0)

    @inbounds for (idx, cell) in enumerate(indices)
        row, col = Tuple(cell)
        
        for shift in neighbor_shifts
            neighbor_row = mod1(row + shift[1], Lsize)
            neighbor_col = mod1(col + shift[2], Lsize)
            neighbor_value = Lattice[neighbor_row, neighbor_col]

            if neighbor_value == predator_species
                death_rates[idx] += 1 
            elseif neighbor_value == 0
                birth_rates[idx] += 1
            end
        end
    end
end

# 사전 할당
max_cells = Lsize * Lsize
birth_rates = [zeros(Float64, max_cells) for _ in 1:nthreads()]
death_rates = [zeros(Float64, max_cells) for _ in 1:nthreads()]

@threads for i in 1:L
    tid = threadid()
    Lattice = @view data_Lattice[i, :, :]
    Trace = @view data_Trace[i, :, :]

    # 종별 인덱스 찾기
    A_indices = findall(x -> x == 1, Lattice)
    B_indices = findall(x -> x == 2, Lattice)
    C_indices = findall(x -> x == 3, Lattice)

    # HSI 계산
    compute_HSI!(birth_rates[tid], death_rates[tid], Lattice, A_indices, neighbor_shifts, Lsize, 2, 3)
    A_birth, A_death = view(birth_rates[tid], 1:length(A_indices)), view(death_rates[tid], 1:length(A_indices))
    
    compute_HSI!(birth_rates[tid], death_rates[tid], Lattice, B_indices, neighbor_shifts, Lsize, 3, 1)
    B_birth, B_death = view(birth_rates[tid], 1:length(B_indices)), view(death_rates[tid], 1:length(B_indices))
    
    compute_HSI!(birth_rates[tid], death_rates[tid], Lattice, C_indices, neighbor_shifts, Lsize, 1, 2)
    C_birth, C_death = view(birth_rates[tid], 1:length(C_indices)), view(death_rates[tid], 1:length(C_indices))

    A_T = @view Trace[A_indices]
    B_T = @view Trace[B_indices]
    C_T = @view Trace[C_indices]

    # 각 종의 최대 나이 계산
    max_A_T = isempty(A_T) ? 0 : Int(maximum(A_T))
    max_B_T = isempty(B_T) ? 0 : Int(maximum(B_T))
    max_C_T = isempty(C_T) ? 0 : Int(maximum(C_T))

    # 평균값 계산
    A_MM_death[i] = [mean(A_death[A_T .== j]) for j in 1:max_A_T if any(A_T .== j)]
    A_MM_birth[i] = [mean(A_birth[A_T .== j]) for j in 1:max_A_T if any(A_T .== j)]
    
    B_MM_death[i] = [mean(B_death[B_T .== j]) for j in 1:max_B_T if any(B_T .== j)]
    B_MM_birth[i] = [mean(B_birth[B_T .== j]) for j in 1:max_B_T if any(B_T .== j)]
    
    C_MM_death[i] = [mean(C_death[C_T .== j]) for j in 1:max_C_T if any(C_T .== j)]
    C_MM_birth[i] = [mean(C_birth[C_T .== j]) for j in 1:max_C_T if any(C_T .== j)]

    @printf("rn=%d, process = %d/%d\n", rn, i, L)
end

# 애니메이션 생성
if !isempty(A_MM_death) && !isempty(B_MM_death) && !isempty(C_MM_death)
    # y축 범위 계산
    y_max_A = maximum(maximum.(A_MM_death)) + 0.5
    y_max_B = maximum(maximum.(B_MM_death)) + 0.5
    y_max_C = maximum(maximum.(C_MM_death)) + 0.5
    y_max_NumS = maximum(maximum.(data_NumS)) + 1000

    anim = @animate for idx in 1:L
        A_MM_d = A_MM_death[idx]
        B_MM_d = B_MM_death[idx]
        C_MM_d = C_MM_death[idx]
        A_MM_b = A_MM_birth[idx]
        B_MM_b = B_MM_birth[idx]
        C_MM_b = C_MM_birth[idx]
        NumS = data_NumS[idx, :]

        p1 = plot(1:length(A_MM_d), A_MM_d, label = "A 사망률", color = :red, markershape = :circle, ylims = (0, y_max_A))
        plot!(p1, 1:length(A_MM_b), A_MM_b, label = "A 출생률", color = :pink, markershape = :square)
        
        p2 = plot(1:length(B_MM_d), B_MM_d, label = "B 사망률", color = :green, markershape = :circle, ylims = (0, y_max_B))
        plot!(p2, 1:length(B_MM_b), B_MM_b, label = "B 출생률", color = :lightgreen, markershape = :square)
        
        p3 = plot(1:length(C_MM_d), C_MM_d, label = "C 사망률", color = :blue, markershape = :circle, ylims = (0, y_max_C))
        plot!(p3, 1:length(C_MM_b), C_MM_b, label = "C 출생률", color = :lightblue, markershape = :square)
        
        p4 = plot(1:length(NumS), NumS, label = "NumS", color = :purple, markershape = :circle, ylims = (0, y_max_NumS))

        plot(p1, p2, p3, p4, layout = (2, 2), size = (800, 600), title = ["A 종" "B 종" "C 종" "NumS"])
    end

    mp4(anim, "species_animation.mp4", fps = 30)  # gif 대신 mp4로 저장
    println("애니메이션이 'species_animation.gif'로 저장되었습니다.")
else
    println("애니메이션을 생성할 데이터가 충분하지 않습니다.")
end