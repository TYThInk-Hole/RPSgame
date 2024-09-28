using Random, HDF5, Printf, Statistics, Plots
using Base.Threads
using Plots: @animate

rn = 1;
Lsize = 200;
neighbor_shifts = [1 0; -1 0; 0 1; 0 -1]

file_dir = "/home/ty/Desktop/yoonD/RPS/RPS_intrat_$(rn).h5"
# file_dir = "/Volumes/yoonD/RPS/RPS_intra_$(rn).h5"
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
function compute_HSI!(birth_rates, death_rates, Lattice, indices, neighbor_shifts, Lsize, prey_species, predator_species, same_species)
    @inbounds for cell in indices
        row, col = cell.I
        birth_count = death_count = 0
        
        for shift in neighbor_shifts
            neighbor_row = mod1(row + shift[1], Lsize)
            neighbor_col = mod1(col + shift[2], Lsize)
            neighbor_value = Lattice[neighbor_row, neighbor_col]

            death_count += (neighbor_value == predator_species)
            death_count += (neighbor_value == same_species)
            birth_count += (neighbor_value == 0)
        end
        
        birth_rates[row, col] = birth_count
        death_rates[row, col] = death_count
    end
end

# 결과를 저장할 구조체 정의
struct SpeciesData
    MM_death::Vector{Vector{Float64}}
    MM_birth::Vector{Vector{Float64}}
end

# 각 종별 데이터 저장
species_data = [SpeciesData(Vector{Vector{Float64}}(undef, L), Vector{Vector{Float64}}(undef, L)) for _ in 1:3]

# 스레드별 작업 함수
function process_chunk(chunk_start, chunk_end)
    local_birth_rates = zeros(Int, Lsize, Lsize)
    local_death_rates = zeros(Int, Lsize, Lsize)

    for i in chunk_start:chunk_end
        Lattice = @view data_Lattice[i, :, :]
        Trace = @view data_Trace[i, :, :]

        for (species, prey, predator, same) in [(1, 2, 3, 1), (2, 3, 1, 2), (3, 1, 2, 3)]
            indices = findall(x -> x == species, Lattice)
            compute_HSI!(local_birth_rates, local_death_rates, Lattice, indices, neighbor_shifts, Lsize, prey, predator, same)
            
            T = Trace[indices]
            max_T = isempty(T) ? 0 : Int(maximum(T))

            death = Float64[]
            birth = Float64[]
            for j in 1:max_T
                mask = T .== j
                if any(mask)
                    push!(death, mean(local_death_rates[indices[mask]]))
                    push!(birth, mean(local_birth_rates[indices[mask]]))
                end
            end

            species_data[species].MM_death[i] = death
            species_data[species].MM_birth[i] = birth
        end

        @printf("rn=%d, process = %d/%d\n", rn, i, L)
    end
end

# 멀티스레딩 실행
@threads for t in 1:nthreads()
    chunk_size = cld(L, nthreads())
    chunk_start = (t - 1) * chunk_size + 1
    chunk_end = min(t * chunk_size, L)
    process_chunk(chunk_start, chunk_end)
end

# 결과 추출
A_MM_death, A_MM_birth = species_data[1].MM_death, species_data[1].MM_birth
B_MM_death, B_MM_birth = species_data[2].MM_death, species_data[2].MM_birth
C_MM_death, C_MM_birth = species_data[3].MM_death, species_data[3].MM_birth

# 애니메이션 생성
if !isempty(A_MM_death) && !isempty(B_MM_death) && !isempty(C_MM_death)
    # y축 범위 계산
    y_max_A = maximum(maximum(x) for x in A_MM_death if !isempty(x); init=0) + 0.5
    y_max_B = maximum(maximum(x) for x in B_MM_death if !isempty(x); init=0) + 0.5
    y_max_C = maximum(maximum(x) for x in C_MM_death if !isempty(x); init=0) + 0.5
    y_max_NumS = maximum(maximum.(data_NumS)) + 1000

    anim = @animate for idx in 1:L
        A_MM_d = A_MM_death[idx]
        B_MM_d = B_MM_death[idx]
        C_MM_d = C_MM_death[idx]
        A_MM_b = A_MM_birth[idx]
        B_MM_b = B_MM_birth[idx]
        C_MM_b = C_MM_birth[idx]
        NumS = data_NumS[idx, :]

        p1 = scatter(1:length(A_MM_d), A_MM_d, label = "A death_rate", color = :darkred, markershape = :circle, ylims = (0, y_max_A), legend = :topright)
        scatter!(p1, 1:length(A_MM_b), A_MM_b, label = "A birth_rate", color = :pink, markershape = :square)
        
        p2 = scatter(1:length(B_MM_d), B_MM_d, label = "B death_rate", color = :darkgreen, markershape = :circle, ylims = (0, y_max_B), legend = :topright)
        scatter!(p2, 1:length(B_MM_b), B_MM_b, label = "B birth_rate", color = :lightgreen, markershape = :square)
        
        p3 = scatter(1:length(C_MM_d), C_MM_d, label = "C death_rate", color = :darkblue, markershape = :circle, ylims = (0, y_max_C), legend = :topright)
        scatter!(p3, 1:length(C_MM_b), C_MM_b, label = "C birth_rate", color = :lightblue, markershape = :square)
        
        p4 = scatter(1:length(NumS), NumS, label = "NumS", color = :purple, markershape = :circle, ylims = (0, y_max_NumS), legend = :topright)

        plot(p1, p2, p3, p4, layout = (4, 1), size = (1200, 800), title = ["A species" "B species" "C species" "NumS"])
        @printf("make animation, rn=%d, process = %d/%d\n", rn, idx, L)
    end

    mp4(anim, "1_species_animation.mp4", fps = 30)
    println("애니메이션이 '1_species_animation.mp4'로 저장되었습니다.")
else
    println("애니메이션을 생성할 데이터가 충분하지 않습니다.")
end

