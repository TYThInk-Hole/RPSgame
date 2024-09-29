using Random, HDF5, Printf, Statistics, Plots
using Base.Threads
using Plots: @animate

rn = 1;
Lsize = 200;
neighbor_shifts = [1 0; -1 0; 0 1; 0 -1]  # 상, 하, 좌, 우 방향만 포함

# file_dir = "/home/ty/Desktop/yoonD/RPS/RPS_intrat_$(rn).h5"
file_dir = "/Volumes/yoonD/RPS/RPS_intrat_$(rn).h5"
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
        birth_rate = death_rate = 0

        for shift in neighbor_shifts
            neighbor_row = mod1(row + shift[1], Lsize)
            neighbor_col = mod1(col + shift[2], Lsize)
            neighbor_value = Lattice[neighbor_row, neighbor_col]

            if neighbor_value == predator_species || neighbor_value == same_species
                death_rate += 1
            elseif neighbor_value == 0
                birth_rate += 1
            end
        end

        birth_rates[row, col] = birth_rate
        death_rates[row, col] = death_rate
    end
end

# 결과를 저장할 구조체 정의
if !@isdefined(SpeciesData)
    mutable struct SpeciesData
        MM_death::Vector{Vector{Float64}}
        MM_birth::Vector{Vector{Float64}}
        death_counts::Vector{Vector{Int}}
        birth_counts::Vector{Vector{Int}}
    end
end

# 각 종별 데이터 저장
species_data = [SpeciesData(Vector{Vector{Float64}}(undef, L), 
                            Vector{Vector{Float64}}(undef, L),
                            Vector{Vector{Int}}(undef, L),
                            Vector{Vector{Int}}(undef, L)) for _ in 1:3]

# 스레드별 작업 함수
function process_chunk(chunk_start, chunk_end)
    local_birth_rates = zeros(Int, Lsize, Lsize)
    local_death_rates = zeros(Int, Lsize, Lsize)
    
    for i in chunk_start:chunk_end
        fill!(local_birth_rates, 0)
        fill!(local_death_rates, 0)

        Lattice = @view data_Lattice[i, :, :]
        Trace = @view data_Trace[i, :, :]

        for (species, prey, predator, same) in [(1, 2, 3, 1), (2, 3, 1, 2), (3, 1, 2, 3)]
            indices = findall(x -> x == species, Lattice)
            compute_HSI!(local_birth_rates, local_death_rates, Lattice, indices, neighbor_shifts, Lsize, prey, predator, same)
                
            T = Trace[indices]
            max_T = isempty(T) ? 0 : Int(maximum(T))

            death = Float64[]
            birth = Float64[]
            death_counts = zeros(Int, max_T)
            birth_counts = zeros(Int, max_T)
            for j in 1:max_T
                mask = T .== j
                if any(mask)
                    selected_indices = indices[mask]
                    selected_death_rates = local_death_rates[selected_indices]
                    selected_birth_rates = local_birth_rates[selected_indices]
                    push!(death, mean(selected_death_rates))
                    push!(birth, mean(selected_birth_rates))
                    death_counts[j] = sum(selected_death_rates)
                    birth_counts[j] = sum(selected_birth_rates)
                else
                    push!(death, NaN)
                    push!(birth, NaN)
                end
            end

            species_data[species].MM_death[i] = death
            species_data[species].MM_birth[i] = birth
            species_data[species].death_counts[i] = death_counts
            species_data[species].birth_counts[i] = birth_counts
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

# 사전 계산
y_max_A = maximum(maximum(x) for x in A_MM_death if !isempty(x); init=0) + 0.5
y_max_B = maximum(maximum(x) for x in B_MM_death if !isempty(x); init=0) + 0.5
y_max_C = maximum(maximum(x) for x in C_MM_death if !isempty(x); init=0) + 0.5
y_max_NumS = maximum(maximum.(data_NumS)) + 1000

# 사전 계산 부분 수정
function process_species_data(MM_death, MM_birth, max_age)
    data = Vector{NamedTuple{(:MM_d, :MM_b, :death_counts, :birth_counts, :death_avg, :birth_avg), 
                             Tuple{Vector{Float64}, Vector{Float64}, Vector{Int}, Vector{Int}, Vector{Float64}, Vector{Float64}}}}(undef, L)
    
    for idx in 1:L
        MM_d = replace(x -> x == 0 ? NaN : x, MM_death[idx])
        MM_b = replace(x -> x == 0 ? NaN : x, MM_birth[idx])
        
        death_counts = zeros(Int, Int(max_age))
        birth_counts = zeros(Int, Int(max_age))
        death_sum = zeros(Float64, Int(max_age))
        birth_sum = zeros(Float64, Int(max_age))
        
        for (age, d, b) in zip(1:length(MM_d), MM_d, MM_b)
            if !isnan(d) && age <= max_age
                death_counts[age] += 1
                death_sum[age] += d
            end
            if !isnan(b) && age <= max_age
                birth_counts[age] += 1
                birth_sum[age] += b
            end
        end
        
        death_avg = [count > 0 ? sum / count : NaN for (sum, count) in zip(death_sum, death_counts)]
        birth_avg = [count > 0 ? sum / count : NaN for (sum, count) in zip(birth_sum, birth_counts)]
        
        data[idx] = (MM_d=MM_d, MM_b=MM_b, death_counts=death_counts, birth_counts=birth_counts,
                     death_avg=death_avg, birth_avg=birth_avg)
    end
    return data
end

# max_age 계산 시 정수로 변환
max_age = Int(maximum(maximum.(data_NumS)))
age_classes = 1:max_age

A_data = process_species_data(A_MM_death, A_MM_birth, max_age)
B_data = process_species_data(B_MM_death, B_MM_birth, max_age)
C_data = process_species_data(C_MM_death, C_MM_birth, max_age)

# 애니메이션 생성 부분 수정
if !isempty(A_MM_death) && !isempty(B_MM_death) && !isempty(C_MM_death)
    anim = @animate for idx in 1:100:L
        A = A_data[idx]
        B = B_data[idx]
        C = C_data[idx]
        NumS = data_NumS[idx, :]

        p1 = scatter(1:length(A.MM_d), A.MM_d, label = "A death_rate", color = :red, markershape = :circle, ylims = (0, y_max_A), markersize = 3, legend = :topright)
        scatter!(p1, 1:length(A.MM_b), A.MM_b, label = "A birth_rate", color = :pink, markershape = :square, markersize = 3)
        
        p2 = scatter(1:length(species_data[1].death_counts[idx]), 
                     replace(species_data[1].death_counts[idx], 0 => NaN), 
                     label = "A death_count", color = :red, markershape = :circle, markersize = 3, legend = :topright)
        scatter!(p2, 1:length(species_data[1].birth_counts[idx]), 
                 replace(species_data[1].birth_counts[idx], 0 => NaN), 
                 label = "A birth_count", color = :blue, markershape = :square, markersize = 3)

        p3 = scatter(1:length(B.MM_d), B.MM_d, label = "B death_rate", color = :green, markershape = :circle, ylims = (0, y_max_B), markersize = 3, legend = :topright)
        scatter!(p3, 1:length(B.MM_b), B.MM_b, label = "B birth_rate", color = :lightgreen, markershape = :square, markersize = 3)
        
        p4 = scatter(1:length(species_data[2].death_counts[idx]), 
                     replace(species_data[2].death_counts[idx], 0 => NaN), 
                     label = "B death_count", color = :red, markershape = :circle, markersize = 3, legend = :topright)
        scatter!(p4, 1:length(species_data[2].birth_counts[idx]), 
                 replace(species_data[2].birth_counts[idx], 0 => NaN), 
                 label = "B birth_count", color = :blue, markershape = :square, markersize = 3)

        p5 = scatter(1:length(C.MM_d), C.MM_d, label = "C death_rate", color = :blue, markershape = :circle, ylims = (0, y_max_C), markersize = 3, legend = :topright)
        scatter!(p5, 1:length(C.MM_b), C.MM_b, label = "C birth_rate", color = :lightblue, markershape = :square, markersize = 3)
        
        p6 = scatter(1:length(species_data[3].death_counts[idx]), 
                     replace(species_data[3].death_counts[idx], 0 => NaN), 
                     label = "C death_count", color = :red, markershape = :circle, markersize = 3, legend = :topright)
        scatter!(p6, 1:length(species_data[3].birth_counts[idx]), 
                 replace(species_data[3].birth_counts[idx], 0 => NaN), 
                 label = "C birth_count", color = :blue, markershape = :square, markersize = 3)

        p7 = scatter(1:length(NumS), NumS, label = "NumS", color = :purple, markershape = :circle, ylims = (0, y_max_NumS), markersize = 3, legend = :topright)

        plot(p1, p2, p3, p4, p5, p6, p7, layout = (7, 1), size = (2000, 1200), title = ["A species" "A Age-class counts" "B species" "B Age-class counts" "C species" "C Age-class counts" "NumS"])

        @printf("애니메이션 생성 중, rn=%d, 진행 = %d/%d\n", rn, idx, L)
    end

    mp4(anim, "species_animation_intercoex.mp4", fps = 30)
    println("애니메이션이 'species_animation.mp4'로 저장되었습니다.")
else
    println("애니메이션을 생성할 데이터가 충분하지 않습니다.")
end