using HDF5
using LinearAlgebra
using Statistics
using Printf
using Base.Threads

# Function to compute cosine similarity between two vectors or matrices
function cosine_similarity(a, b)
    numerator = dot(a, b)
    denominator = norm(a) * norm(b)
    return denominator == 0 ? 0.0 : numerator / denominator
end

# Main function
function run_simulation()
    # Parameters
    inte = 1:100
    j_vals = 0:0.01:1
    intra = 0.1:0.1:1.5
    intra2, intra3 = 1.0, 0.5
    st = 200
    mmat = zeros(100, 101)

    for i in 1:length(intra)
        intra1 = intra[i]
        S_probability = zeros(100, 101, 2)

        for ii in [10]  # MATLAB code only runs for ii = 10
            mmat[ii, 1] = inte[ii]
            intraval = inte[ii]

            Threads.@threads for rn in 1:100  # 병렬 처리 추가
                # Construct dataset paths
                intra_path = @sprintf("/intra_%.3f_%.3f_%.3f", intra1, intra2, intra3)
                histo_path = intra_path * "/Histogram/$(rn)"
                nums_path = intra_path * "/NumS/$(rn)"

                # Read data from HDF5 files
                h5_filename = "/Volumes/yoondata/RPS/intra/RPS_intra.h5"
                data1 = h5read(h5_filename, histo_path)[2:end, :, :]
                NumS = h5read(h5_filename, nums_path)[2:end, :]

                L = size(NumS, 1)
                cos_vals = zeros(L - intraval, 3)

                for idx in 1:(L - intraval)
                    # Extract and concatenate data for cosine similarity
                    vec_a1 = vcat(data1[idx, :, 1], data1[idx, :, 4])
                    vec_a2 = vcat(data1[idx + intraval, :, 1], data1[idx + intraval, :, 4])
                    a_similarity = cosine_similarity(vec_a1, vec_a2)

                    vec_b1 = vcat(data1[idx, :, 2], data1[idx, :, 5])
                    vec_b2 = vcat(data1[idx + intraval, :, 2], data1[idx + intraval, :, 5])
                    b_similarity = cosine_similarity(vec_b1, vec_b2)

                    vec_c1 = vcat(data1[idx, :, 3], data1[idx, :, 6])
                    vec_c2 = vcat(data1[idx + intraval, :, 3], data1[idx + intraval, :, 6])
                    c_similarity = cosine_similarity(vec_c1, vec_c2)

                    # Store the similarities
                    cos_vals[idx, :] = [a_similarity, b_similarity, c_similarity]
                end

                # Calculate distance matrices
                dist_a = hcat(norm.(cos_vals[st:end, 1] .- cos_vals[st:end, 2], 2),
                              norm.(cos_vals[st:end, 1] .- cos_vals[st:end, 3], 2))
                dist_b = hcat(norm.(cos_vals[st:end, 1] .- cos_vals[st:end, 2], 2),
                              norm.(cos_vals[st:end, 2] .- cos_vals[st:end, 3], 2))
                dist_c = hcat(norm.(cos_vals[st:end, 1] .- cos_vals[st:end, 3], 2),
                              norm.(cos_vals[st:end, 2] .- cos_vals[st:end, 3], 2))

                # 평균 거리가 0.10을 초과하는 최소 인덱스 찾기
                md_candidates = [
                    findfirst(x -> x > 0.10, vec(mean(dist_a, dims=2))),
                    findfirst(x -> x > 0.10, vec(mean(dist_b, dims=2))),
                    findfirst(x -> x > 0.10, vec(mean(dist_c, dims=2)))
                ]
                # md_candidates에서 nothing이 아닌 요소만 필터링
                valid_md_candidates = filter(!isnothing, md_candidates)
                md = isempty(valid_md_candidates) ? st : minimum(valid_md_candidates)

                th = zeros(length(j_vals))

                for (k, dth) in enumerate(j_vals)
                    success_e, success_s = 0, 0

                    for species in 1:3
                        cos_condition = any(cos_vals[md:end, species] .< dth)
                        nums_condition = any(NumS[2:end, species] .== 0)

                        if cos_condition && nums_condition
                            success_e += 1
                        elseif !cos_condition && !nums_condition
                            success_s += 1
                        end
                    end

                    S_probability[rn, k, :] = [success_e, success_s]

                    if success_e + success_s == 3
                        th[k] = dth
                    end
                end

                max_th = maximum(th)
                println(@sprintf("rn = %d, intraval = %d, max_th = %g", rn, intraval, max_th))
                mmat[ii, rn + 1] = max_th
            end
        end

        # Save S_probability to HDF5 file
        h5_filename = "/Volumes/yoondata/RPS/intra/RPS_intra.h5"
        dataset_path = @sprintf("/intra_%.3f_%.3f_%.3f/P", intra1, intra2, intra3)
        h5open(h5_filename, "cw") do file
            write(file, dataset_path, S_probability)
        end
    end
end

# Run the simulation
run_simulation()
