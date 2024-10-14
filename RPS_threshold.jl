using HDF5
using LinearAlgebra
using Statistics
using Plots
using Base.Threads
using Printf

# Function to compute cosine similarity
function cosine_similarity(A, B)
    dot_product = dot(A, B)
    norm_A = norm(A)
    norm_B = norm(B)
    if norm_A == 0 || norm_B == 0
        return 0.0
    else
        return dot_product / (norm_A * norm_B)
    end
end

# First code block
function compute_mmat()
    inte = 20:20
    j_values = 0:0.01:1
    intra1_values = 0.1:0.1:1.5
    intra2 = 1.0
    intra3 = 0.5
    st = 200
    mmat = zeros(Float64, 100, 100 + 1)

    for (ii, intraval) in enumerate(inte)
        mmat[ii, 1] = intraval
        @threads for rn in 1:100
            for intra1 in intra1_values
                intra1_str = @sprintf("%.3f", intra1)
                intra2_str = @sprintf("%.3f", intra2)
                intra3_str = @sprintf("%.3f", intra3)
                # File paths and dataset names
                h5_filename = "/home/ty/Desktop/yoonD/RPS/intra/RPS_intra.h5"
                histogram_path = "/intra_$(intra1_str)_$(intra2_str)_$(intra3_str)/Histogram/$(rn)"
                nums_path = "/intra_$(intra1_str)_$(intra2_str)_$(intra3_str)/NumS/$(rn)"

                # Open HDF5 file once per thread
                h5file = h5open(h5_filename, "r")
                try
                    if haskey(h5file, histogram_path)
                        data = read(h5file[histogram_path])
                        data1 = data[2:end, :, :]
                        data2 = read(h5file[nums_path])
                        NumS = data2[2:end, :]
                        L = size(NumS, 1)
                        cos_sim = zeros(L - intraval, 3)

                        # Compute cosine similarities
                        for i in 1:(L - intraval)
                            a = cosine_similarity(vcat(data1[i, :, 1], data1[i, :, 4]),
                                                  vcat(data1[i + intraval, :, 1], data1[i + intraval, :, 4]))
                            b = cosine_similarity(vcat(data1[i, :, 2], data1[i, :, 5]),
                                                  vcat(data1[i + intraval, :, 2], data1[i + intraval, :, 5]))
                            c = cosine_similarity(vcat(data1[i, :, 3], data1[i, :, 6]),
                                                  vcat(data1[i + intraval, :, 3], data1[i + intraval, :, 6]))
                            cos_sim[i, :] = [a, b, c]
                        end

                        # Compute distances
                        dist_a = hcat(norm.(cos_sim[st:end, 1] .- cos_sim[st:end, 2]),
                                      norm.(cos_sim[st:end, 1] .- cos_sim[st:end, 3]))
                        dist_b = hcat(norm.(cos_sim[st:end, 1] .- cos_sim[st:end, 2]),
                                      norm.(cos_sim[st:end, 2] .- cos_sim[st:end, 3]))
                        dist_c = hcat(norm.(cos_sim[st:end, 1] .- cos_sim[st:end, 3]),
                                      norm.(cos_sim[st:end, 2] .- cos_sim[st:end, 3]))

                        # Find md
                        mean_dist_a = mean(dist_a, dims=2)
                        mean_dist_b = mean(dist_b, dims=2)
                        mean_dist_c = mean(dist_c, dims=2)

                        md_candidates = [
                            findfirst(x -> x > 0.12, vec(mean_dist_a)),
                            findfirst(x -> x > 0.12, vec(mean_dist_b)),
                            findfirst(x -> x > 0.12, vec(mean_dist_c))
                        ]
                        md_candidates = filter(!isnothing, md_candidates)
                        if isempty(md_candidates)
                            md = st
                        else
                            md = minimum(md_candidates)
                        end

                        # Compute thresholds
                        th = zeros(length(j_values))
                        for (k, dth) in enumerate(j_values)
                            success = 0
                            for species in 1:3
                                cond1 = (count(<(dth), cos_sim[md:end, species]) != 0) && (count(==(0), NumS[2:end, species]) != 0)
                                cond2 = (count(<(dth), cos_sim[md:end, species]) == 0) && (count(==(0), NumS[2:end, species]) == 0)
                                if cond1 || cond2
                                    success += 1
                                end
                            end
                            if success == 3
                                th[k] = dth
                            end
                        end
                        max_th = maximum(th)
                        @printf("rn = %d, intraval = %d, max_th = %g\n", rn, intraval, max_th)
                        mmat[ii, rn + 1] = max_th
                    else
                        @warn "Dataset not found for path: $histogram_path"
                    end
                finally
                    close(h5file)
                end
            end
        end
    end

    # Save mmat to HDF5 if needed
    # h5write("/Volumes/yoondata/ERPS_intra/ERPS_intra_test.h5", "/CSim", mmat)
end

# Second code block
function compute_S_probability()
    inte = 20:20
    j_values = 0:0.01:1
    intra2 = 1.0
    intra3 = 0.5

    for intra1_factor in 1:5
        intra1 = intra1_factor * 0.2
        st = 200
        mmat = zeros(Float64, 100, 100 + 1)
        S_probability = zeros(Float64, 100, 101, 2)

        for ii in 10:length(inte)
            mmat[ii, 1] = inte[ii]
            intraval = inte[ii]

            @threads for rn in 1:100
                # File paths and dataset names
                intra1_str = @sprintf("%.3f", intra1)
                intra2_str = @sprintf("%.3f", intra2)
                intra3_str = @sprintf("%.3f", intra3)
                h5_filename = "/home/ty/Desktop/yoonD/RPS/intra/RPS_intra.h5"
                histogram_path = "/intra_$(intra1_str)_$(intra2_str)_$(intra3_str)/Histogram/$(rn)"
                nums_path = "/intra_$(intra1_str)_$(intra2_str)_$(intra3_str)/NumS/$(rn)"

                # Open HDF5 file once per thread
                h5file = h5open(h5_filename, "r")
                try
                    data = read(h5file[histogram_path])
                    data1 = data[2:end, :, :]
                    data2 = read(h5file[nums_path])
                    NumS = data2[2:end, :]
                    L = size(NumS, 1)
                    cos_sim = zeros(L - intraval, 3)

                    # Compute cosine similarities
                    for i in 1:(L - intraval)
                        a = cosine_similarity(vcat(data1[i, :, 1], data1[i, :, 2]),
                                              vcat(data1[i + intraval, :, 1], data1[i + intraval, :, 2]))
                        b = cosine_similarity(vcat(data1[i, :, 1], data1[i, :, 3]),
                                              vcat(data1[i + intraval, :, 1], data1[i + intraval, :, 3]))
                        c = cosine_similarity(vcat(data1[i, :, 1], data1[i, :, 4]),
                                              vcat(data1[i + intraval, :, 1], data1[i + intraval, :, 4]))
                        cos_sim[i, :] = [a, b, c]
                    end

                    # Compute distances
                    dist_a = hcat(norm.(cos_sim[st:end, 1] .- cos_sim[st:end, 2]),
                                  norm.(cos_sim[st:end, 1] .- cos_sim[st:end, 3]))
                    dist_b = hcat(norm.(cos_sim[st:end, 1] .- cos_sim[st:end, 2]),
                                  norm.(cos_sim[st:end, 2] .- cos_sim[st:end, 3]))
                    dist_c = hcat(norm.(cos_sim[st:end, 1] .- cos_sim[st:end, 3]),
                                  norm.(cos_sim[st:end, 2] .- cos_sim[st:end, 3]))

                    # Find md
                    mean_dist_a = mean(dist_a, dims=2)
                    mean_dist_b = mean(dist_b, dims=2)
                    mean_dist_c = mean(dist_c, dims=2)

                    md_candidates = [
                        findfirst(x -> x > 0.10, vec(mean_dist_a)),
                        findfirst(x -> x > 0.10, vec(mean_dist_b)),
                        findfirst(x -> x > 0.10, vec(mean_dist_c))
                    ]
                    md_candidates = filter(!isnothing, md_candidates)
                    if isempty(md_candidates)
                        md = st
                    else
                        md = minimum(md_candidates)
                    end

                    # Compute thresholds
                    th = zeros(length(j_values))
                    for (k, dth) in enumerate(j_values)
                        success_e = 0
                        success_s = 0
                        for species in 1:3
                            cond_e = (count(<(dth), cos_sim[md:end, species]) != 0) && (count(==(0), NumS[2:end, species]) != 0)
                            cond_s = (count(<(dth), cos_sim[md:end, species]) == 0) && (count(==(0), NumS[2:end, species]) == 0)
                            if cond_e
                                success_e += 1
                            elseif cond_s
                                success_s += 1
                            end
                        end
                        S_probability[rn, k, :] = [success_e, success_s]
                        if success_e + success_s == 3
                            th[k] = dth
                        end
                    end
                    max_th = maximum(th)
                    @printf("rn = %d, intraval = %d, max_th = %g\n", rn, intraval, max_th)
                    mmat[ii, rn + 1] = max_th
                finally
                    close(h5file)
                end
            end
        end

        # Save S_probability to HDF5
        dataset_path = "/intra_$(intra1_str)_$(intra2_str)_$(intra3_str)/P"
        h5write(h5_filename, dataset_path, S_probability)
    end
end

# Plotting code
function plot_S_probability()
    j_values = 0:0.01:1
    ii_values = 0.1:0.1:1.5
    intra2 = 1.0
    intra3 = 0.5

    plot()
    for (i, intra1) in enumerate(ii_values)
        intra1_str = @sprintf("%.3f", intra1)
        intra2_str = @sprintf("%.3f", intra2)
        intra3_str = @sprintf("%.3f", intra3)
        # h5_filename = "/Volumes/yoondata/RPS/intra/RPS_intra.h5"  # Fixed typo in path
        h5_filename = "/home/ty/Desktop/yoonD/RPS/intra/RPS_intra.h5"
        dataset_path = "/intra_$(intra1_str)_$(intra2_str)_$(intra3_str)/P"
        S_probability = h5read(h5_filename, dataset_path)
        S_P = zeros(length(j_values), 3)

        for (jte, j_val) in enumerate(j_values)
            S_P[jte, :] = [
                j_val,
                sum(S_probability[:, jte, 1]) / 200,
                sum(S_probability[:, jte, 2]) / 100
            ]
        end

        plot!(S_P[:, 1], S_P[:, 2] .* S_P[:, 3], label="intra1=$(intra1)")
    end

    xlabel!("j")
    ylabel!("S_P")
    title!("Plot of S_P vs j")
    legend()
    display(current())
end

# Run the functions
compute_mmat()
compute_S_probability()
plot_S_probability()