using HDF5, LinearAlgebra, Random, Printf, Statistics

function cosine_similarity(A, B)
    dotA = sum(A.^2)
    dotB = sum(B.^2)
    return dot(A, B) / sqrt(dotA * dotB)
end

function RPS_simulation(inte, j, intra_range, intra2, intra3)
    st = 200  # 'st' 매개변수 추가

    for intra1_factor in intra_range
        intra1 = intra1_factor

        mmat = zeros(Float64, length(inte), 101)
        S_probability = zeros(Float64, 100, length(j), 2)

        for ii in eachindex(inte)
            mmat[ii, 1] = inte[ii]
            intraval = inte[ii]

            for rn in 1:100
                hdf5_file = "/Volumes/yoonD/RPS/intra/RPS_intra.h5"
                intra_group = "/intra_$(Printf.@sprintf("%.3f", intra1))_$(Printf.@sprintf("%.3f", intra2))_$(Printf.@sprintf("%.3f", intra3))"
                
                data = h5read(hdf5_file, "$(intra_group)/Histogram/$(rn)")
                data1 = data[2:end, :, :]
                data2 = h5read(hdf5_file, "$(intra_group)/NumS/$(rn)")
                NumS = data2[2:end, :]
                L = size(NumS, 1)
                cos_sim = zeros(Float64, L - intraval, 3)

                for i in 1:(L - intraval)
                    a = cosine_similarity([data1[i, :, 1]; data1[i, :, 4]], 
                                          [data1[i + intraval, :, 1]; data1[i + intraval, :, 4]])
                    b = cosine_similarity([data1[i, :, 2]; data1[i, :, 5]], 
                                          [data1[i + intraval, :, 2]; data1[i + intraval, :, 5]])
                    c = cosine_similarity([data1[i, :, 3]; data1[i, :, 6]], 
                                          [data1[i + intraval, :, 3]; data1[i + intraval, :, 6]])
                    cos_sim[i, :] = [a, b, c]  # 여기를 수정했습니다
                end

                dist_a = [norm(cos_sim[st:end, 1] - cos_sim[st:end, 2], 2), norm(cos_sim[st:end, 1] - cos_sim[st:end, 3], 2)]
                dist_b = [norm(cos_sim[st:end, 1] - cos_sim[st:end, 2], 2), norm(cos_sim[st:end, 1] - cos_sim[st:end, 3], 2)]
                dist_c = [norm(cos_sim[st:end, 1] - cos_sim[st:end, 3], 2), norm(cos_sim[st:end, 2] - cos_sim[st:end, 3], 2)]

                md_values = filter(!isnothing, [findfirst(x -> x > 0.1, mean(dist_a, dims=2)),
                                                findfirst(x -> x > 0.1, mean(dist_b, dims=2)),
                                                findfirst(x -> x > 0.1, mean(dist_c, dims=2))])
                
                md = isempty(md_values) ? nothing : minimum(md_values)

                if isnothing(md)
                    println("Warning: md is nothing for rn = $rn, intraval = $intraval")
                    continue
                end

                th = zeros(Float64, length(j))

                for k in eachindex(j)
                    dth = j[k]
                    success_e = 0
                    success_s = 0

                    for species in 1:3
                        if !isempty(findall(x -> x < dth, cos_sim[md:end, species])) && !isempty(findall(x -> x == 0, NumS[2:end, species]))
                            success_e += 1
                        elseif isempty(findall(x -> x < dth, cos_sim[md:end, species])) && isempty(findall(x -> x == 0, NumS[2:end, species]))
                            success_s += 1
                        end
                    end

                    S_probability[rn, k, :] = [success_e, success_s]

                    if success_e + success_s == 3
                        th[k] = dth
                    end
                end

                println("rn = $rn, intraval = $intraval, max_th = $(maximum(th))")
                mmat[ii, rn + 1] = maximum(th)
            end
        end

        h5open("/Volumes/yoonD/RPS/intra/RPS_intra.h5", "cw") do file
            group_name = "/intra_$(Printf.@sprintf("%.3f", intra1))_$(Printf.@sprintf("%.3f", intra2))_$(Printf.@sprintf("%.3f", intra3))"
            
            if haskey(file, "$group_name/P")
                delete_object(file, "$group_name/P")
            end
            write(file, "$group_name/P", S_probability)
            
            if haskey(file, "$group_name/mmat")
                delete_object(file, "$group_name/mmat")
            end
            write(file, "$group_name/mmat", mmat)
        end
    end
end

# 예제 매개변수
inte = 20#1:100
j = 0:0.01:1
intra_range = 1.1:0.1:1.2
intra2 = 1.000
intra3 = 0.500

# 시뮬레이션 실행
RPS_simulation(inte, j, intra_range, intra2, intra3)
