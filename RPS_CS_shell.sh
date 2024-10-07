#!/bin/bash

# rn 범위와 청크 크기 설정
start_rn=1
end_rn=100
chunk_size=10

# Julia 함수를 위한 다른 매개변수 설정
Lsize=200
reproduction_rate=2.0
selection_rate=2.0
mobility=30
intra1_start=0.2
intra1_end=0.8
intra2_start=0.2
intra2_end=0.8
intra3_start=0.2
intra3_end=0.8
ext=2
para=6.25

# rn 범위를 청크로 나누어 반복
for ((start=$start_rn; start<=$end_rn; start+=$chunk_size))
do
    end=$((start + $chunk_size - 1))

    if ((end > end_rn)); then
        end=$end_rn
    fi

    # rn 값의 청크로 Julia 스크립트 호출
    echo "rn = [$start, $end]에 대해 Julia 실행 중..."
    julia -t 10 RPS_CS_test.jl $Lsize $reproduction_rate $selection_rate $mobility $intra1_start $intra1_end $intra2_start $intra2_end $intra3_start $intra3_end $start $end

    # Julia가 성공적으로 완료되었는지 확인
    if [ $? -ne 0 ]; then
        echo "rn = [$start, $end]에 대한 Julia 실행이 실패했습니다. 중지합니다."
        exit 1
    fi
done

echo "모든 Julia 작업이 성공적으로 완료되었습니다!"