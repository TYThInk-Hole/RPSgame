#!/bin/bash

# rn 범위와 청크 크기 설정
start_rn=1
end_rn=1000
chunk_size=50

# Julia 함수를 위한 다른 매개변수 설정
# 이 매개변수들은 명령줄을 통해 전달되거나 여기에 하드코딩될 수 있습니다.
Lsize=200
reproduction_rate=2.0
selection_rate=2.0
mobility=30
intra1=1.9
intra2=2.0
intra3=6.5
intra4=1.3
intra5=0.7
ext=2
para=6.25

# rn 범위를 청크로 나누어 반복
for ((start=$start_rn; start<$end_rn; start+=$chunk_size))
do
    end=$((start + $chunk_size - 1))

    if ((end >= end_rn)); then
        end=$((end_rn - 1))
    fi

    # rn 값의 청크로 Julia 스크립트 호출
    echo "rn = [$start, $end]에 대해 Julia 실행 중..."
    julia -t 50 ERPS_CS_test.jl $Lsize $reproduction_rate $selection_rate $mobility $intra1 $intra2 $intra3 $intra4 $intra5 $ext $para $start $end

    # Julia가 성공적으로 완료되었는지 확인
    if [ $? -ne 0 ]; then
        echo "rn = [$start, $end]에 대한 Julia 실행이 실패했습니다. 중지합니다."
        exit 1
    fi
done

echo "모든 Julia 작업이 성공적으로 완료되었습니다!"