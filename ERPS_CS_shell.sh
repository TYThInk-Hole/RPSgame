#!/bin/bash

# 사용자로부터 모든 값을 한 번에 입력 받기
echo "n_threads start_rn end_rn chunk_size intra1 intra2 intra3_start intra3_end intra3_step intra4 intra5 ext para 값을 순서대로 입력하세요 (공백으로 구분):"
read n_threads start_rn end_rn chunk_size intra1 intra2 intra3_start intra3_end intra3_step intra4 intra5 ext para

if [ ! -f RPSLS_CS_input_log.txt ]; then
    touch ERPS_CS_input_log.txt
fi

echo "$n_threads $start_rn $end_rn $chunk_size $intra1 $intra2 $intra3_start $intra3_end $intra3_step $intra4 $intra5 $ext $para" >> ERPS_CS_input_log.txt
# Julia 함수를 위한 다른 매개변수 설정
Lsize=200
reproduction_rate=2.0
selection_rate=2.0
mobility=30
para=3.0

# 입력값 확인
echo "입력된 값:"
echo "start_rn: $start_rn"
echo "end_rn: $end_rn"
echo "chunk_size: $chunk_size"
echo "intra1: $intra1"
echo "intra2: $intra2"
echo "intra3_start: $intra3_start"
echo "intra3_end: $intra3_end"
echo "intra3_step: $intra3_step"
echo "intra4: $intra4"
echo "intra5: $intra5"
echo "n_threads: $n_threads"
echo "ext: $ext"

# rn 범위를 청크로 나누어 반복
for ((start=$start_rn; start<$end_rn; start+=$chunk_size))
do
    end=$((start + $chunk_size - 1))

    if ((end >= end_rn)); then
        end=$end_rn
    fi

    # rn 값의 청크로 Julia 스크립트 호출
    echo "rn = [$start, $end]에 대해 Julia 실행 중..."
    julia -t $n_threads ERPS_CS_test.jl $Lsize $reproduction_rate $selection_rate $mobility $intra1 $intra2 $intra3_start $intra3_end $intra3_step $intra4 $intra5 $ext $para $start $end

    # Julia가 성공적으로 완료되었는지 확인
    if [ $? -ne 0 ]; then
        echo "rn = [$start, $end]에 대한 Julia 실행이 실패했습니다. 중지합니다."
        exit 1
    fi
done

echo "모든 Julia 작업이 성공적으로 완료되었습니다!"