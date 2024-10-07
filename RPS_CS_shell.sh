#!/bin/bash

# Julia 함수를 위한 매개변수 설정
Lsize=200
reproduction_rate=2.0
selection_rate=2.0
mobility=30
intra1_start=0.2
intra1_end=0.2
intra2_start=0.2
intra2_end=0.2
intra3_start=0.2
intra3_end=0.2
ext=2
para=6.25
start_rn=1
end_rn=100

echo "Julia 실행 중..."
julia -t 10 RPS_CS_test.jl $Lsize $reproduction_rate $selection_rate $mobility $intra1_start $intra1_end $intra2_start $intra2_end $intra3_start $intra3_end $ext $para $start_rn $end_rn

# Julia가 성공적으로 완료되었는지 확인
if [ $? -ne 0 ]; then
    echo "Julia 실행이 실패했습니다. 중지합니다."
    exit 1
fi

echo "모든 Julia 작업이 성공적으로 완료되었습니다!"