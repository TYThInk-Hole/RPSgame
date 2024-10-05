#!/bin/bash

# Set the total range for rn and the size of each chunk
start_rn=0
end_rn=1000
chunk_size=50

# Get the other parameters for the Julia function
# These parameters could be passed via the command line or hardcoded here.
Lsize=200
reproduction_rate=2.0
selection_rate=2.0
mobility=30
intra1=0.1
intra2=0.1
intra3=0.1
ext=1
para=1.0

# Loop over the rn range in chunks
for ((start=$start_rn; start<$end_rn; start+=$chunk_size))
do
    end=$((start + $chunk_size - 1))

    if ((end > end_rn)); then
        end=$end_rn
    fi

    # Call the Julia script with the chunk of rn values
    echo "Running Julia for rn = [$start, $end]..."
    julia --threads=auto RPS_CS_test.jl $Lsize $reproduction_rate $selection_rate $mobility $intra1 $intra2 $intra3 $ext $para $start $end

    # Check if Julia completed successfully
    if [ $? -ne 0 ]; then
        echo "Julia execution failed for rn = [$start, $end]. Stopping."
        exit 1
    fi
done

echo "All Julia jobs completed successfully!"