#!/usr/bin/env bash

set -euo pipefail

RUNS=20

BIN_FFTW3="./build/bench_fftw3"
BIN_PFFFT="./build/bench_pffft"

FFT_SIZES=(512 1024 2048 4096)

echo "==================================================="
echo " Benchmark Suite ($RUNS runs per FFT size)"
echo "==================================================="

run_and_parse() {
    local binary="$1"
    local name="$2"

    echo
    echo "===== $name ====="

    for size in "${FFT_SIZES[@]}"; do
        local total=0

        printf "FFT Size %d " "$size"

        for ((i=1; i<=RUNS; i++)); do

            # Extract ONLY the matching FFT size result
            val=$(
                "$binary" |
                awk -v target="$size" '
                    /Running/ {
                        current_size = $(NF-1)
                    }

                    /Average Time per FFT:/ {
                        if (current_size == target) {
                            print $(NF-1)
                        }
                    }
                '
            )

            if [[ -z "$val" ]]; then
                echo
                echo "ERROR: Failed to parse benchmark output for size $size"
                exit 1
            fi

            total=$(echo "$total + $val" | bc -l)

            printf "."
        done

        avg=$(echo "scale=8; $total / $RUNS" | bc -l)

        printf "\n-> Stable Average (%d): %.5f us\n\n" "$size" "$avg"
    done
}

run_and_parse "$BIN_FFTW3" "FFTW3"
run_and_parse "$BIN_PFFFT" "PFFFT"

echo "==================================================="
echo " Benchmark complete."
echo "==================================================="
