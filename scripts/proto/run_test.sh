#!/bin/bash
# usage: run_test.sh NAME KIND "gen args" "namelist extras" [analysis args]
# BIN is the spmm_bond executable (default: spmm_bond on PATH); runs go in scripts/proto/runs/
P=$(dirname "$(readlink -f "$0")"); name=$1; kind=$2; gen=$3; nml=$4; shift 4
mkdir -p $P/runs/$name && cd $P/runs/$name || exit 1
python3 $P/gen_shapes.py $kind shape.dat $gen > gen.log 2>&1 || { echo "gen failed"; exit 1; }
printf '&spmm\n do_melt=.false., check=.true., prefix="t",\n %s\n/\n' "$nml" > spmm.nml
${BIN:-spmm_bond} > log.txt 2>&1 || { echo "RUN FAILED: $(tail -2 log.txt | cut -c1-200)"; exit 1; }
echo "  [$name] $(tail -1 log.txt | sed 's/  */ /g' | cut -c1-110) invfail=$(grep -c INVARIANT log.txt)"
echo "  $(python3 $P/analyze.py $kind t_init.out t_final.out "$@" 2>&1 | tail -1)"
