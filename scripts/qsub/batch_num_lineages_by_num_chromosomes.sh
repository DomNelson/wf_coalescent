#!/bin/bash
set -e

SAMPLE_SIZE=10000
NE=10000

for NC in $( seq 1 22 ); do python num_lineages.py --srun --time 6:00:00 --mem 7000 --nchroms $NC --Ne $NE --sample_size $SAMPLE_SIZE --out_dir ~/dnelson_projects/wf_coalescent/results/ --model dtwf; sleep 1; done

for NC in $( seq 1 22 ); do python num_lineages.py --srun --time 6:00:00 --mem 7000 --nchroms $NC --Ne $NE --sample_size $SAMPLE_SIZE --out_dir ~/dnelson_projects/wf_coalescent/results/ --model hudson; sleep 1; done
