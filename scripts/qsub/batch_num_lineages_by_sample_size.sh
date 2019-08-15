#!/bin/bash
set -e

python num_lineages.py --srun --time 24:00:00 --mem 7000 --nchroms 22 --Ne 10000 --sample_size 500 --out_dir ~/dnelson_projects/wf_coalescent/results/ --model dtwf --growth_rate 0
sleep 1
python num_lineages.py --srun --time 24:00:00 --mem 7000 --nchroms 22 --Ne 10000 --sample_size 1000 --out_dir ~/dnelson_projects/wf_coalescent/results/ --model dtwf --growth_rate 0
sleep 1
python num_lineages.py --srun --time 24:00:00 --mem 7000 --nchroms 22 --Ne 10000 --sample_size 2000 --out_dir ~/dnelson_projects/wf_coalescent/results/ --model dtwf --growth_rate 0
sleep 1
python num_lineages.py --srun --time 24:00:00 --mem 7000 --nchroms 22 --Ne 10000 --sample_size 5000 --out_dir ~/dnelson_projects/wf_coalescent/results/ --model dtwf --growth_rate 0
sleep 1
python num_lineages.py --srun --time 24:00:00 --mem 7000 --nchroms 22 --Ne 10000 --sample_size 10000 --out_dir ~/dnelson_projects/wf_coalescent/results/ --model dtwf --growth_rate 0
sleep 1

python num_lineages.py --srun --time 24:00:00 --mem 7000 --nchroms 22 --Ne 10000 --sample_size 500 --out_dir ~/dnelson_projects/wf_coalescent/results/ --model hudson --growth_rate 0
sleep 1
python num_lineages.py --srun --time 24:00:00 --mem 7000 --nchroms 22 --Ne 10000 --sample_size 1000 --out_dir ~/dnelson_projects/wf_coalescent/results/ --model hudson --growth_rate 0
sleep 1
python num_lineages.py --srun --time 24:00:00 --mem 7000 --nchroms 22 --Ne 10000 --sample_size 2000 --out_dir ~/dnelson_projects/wf_coalescent/results/ --model hudson --growth_rate 0
sleep 1
python num_lineages.py --srun --time 24:00:00 --mem 7000 --nchroms 22 --Ne 10000 --sample_size 5000 --out_dir ~/dnelson_projects/wf_coalescent/results/ --model hudson --growth_rate 0
sleep 1
python num_lineages.py --srun --time 24:00:00 --mem 7000 --nchroms 22 --Ne 10000 --sample_size 10000 --out_dir ~/dnelson_projects/wf_coalescent/results/ --model hudson --growth_rate 0
