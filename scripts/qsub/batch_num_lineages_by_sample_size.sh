#!/bin/bash
set -e

# # DTWF sims
# python num_lineages.py --srun --time 24:00:00 --mem 7000 --nchroms 22 --Ne 10000 --sample_size 500 --out_dir ~/dnelson_projects/wf_coalescent/results/masked_num_lineages/ --model dtwf --growth_rate 0.001
# sleep 1
# python num_lineages.py --srun --time 24:00:00 --mem 7000 --nchroms 22 --Ne 10000 --sample_size 1000 --out_dir ~/dnelson_projects/wf_coalescent/results/masked_num_lineages/ --model dtwf --growth_rate 0.001
# sleep 1
# python num_lineages.py --srun --time 24:00:00 --mem 7000 --nchroms 22 --Ne 10000 --sample_size 2000 --out_dir ~/dnelson_projects/wf_coalescent/results/masked_num_lineages/ --model dtwf --growth_rate 0.001
# sleep 1
# python num_lineages.py --srun --time 24:00:00 --mem 7000 --nchroms 22 --Ne 10000 --sample_size 5000 --out_dir ~/dnelson_projects/wf_coalescent/results/masked_num_lineages/ --model dtwf --growth_rate 0.001
# sleep 1
# python num_lineages.py --srun --time 24:00:00 --mem 7000 --nchroms 22 --Ne 10000 --sample_size 10000 --out_dir ~/dnelson_projects/wf_coalescent/results/masked_num_lineages/ --model dtwf --growth_rate 0.001
# sleep 1
#
# # Hudson sims
# python num_lineages.py --srun --time 24:00:00 --mem 7000 --nchroms 22 --Ne 10000 --sample_size 500 --out_dir ~/dnelson_projects/wf_coalescent/results/masked_num_lineages/ --model hudson --growth_rate 0.001
# sleep 1
# python num_lineages.py --srun --time 24:00:00 --mem 7000 --nchroms 22 --Ne 10000 --sample_size 1000 --out_dir ~/dnelson_projects/wf_coalescent/results/masked_num_lineages/ --model hudson --growth_rate 0.001
# sleep 1
# python num_lineages.py --srun --time 24:00:00 --mem 7000 --nchroms 22 --Ne 10000 --sample_size 2000 --out_dir ~/dnelson_projects/wf_coalescent/results/masked_num_lineages/ --model hudson --growth_rate 0.001
# sleep 1
# python num_lineages.py --srun --time 24:00:00 --mem 7000 --nchroms 22 --Ne 10000 --sample_size 5000 --out_dir ~/dnelson_projects/wf_coalescent/results/masked_num_lineages/ --model hudson --growth_rate 0.001
# sleep 1
# python num_lineages.py --srun --time 24:00:00 --mem 7000 --nchroms 22 --Ne 10000 --sample_size 10000 --out_dir ~/dnelson_projects/wf_coalescent/results/masked_num_lineages/ --model hudson --growth_rate 0.001

# Hybrid sims - 100 generations
python num_lineages.py --srun --time 24:00:00 --mem 7000 --nchroms 22 --Ne 10000 --sample_size 500 --out_dir ~/dnelson_projects/wf_coalescent/results/masked_num_lineages/ --model hybrid --growth_rate 0.001 --hybrid_switch_time 100
sleep 1
python num_lineages.py --srun --time 24:00:00 --mem 7000 --nchroms 22 --Ne 10000 --sample_size 1000 --out_dir ~/dnelson_projects/wf_coalescent/results/masked_num_lineages/ --model hybrid --growth_rate 0.001 --hybrid_switch_time 100
sleep 1
python num_lineages.py --srun --time 24:00:00 --mem 7000 --nchroms 22 --Ne 10000 --sample_size 2000 --out_dir ~/dnelson_projects/wf_coalescent/results/masked_num_lineages/ --model hybrid --growth_rate 0.001 --hybrid_switch_time 100
sleep 1
python num_lineages.py --srun --time 24:00:00 --mem 7000 --nchroms 22 --Ne 10000 --sample_size 5000 --out_dir ~/dnelson_projects/wf_coalescent/results/masked_num_lineages/ --model hybrid --growth_rate 0.001 --hybrid_switch_time 100
sleep 1
python num_lineages.py --srun --time 24:00:00 --mem 7000 --nchroms 22 --Ne 10000 --sample_size 10000 --out_dir ~/dnelson_projects/wf_coalescent/results/masked_num_lineages/ --model hybrid --growth_rate 0.001 --hybrid_switch_time 100
sleep 1

# Hybrid sims - 1000 generations
python num_lineages.py --srun --time 24:00:00 --mem 7000 --nchroms 22 --Ne 10000 --sample_size 500 --out_dir ~/dnelson_projects/wf_coalescent/results/masked_num_lineages/ --model hybrid --growth_rate 0.001 --hybrid_switch_time 1000
sleep 1
python num_lineages.py --srun --time 24:00:00 --mem 7000 --nchroms 22 --Ne 10000 --sample_size 1000 --out_dir ~/dnelson_projects/wf_coalescent/results/masked_num_lineages/ --model hybrid --growth_rate 0.001 --hybrid_switch_time 1000
sleep 1
python num_lineages.py --srun --time 24:00:00 --mem 7000 --nchroms 22 --Ne 10000 --sample_size 2000 --out_dir ~/dnelson_projects/wf_coalescent/results/masked_num_lineages/ --model hybrid --growth_rate 0.001 --hybrid_switch_time 1000
sleep 1
python num_lineages.py --srun --time 24:00:00 --mem 7000 --nchroms 22 --Ne 10000 --sample_size 5000 --out_dir ~/dnelson_projects/wf_coalescent/results/masked_num_lineages/ --model hybrid --growth_rate 0.001 --hybrid_switch_time 1000
sleep 1
python num_lineages.py --srun --time 24:00:00 --mem 7000 --nchroms 22 --Ne 10000 --sample_size 10000 --out_dir ~/dnelson_projects/wf_coalescent/results/masked_num_lineages/ --model hybrid --growth_rate 0.001 --hybrid_switch_time 1000
