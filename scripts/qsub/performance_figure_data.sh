## TEST
# python ../performance_comparison.py --qsub --models dtwf --walltime 20:00 --nprocs 1 --Ne 50 --sample_size 30 --max_chroms 2 --replicates 2 --outfile ~/nelson_projects/wf_coalescent/results/performance/test_performance.npz


## Ne 10000, samples 1000 - ** Only output paths matter, not file names **
python ../performance_comparison.py --qsub --models dtwf --walltime 30:00:00 --nprocs 1 --Ne 10000 --sample_size 1000 --max_chroms 22 --replicates 10 --outfile ~/nelson_projects/wf_coalescent/results/performance/dtwf_Ne10000_samples1000.npz
sleep 1

python ../performance_comparison.py --qsub --models hudson --walltime 30:00:00 --nprocs 1 --Ne 10000 --sample_size 1000 --max_chroms 22 --replicates 10 --outfile ~/nelson_projects/wf_coalescent/results/performance/hudson_Ne10000_samples1000.npz
sleep 1

python ../performance_comparison.py --qsub --models hybrid --walltime 30:00:00 --nprocs 1 --Ne 10000 --sample_size 1000 --max_chroms 22 --replicates 10 --outfile ~/nelson_projects/wf_coalescent/results/performance/hybrid100_Ne10000_samples1000.npz --hybrid_wf_gens 100
sleep 1

python ../performance_comparison.py --qsub --models hybrid --walltime 30:00:00 --nprocs 1 --Ne 10000 --sample_size 1000 --max_chroms 22 --replicates 10 --outfile ~/nelson_projects/wf_coalescent/results/performance/hybrid1000_Ne10000_samples1000.npz --hybrid_wf_gens 1000
sleep 1


## Ne 10000, samples 10000 - ** Only output paths matter, not file names **
python ../performance_comparison.py --qsub --models dtwf --walltime 48:00:00 --nprocs 2 --Ne 10000 --sample_size 10000 --max_chroms 22 --replicates 10 --outfile ~/nelson_projects/wf_coalescent/results/performance/dtwf_Ne10000_samples1000.npz
sleep 1

python ../performance_comparison.py --qsub --models hudson --walltime 48:00:00 --nprocs 2 --Ne 10000 --sample_size 10000 --max_chroms 22 --replicates 10 --outfile ~/nelson_projects/wf_coalescent/results/performance/hudson_Ne10000_samples1000.npz
sleep 1

python ../performance_comparison.py --qsub --models hybrid --walltime 48:00:00 --nprocs 2 --Ne 10000 --sample_size 10000 --max_chroms 22 --replicates 10 --outfile ~/nelson_projects/wf_coalescent/results/performance/hybrid100_Ne10000_samples1000.npz --hybrid_wf_gens 100
sleep 1

python ../performance_comparison.py --qsub --models hybrid --walltime 48:00:00 --nprocs 2 --Ne 10000 --sample_size 10000 --max_chroms 22 --replicates 10 --outfile ~/nelson_projects/wf_coalescent/results/performance/hybrid1000_Ne10000_samples1000.npz --hybrid_wf_gens 1000
sleep 1

python ../performance_comparison.py --qsub --models hybrid --walltime 48:00:00 --nprocs 2 --Ne 10000 --sample_size 10000 --max_chroms 22 --replicates 10 --outfile ~/nelson_projects/wf_coalescent/results/performance/hybrid1000_Ne10000_samples1000.npz --hybrid_wf_gens 10000
sleep 1


## Ne 500, samples 500 - ** Only output paths matter, not file names **
python ../performance_comparison.py --qsub --models dtwf --walltime 12:00:00 --nprocs 1 --Ne 500 --sample_size 500 --max_chroms 22 --replicates 10 --outfile ~/nelson_projects/wf_coalescent/results/performance/dtwf_Ne500_samples500.npz
sleep 1

python ../performance_comparison.py --qsub --models hudson --walltime 12:00:00 --nprocs 1 --Ne 500 --sample_size 500 --max_chroms 22 --replicates 10 --outfile ~/nelson_projects/wf_coalescent/results/performance/hudson_Ne500_samples500.npz
sleep 1

python ../performance_comparison.py --qsub --models hybrid --walltime 12:00:00 --nprocs 1 --Ne 500 --sample_size 500 --max_chroms 22 --replicates 10 --outfile ~/nelson_projects/wf_coalescent/results/performance/hybrid1000_Ne500_samples500.npz --hybrid_wf_gens 1000
sleep 1

python ../performance_comparison.py --qsub --models hybrid --walltime 12:00:00 --nprocs 1 --Ne 500 --sample_size 500 --max_chroms 22 --replicates 10 --outfile ~/nelson_projects/wf_coalescent/results/performance/hybrid100_Ne500_samples500.npz --hybrid_wf_gens 100
sleep 1
