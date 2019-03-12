python performance_comparison.py --qsub --models dtwf --walltime 2:00:00:00 --nprocs 3 --Ne 500 --sample_size 500 --max_chroms 22 --replicates 10 --outfile ~/nelson_projects/wf_coalescent/results/performance/dtwf_Ne500_samples500.npz

python performance_comparison.py --qsub --models dtwf --walltime 2:00:00:00 --nprocs 3 --Ne 10000 --sample_size 1000 --max_chroms 22 --replicates 10 --outfile ~/nelson_projects/wf_coalescent/results/performance/dtwf_Ne1000_samples1000.npz

python performance_comparison.py --qsub --models hudson --walltime 2:00:00:00 --nprocs 3 --Ne 500 --sample_size 500 --max_chroms 22 --replicates 10 --outfile ~/nelson_projects/wf_coalescent/results/performance/hudson_Ne500_samples500.npz

python performance_comparison.py --qsub --models hudson --walltime 2:00:00:00 --nprocs 3 --Ne 10000 --sample_size 1000 --max_chroms 22 --replicates 10 --outfile ~/nelson_projects/wf_coalescent/results/performance/hudson_Ne10000_samples1000.npz

python performance_comparison.py --qsub --models hybrid --walltime 2:00:00:00 --nprocs 3 --Ne 500 --sample_size 500 --max_chroms 22 --replicates 10 --outfile ~/nelson_projects/wf_coalescent/results/performance/hybrid1000_Ne500_samples500.npz --hybrid_wf_gens 1000

python performance_comparison.py --qsub --models hybrid --walltime 2:00:00:00 --nprocs 3 --Ne 10000 --sample_size 1000 --max_chroms 22 --replicates 10 --outfile ~/nelson_projects/wf_coalescent/results/performance/hybrid1000_Ne10000_samples1000.npz --hybrid_wf_gens 1000

python performance_comparison.py --qsub --models hybrid --walltime 2:00:00:00 --nprocs 3 --Ne 500 --sample_size 500 --max_chroms 22 --replicates 10 --outfile ~/nelson_projects/wf_coalescent/results/performance/hybrid100_Ne500_samples500.npz --hybrid_wf_gens 100

python performance_comparison.py --qsub --models hybrid --walltime 2:00:00:00 --nprocs 3 --Ne 10000 --sample_size 1000 --max_chroms 22 --replicates 10 --outfile ~/nelson_projects/wf_coalescent/results/performance/hybrid100_Ne10000_samples1000.npz --hybrid_wf_gens 100
