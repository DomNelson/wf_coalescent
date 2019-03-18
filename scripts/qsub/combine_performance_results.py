import sys, os
sys.path.append('/home/dnelson/project/wf_coalescent/scripts/')
import performance_comparison as pc
import numpy as np

results_dir = '/home/dnelson/project/wf_coalescent/results/performance/abacus/'

out_files_10000 = []
out_files_500 = []
for entry in list(os.scandir(results_dir)):
    if not entry.is_file():
        continue
        
    name, ext = os.path.splitext(entry.name)
    if ext == '.npz':
	## Crude way of extracting specific parameters
        if '500' in name:
            out_files_500.append(os.path.join(results_dir, entry.name))
        elif '10000' in name:
            out_files_10000.append(os.path.join(results_dir, entry.name))

results_10000 = {}
for f in out_files_10000:
    loaded = np.load(f)
    results_10000.update(loaded)
    
results_500 = {}
for f in out_files_500:
    loaded = np.load(f)
    results_500.update(loaded)

plotfile_500 = '/home/dnelson/temp/perf_500.pdf'
pc.plot_times(plotfile_500, **results_500)
plotfile_500 = '/home/dnelson/temp/perf_500.png'
pc.plot_times(plotfile_500, **results_500)

plotfile_10000 = '/home/dnelson/temp/perf_10000.pdf'
pc.plot_times(plotfile_10000, **results_10000)


