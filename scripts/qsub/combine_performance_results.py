import sys, os
sys.path.append(os.path.expanduser('~/project/wf_coalescent/scripts/'))
import performance_comparison as pc
import numpy as np
from pathlib import Path

results_dir = os.path.expanduser(sys.argv[1])
key = sys.argv[2]

out_files = []
for entry in Path(results_dir).glob(key + '*.npz'):
    out_files.append(entry)

last_length = np.array([-1])
results = {}
for f in out_files:
    print("Loading:", f)
    loaded = np.load(f)
    if loaded['lengths'].shape[0] != 22:
        print("Incomplete - skipping")
        continue

    ## Check all loaded results are for the same chrom lengths
    if -1 in last_length:
        last_length = loaded['lengths']
    if len(set(last_length).intersection(set(loaded['lengths']))) < len(last_length) or \
            len(last_length) != len(loaded['lengths']):
        print("Length mismatch!")
        print("Loaded:", last_length)
        print("Discarding", loaded['lengths'])
        continue

    ## Add new keys or append to existing results
    for key, value in loaded.items():
        if key in ['num_chroms', 'lengths', 'args']:
            results[key] = value
            print("Adding", key)
            continue

        if key in results:
            x = results[key]
            results[key] = np.hstack([x, value])
            print("Appending", key)
        else:
            results[key] = value

# import IPython; IPython.embed()
    
# last_length = np.array([-1])
# results_100000 = {}
# for f in out_files_100000:
#     print("Loading:", f)
#     loaded = np.load(f)
#     if last_length.any() == -1:
#         last_length = loaded['lengths']
#     if last_length.all() != loaded['lengths'].all():
#         print("Length mismatch!", last_length, loaded['lengths'])
#         sys.exit()
#     results_100000.update(loaded)

# plotfile_100000 = os.path.expanduser('~/temp/perf_100000.pdf')
# pc.plot_times(plotfile_100000, **results_100000)
# plotfile_500 = '/home/dnelson/temp/perf_500.png'
# pc.plot_times(plotfile_500, **results_500)

plotfile_10000 = os.path.expanduser('~/temp/perf_10000.pdf')

try:
    pc.plot_times(plotfile_10000, **results)
except:
    import IPython; IPython.embed()
    raise


