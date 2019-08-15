import numpy as np
import sys, os
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import seaborn

if len(sys.argv) != 2:
    print("Usage: python num_segs_dist.py PLOT_OUTFILE")
    sys.exit()

outfile = os.path.expanduser(sys.argv[1])

n_pairs = 50000
hs = [sum(np.random.normal(loc=2.135, scale=0.942) for i in range(22)) for i in range(n_pairs)]
hc = [sum(np.random.normal(loc=0.941, scale=1.020) for i in range(22)) for i in range(n_pairs)]
h2c = [sum(np.random.normal(loc=0.336, scale=0.445) for i in range(22)) for i in range(n_pairs)]

fig, ax = plt.subplots(1, 1)

counts, bin_edges = np.histogram(hs, bins=range(60))
counts = counts / n_pairs
ax.plot(bin_edges[:-1], counts, label="HS")
ax.fill_between(bin_edges[:-1], counts, alpha=0.5)

counts, bin_edges = np.histogram(hc, bins=range(60))
counts = counts / n_pairs
ax.plot(bin_edges[:-1], counts, label="HC")
ax.fill_between(bin_edges[:-1], counts, alpha=0.5)

counts, bin_edges = np.histogram(h2c, bins=range(60))
counts = counts / n_pairs
ax.plot(bin_edges[:-1], counts, label="H2C")
ax.fill_between(bin_edges[:-1], counts, alpha=0.5)


# out = ax.hist(hs, bins=range(60), alpha=0.5)
# out = ax.hist(hc, bins=range(60), alpha=0.5)
# out = ax.hist(h2c, bins=range(60), alpha=0.5)

fig.legend(bbox_to_anchor=(0.9, 0.85))
ax.set_xlabel("Number of shared segments")
ax.set_ylabel("Simulated probability")
ax.set_title("Simulated distribution of number of segments given\npedigree relationship")

fig.savefig(outfile)
