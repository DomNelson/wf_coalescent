# %matplotlib inline
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys, os
sys.path.append('/home/dnelson/project/wf_coalescent/scripts/')
sys.path.append('/home/dnelson/project/msprime/')
import numpy as np
import scipy.sparse
import argparse
from tqdm import tqdm


def plot_ibd_df(df, min_length=1, ax=None, max_ca_time=10):
    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 4))

    df = df[df['len'] >= min_length]

    df = df.assign(count=np.ones(df.shape[0]))
    warned = False
    for ind in tqdm(set(df['ind1'].values)):
        ind_df = df[df['ind1'] == ind][['ind2', 'len', 'count']]
        df2 = ind_df.groupby('ind2').sum()
        df2 = df2.reset_index()

        df2 = df2.assign(colours=np.zeros(df2.shape[0]))

        num_segments = df2['count'].values
        total_IBD = df2['len'].values * 1e6
#         print(num_segments.max())
#         colours = df2['tmrca'].values
        colours = np.zeros(df2.shape[0])

        ax.scatter(total_IBD, num_segments, c=colours, s=5,
                cmap='viridis', vmin=0, vmax=max_ca_time)

    ax.set_ylabel('Number of IBD segments')
    ax.set_xlabel('')
    if ax.is_last_row():
        ax.set_xlabel('Total IBD (base pairs)')
    ax.set_xscale('log')
#     ax.set_xlim((1, 4e3))

    return ax


def main(args):
    ibd_file = os.path.expanduser(args.ibd_file)
    plot_file = os.path.expanduser(args.plot_file)

    print("Loading IBD...")
    df = pd.read_csv(ibd_file, nrows=args.nrows, delim_whitespace=True, header=None,
		    usecols=(1, 3, 5, 6, 10), names=("ind1", "ind2", "start", "end", "len"))

    print("Plotting...")
    fig, ax = plt.subplots()
    ax_out = plot_ibd_df(df, min_length=args.min_length, ax=ax)
    fig.savefig(plot_file)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--ibd_file", required=True)
    parser.add_argument("--nrows", type=int)
    parser.add_argument("--min_length", default=1, type=float)
    parser.add_argument("--plot_file", required=True)

    args = parser.parse_args()
    main(args)
