import sys, os
sys.path.append("../../msprime")
sys.path.append("../../msprime/lib/subprojects/git-submodules/tskit/python/")
import msprime
import numpy as np
import scipy.sparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import argparse


def expected_num_segs(t, K, L):
    return (K + 2 * L * t) / 2 ** (2 * t - 1)


def expected_total_length(t, L):
    return L / 2 ** (2 * t - 1)


def ibd_list_to_df(ibd_list, ca_times=None):
    cols = ["ind1", "ind2", "start", "end"]

    df = pd.DataFrame(np.array(ibd_list), columns=cols)
    df['len'] = df['end'] - df['start']

    return df


def plot_expected_ibd(max_ca_time, length_in_morgans, ax, K=22):
    ## Get theory points
    K = 22
    L = length_in_morgans
    times = range(1, max_ca_time + 1)

    expected_x = [expected_total_length(t, L) * 1e8 for t in times]
    expected_y = [expected_num_segs(t, K, L) for t in times]

    ax.scatter(expected_x, expected_y, s=20, c='black', label='Expected')


def plot_ibd_df(df, ca_times=None, min_length=1e6, ax=None, max_ca_time=10):
    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 4))

    df = df[df['len'] >= min_length]

    df = df.assign(count=np.ones(df.shape[0]))
    warned = False
    for ind in set(df['ind1'].values):
        ind_df = df[df['ind1'] == ind][['ind2', 'len', 'count']]
        df2 = ind_df.groupby('ind2').sum()
        df2 = df2.reset_index()

        df2 = df2.assign(colours=np.zeros(df2.shape[0]))
        df2['tmrca'] = 1e-9
        if ca_times is not None:
            get_ca_time = lambda x: ca_times[sorted([ind, x])[0], sorted([ind, x])[1]]
            df2 = df2.assign(
                    tmrca=df2['ind2'].apply(get_ca_time), axis=1
                    )

        df2 = df2[df2['tmrca'] <= max_ca_time]

        ## We can have TMRCA of 0 if max_ca_time is different than the value used
        ## to generate the IBD array
        if df2['tmrca'].min() == 0:
            if warned is False:
                print("Warning - filtering TMRCA >", max_ca_time)
                warned = True
            df2 = df2[df2['tmrca'] > 0]

        num_segments = df2['count'].values
        total_IBD = df2['len'].values
        colours = df2['tmrca'].values

        ax.scatter(total_IBD, num_segments, c=colours, s=5,
                cmap='viridis', vmin=0, vmax=max_ca_time)

    ax.set_ylabel('Number of IBD segments')
    ax.set_xlabel('')
    if ax.is_last_row():
        ax.set_xlabel('Total IBD (base pairs)')
    ax.set_xscale('log')
    ax.set_xlim((0.7 * min_length, 4e9))

    return ax


def main(args):
    ts_file = os.path.expanduser(args.dtwf_ts_file)
    ts = msprime.load(ts_file)

    ## Set up plot axes
    fig, ax_arr = plt.subplots(3, 1, figsize=(7, 7), sharex=True, sharey=True)

    ## Load Genizon IBD and plot
    ibd_file = os.path.expanduser(args.genizon_ibd_file)
    loaded = np.load(ibd_file)
    ibd_array = loaded['ibd_array']
    ibd_df = ibd_list_to_df(ibd_array)
    # nrows = 100000
    # ibd_df = pd.read_csv(ibd_file, nrows=nrows, delim_whitespace=True,
    #         header=None, usecols=(1, 3, 5, 6, 10),
    #         names=("ind1", "ind2", "start", "end", "len"))
    # ibd_df['len'] = ibd_df['len'] * 1e6 # Convert cM to base pairs
    plot_ibd_df(ibd_df, ax=ax_arr[0])
    ax_arr[0].set_title('Genizon Data')

    ## Get DTWF common ancestor times for IBD segments
    print("Loading DTWF")
    ca_times_dtwf = scipy.sparse.load_npz(
                os.path.expanduser(args.dtwf_ca_file)).tolil()

    ## Load DTWF IBD and plot
    ibd_file = os.path.expanduser(args.dtwf_ibd_file)
    loaded = np.load(ibd_file)
    ibd_array = loaded['ibd_array']
    ibd_df = ibd_list_to_df(ibd_array)
    plot_ibd_df(ibd_df, ca_times_dtwf, ax=ax_arr[1])
    ax_arr[1].set_title('msprime (WF)')

    ## Get Hudson common ancestor times for IBD segments
    print("Loading Hudson")
    ts_file = os.path.expanduser(args.hudson_ts_file)
    ts = msprime.load(ts_file)
    ca_times_hudson = scipy.sparse.load_npz(
            os.path.expanduser(args.hudson_ca_file)).tolil()

    ## Load Hudson IBD and plot
    ibd_file = os.path.expanduser(args.hudson_ibd_file)
    loaded = np.load(ibd_file)
    ibd_array = loaded['ibd_array']
    ibd_df = ibd_list_to_df(ibd_array)
    plot_ibd_df(ibd_df, ca_times_hudson, ax=ax_arr[2])
    ax_arr[2].set_title('msprime (Hudson)')

    ## Plot expected IBD cluster means
    max_theory_ca_time = 5
    K = 22
    length_in_morgans = ts.get_sequence_length() / 1e8
    plot_expected_ibd(max_theory_ca_time, length_in_morgans, ax=ax_arr[0], K=K)
    plot_expected_ibd(max_theory_ca_time, length_in_morgans, ax=ax_arr[1], K=K)
    plot_expected_ibd(max_theory_ca_time, length_in_morgans, ax=ax_arr[2], K=K)

    ## Legend only on upper plot
    ax_arr[0].legend()

    ## Jump through a few hoops to set colourbars
    sm = plt.cm.ScalarMappable(cmap='viridis',
            norm=plt.Normalize(vmin=0, vmax=args.max_ibd_time_gens))
    sm._A = []
    cbar = fig.colorbar(sm, aspect=50, ax=ax_arr.ravel().tolist())
    cbar.ax.get_yaxis().labelpad = 5
    cbar.ax.set_ylabel('TMRCA (generations)', rotation=90)

    fig.savefig(os.path.expanduser(args.outfile))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--dtwf_ibd_file", required=True)
    parser.add_argument("--dtwf_ca_file", required=True)
    parser.add_argument("--dtwf_ts_file", required=True)
    parser.add_argument("--hudson_ibd_file", required=True)
    parser.add_argument("--hudson_ca_file", required=True)
    parser.add_argument("--hudson_ts_file", required=True)
    parser.add_argument("--genizon_ibd_file", required=True)
    parser.add_argument("--outfile", required=True)
    parser.add_argument("--max_ibd_time_gens", type=int, default=10)
    args = parser.parse_args()

    main(args)
