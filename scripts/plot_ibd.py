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
import functools
from tqdm import tqdm


def expected_num_segs(t, K, L, diploid=False):
    n = (K + 2 * L * t) / 2 ** (2 * t - 1)

    if diploid:
        n *= 2

    return n


def expected_total_length(t, L, diploid=False):
    length = L / 2 ** (2 * t - 1)

    if diploid:
        length *= 2

    return length


def ibd_list_to_df(ibd_list, diploid=False, germline=False):
    ibd_array = np.array(ibd_list)
    cols = ["ind1", "ind2", "start", "end"]
    df = pd.DataFrame(ibd_array, columns=cols)
    df['len'] = df['end'] - df['start']

    if diploid is True:

        if germline:
            contains_decimal = False
            print("Unphasing GERMLINE data")
            for i in tqdm(range(ibd_array.shape[0])):
                for j in range(2):
                    # Weak check that decimal copy numbers are present
                    ind = ibd_array[i, j]
                    if ind != int(ind):
                        contains_decimal = True

                    # Strip GERMLINE decimal copy numbers to form diploid IDs
                    ibd_array[i, j] = int(ind)

            if not contains_decimal:
                raise ValueError("GERMLINE decimal copy numbers not found!")

            ## Messy but we just rebuild the dataframe with the new diploid IDs
            df = pd.DataFrame(ibd_array, columns=cols)
            df['len'] = df['end'] - df['start']

        else:
            print("Merging simulated haploids into unphased diploid")
            inds = list(set(ibd_array[:, 0:2].ravel()))
            np.random.shuffle(inds)

            # Get even number of inds and pair randomly
            diploid_dict = {}
            extra_haploids = len(inds) % 2
            for i in np.arange(0, len(inds) - extra_haploids, 2):
                ind1 = inds[i]
                ind2 = inds[i+1]
                assert ind1 == int(ind1)
                assert ind2 == int(ind2)
                ind1 = int(ind1)
                ind2 = int(ind2)

                diploid_ind = float("{}.{}".format(ind1, ind2))
                diploid_dict[ind1] = diploid_ind
                diploid_dict[ind2] = diploid_ind

            df[['ind1', 'ind2']] = df[['ind1', 'ind2']].replace(
                    to_replace=diploid_dict)

    # import IPython; IPython.embed()
    return df


def get_labels(times, diploid):
    labels_dict = {1: 'Half siblings', 2: 'Half-first cousins',
            3: 'Half-second cousins'}
    if diploid is True:
        labels_dict = {2: 'First cousins', 3: 'Second cousins'}
    
    labels = [labels_dict[t] for t in times if t in labels_dict]

    return labels


def plot_expected_ibd(max_ca_time, length_in_morgans, ax, K=22,
        diploid=False, avuncular=False, grandparental=False, min_ca_time=1,
        label_points=False):
    ## Get theory points
    K = 22
    L = length_in_morgans
    times = range(min_ca_time, max_ca_time + 1)

    expected_x = [expected_total_length(t, L, diploid) * 1e8 for t in times]
    expected_y = [expected_num_segs(t, K, L, diploid) for t in times]

    labels = get_labels(times, diploid)
    for l, x, y in zip(labels, expected_x, expected_y):
        if label_points:
            ax.annotate(l, xy=(x, y), xytext=(-25, 25),
                    textcoords='offset points', ha='right', va='bottom',
                    arrowprops=dict(arrowstyle = '-', connectionstyle='arc3,rad=0'),
                    fontsize=10)

    if avuncular:
        x = (L / 2) * 1e8
        y = (K + 3 * L) / 2
        if label_points:
            ax.annotate('Avuncular', xy=(x, y), xytext=(-30, 25),
                    textcoords='offset points', ha='right', va='bottom',
                    arrowprops=dict(arrowstyle = '-', connectionstyle='arc3,rad=0'))
        expected_x.append(x)
        expected_y.append(y)

    if grandparental:
        x = (L / 2) * 1e8
        y = (K + 2 * L) / 2
        if label_points:
            ax.annotate('Grandparent / grandchild', xy=(x, y), xytext=(35, -58),
                    textcoords='offset points', ha='right', va='bottom',
                    arrowprops=dict(arrowstyle = '-', connectionstyle='arc3,rad=0'))
        expected_x.append(x)
        expected_y.append(y)

    label = 'Expected (non-monogamous)'
    if diploid:
        label = 'Expected (monogamous)'

    ax.scatter(expected_x, expected_y, s=20, facecolors='white',
            edgecolors='black', label=label)


def get_tmrca(ind1, ind2, ca_times):
    if isinstance(ind1, str) and isinstance(ind2, str):
        tmrca = np.inf
        assert '_' in ind1 and '_' in ind2
        # Can't convert str representation of float to int directly
        ind1_copies = [int(float(x)) for x in ind1.split('_')]
        ind2_copies = [int(float(x)) for x in ind2.split('_')]

        for c1 in ind1_copies:
            for c2 in ind2_copies:
                idx1, idx2 = sorted([c1, c2])
                t = ca_times[idx1, idx2]
                if t > 0 and t < tmrca:
                    tmrca = ca_times[idx1, idx2]

    elif isinstance(ind1, float) and isinstance(ind2, float):
        tmrca = np.inf
        ind1 = str(ind1)
        ind2 = str(ind2)
        assert '.' in ind1 and '.' in ind2
        # Can't convert str representation of float to int directly
        ind1_copies = [int(float(x)) for x in ind1.split('.')]
        ind2_copies = [int(float(x)) for x in ind2.split('.')]

        for c1 in ind1_copies:
            for c2 in ind2_copies:
                idx1, idx2 = sorted([c1, c2])
                t = ca_times[idx1, idx2]
                if t > 0 and t < tmrca:
                    tmrca = ca_times[idx1, idx2]

    elif isinstance(ind1, int) and isinstance(ind2, int):
        idx1, idx2 = sorted([ind1, ind2])
        tmrca = ca_times[idx1, idx2]

    # if tmrca <= 0:
    #     import IPython; IPython.embed()

    return tmrca


def plot_ibd_df(df, ca_times=None, min_length=1e6, ax=None, max_ca_time=10):
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
        df2['tmrca'] = 1e-9
        if ca_times is not None:
            p_get_tmrca = functools.partial(
                    get_tmrca, ind2=ind, ca_times=ca_times
                    )
            df2 = df2.assign(
                    tmrca=df2['ind2'].apply(p_get_tmrca), axis=1
                    )

        df2 = df2[df2['tmrca'] <= max_ca_time]

        ## We can have TMRCA of 0 if max_ca_time is different than the value used
        ## to generate the IBD array
        # if df2['tmrca'].min() == 0:
        #     if warned is False:
        #         print("Warning - filtering TMRCA >", max_ca_time)
        #         warned = True
        #     df2 = df2[df2['tmrca'] > 0]

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

    # diploid = args.diploid
    cur_subplot = 0
    num_subplots = 2
    if args.genizon_ibd_file:
        num_subplots = 3

    ## Set up plot axes
    fig, ax_arr = plt.subplots(num_subplots, 1, figsize=(7, 7), sharex=True,
            sharey=True)

    ## Load Genizon IBD and plot
    print("Loading Genizen data")
    if args.genizon_ibd_file:
        ibd_file = os.path.expanduser(args.genizon_ibd_file)
        loaded = np.load(ibd_file)
        ibd_array = loaded['ibd_array']
        ibd_df = ibd_list_to_df(ibd_array, diploid=True, germline=True)
        plot_ibd_df(ibd_df, ax=ax_arr[0])
        ax_arr[cur_subplot].set_title('Genizon Data')
        cur_subplot += 1

    ## Get DTWF common ancestor times for IBD segments
    print("Loading DTWF")
    ca_times_dtwf = scipy.sparse.load_npz(
                os.path.expanduser(args.dtwf_ca_file)).tolil()

    ## Load DTWF IBD and plot
    ibd_file = os.path.expanduser(args.dtwf_ibd_file)
    loaded = np.load(ibd_file)
    ibd_array = loaded['ibd_array']
    ibd_df = ibd_list_to_df(ibd_array, diploid=False)
    plot_ibd_df(ibd_df, ca_times_dtwf, ax=ax_arr[cur_subplot],
            min_length=args.min_length, max_ca_time=args.max_ibd_time_gens)
    ax_arr[cur_subplot].set_title('msprime (WF)')
    cur_subplot += 1

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
    ibd_df = ibd_list_to_df(ibd_array, diploid=False)
    plot_ibd_df(ibd_df, ca_times_hudson, ax=ax_arr[cur_subplot],
            min_length=args.min_length, max_ca_time=args.max_ibd_time_gens)
    ax_arr[cur_subplot].set_title('msprime (Hudson)')
    cur_subplot += 1

    ## Plot expected IBD cluster means
    max_theory_ca_time = 5
    K = 22
    length_in_morgans = ts.get_sequence_length() / 1e8

    plot_expected_ibd(max_theory_ca_time, length_in_morgans, ax=ax_arr[0],
            K=K, diploid=True, avuncular=True, grandparental=True, min_ca_time=2,
            label_points=True)
    plot_expected_ibd(max_theory_ca_time, length_in_morgans, ax=ax_arr[1],
            K=K, diploid=False, label_points=True)
    plot_expected_ibd(max_theory_ca_time, length_in_morgans, ax=ax_arr[2],
            K=K, diploid=False, label_points=False)

    ## Legend only on upper plot
    ax_arr[0].legend(loc='upper left')
    ax_arr[1].legend(loc='upper left')
    ax_arr[2].legend(loc='upper left')

    ## Jump through a few hoops to set colourbars
    sm = plt.cm.ScalarMappable(cmap='viridis',
            norm=plt.Normalize(vmin=0, vmax=args.max_ibd_time_gens))
    sm._A = []
    cbar = fig.colorbar(sm, aspect=50, ax=ax_arr.ravel().tolist())
    cbar.ax.get_yaxis().labelpad = 5
    cbar.ax.set_ylabel('TMRCA (generations)', rotation=90)

    try:
        fig.savefig(os.path.expanduser(args.outfile), dpi=600)
        print("Plotted at high dpi")
    except:
        fig.savefig(os.path.expanduser(args.outfile))
        print("Failed to plot at high dpi - falling back to default")



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--dtwf_ibd_file", required=True)
    parser.add_argument("--dtwf_ca_file", required=True)
    parser.add_argument("--dtwf_ts_file", required=True)
    parser.add_argument("--hudson_ibd_file", required=True)
    parser.add_argument("--hudson_ca_file", required=True)
    parser.add_argument("--hudson_ts_file", required=True)
    parser.add_argument("--genizon_ibd_file")
    parser.add_argument("--outfile", required=True)
    # parser.add_argument("--diploid", action='store_true')
    parser.add_argument("--max_ibd_time_gens", type=int, default=5)
    parser.add_argument('--min_length', type=float, default=5e6)
    args = parser.parse_args()

    main(args)
