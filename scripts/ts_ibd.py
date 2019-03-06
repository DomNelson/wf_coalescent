import sys, os
sys.path.append('/home/dnelson/project/msprime/')
sys.path.append('/home/dnelson/project/msprime/lib/subprojects/git-submodules/tskit/python/')
import msprime
import numpy as np
import scipy.sparse
from collections import defaultdict
import bisect
from itertools import combinations
import dask
import sparse
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
from tqdm import tqdm
from sortedcontainers import sortedlist
import gc
import weakref
# from profilehooks import profile


def expected_num_segs(t, K, L):
    return (K + 2 * L * t) / 2 ** (2 * t - 1)


def expected_total_length(t, L):
    return L / 2 ** (2 * t - 1)


class TSRelatives:
    def __init__(self, max_time, ts=None, ts_file=None):
        if ts is None and ts_file is None:
            print("One of ts or ts_file must be specified")
            raise ValueError

        if ts is None:
            ts = msprime.load(ts_file)

        self.ts = ts
        self.max_time = max_time
        self.bps = list(self.ts.breakpoints())

        self.node_times = self.get_node_times()
        self.ca_times = scipy.sparse.lil_matrix((ts.num_samples, ts.num_samples))
        self.ca_last = scipy.sparse.lil_matrix((ts.num_samples, ts.num_samples))
        self.ca_count = scipy.sparse.lil_matrix((ts.num_samples, ts.num_samples))
        self.ibd_list = []


    def list_edges(self, max_items=100):
        edges = self.ts.edges()
        edges_list = [n for i, n in enumerate(edges) if i < max_items]

        return edges_list


    def list_nodes(self, max_items=100):
        nodes = self.ts.nodes()
        nodes_list = [n for i, n in enumerate(nodes) if i < max_items]

        return nodes_list


    def draw_trees(self, first=0, last=5, locus=None, times=False):
        node_labels = None
        if locus is not None:
            print("Searching for tree spanning position", locus, "- " +\
                    "ignoring tree indices")
            first = 0
            last = self.ts.num_trees

        trees = self.ts.trees()
        for i, tree in enumerate(trees):
            start, end = tree.get_interval()
            if locus is not None:
                if start > locus or end <= locus:
                    continue

            if i >= first and i < last:
                if times is True:
                    make_label = lambda n: str(n) + ': ' + str(self.node_times[n])
                    node_labels = {n: make_label(n) for n in tree.nodes()}

                ## TODO: This logic is pretty ugly
                print("Interval:", start, "-", end)
                print(tree.draw(format='unicode', node_labels=node_labels))
                if locus is not None:
                    break
            
            
    def get_node_times(self):
        node_times = {}
        for node in self.ts.nodes():
            node_times[node.id] = node.time
            
        node_times[-1] = np.inf

        return node_times
    
    
    def get_first_ancestor_below_max_time(self, tree, node):
        if self.node_times[node] > self.max_time:
            return -1
        
        parent = tree.get_parent(node)
        parent_time = self.node_times[parent]
        while parent >= 0:
            if parent_time > self.max_time:
                break
            node = parent
            parent = tree.get_parent(node)
            parent_time = self.node_times[parent]

        # print(node)
        assert self.node_times[node] <= self.max_time
        
        return node
    
    
    def get_samples_below_edge_diff_conservative(self, tree, edge_diff):
        nodes = set()
        samples = set()
        segment, edges_in, edges_out = edge_diff
        assert type(edges_in) == list
        
        for edge in edges_in + edges_out:
            oldest_relevant_anc = self.get_first_ancestor_below_max_time(
                    tree, edge.child)
            if oldest_relevant_anc >= 0:
                nodes.add(oldest_relevant_anc)
            
        oldest_change = 0
        for node in nodes:
            node_time = self.node_times[node]
            if node_time > oldest_change:
                oldest_change = node_time
            samples.update(tree.get_leaves(node))
            
        return list(samples), oldest_change
    
    
    def get_samples_below_edge_diff(self, tree, edge_diff):
        samples = set()
        segment, edges_in, edges_out = edge_diff
        assert type(edges_in) == list
        
        oldest_change = 0
        for edge in edges_in + edges_out:
            node_time = self.node_times[edge.parent]
            if node_time > oldest_change:
                oldest_change = node_time
            s = tree.get_leaves(edge.parent)
            samples.update(s)
            
        return list(samples), oldest_change
    
    
    def get_tree_min_common_ancestor_times(self, tree, edge_diff):
        nodes_in_diff, _ = self.get_samples_below_edge_diff_conservative(tree, edge_diff)
        
        for a, b in combinations(nodes_in_diff, 2):
            if a == b:
                continue
            a, b = sorted([a, b])
            ca = tree.get_mrca(a, b)
            ca_time = self.node_times[ca]
            if ca_time > self.max_time:
                continue
            stored_time = self.ca_times[a, b]
            if stored_time == 0 or ca_time < stored_time:
                ## If stored_time == 0 we're recording the first common ancestor
                self.ca_times[a, b] = ca_time
                self.ca_last[a, b] = ca
                self.ca_count[a, b] = 1
                
            elif ca_time == stored_time and self.ca_last[a, b] != ca:
                ## Here we've found a second common ancestor from the same
                ## generation as the first
                self.ca_count[a, b] = 2
                
                
    def get_all_min_common_ancestor_times(self):
        with tqdm(total=self.ts.num_trees) as pbar:
            for tree, edge_diff in zip(self.ts.trees(), self.ts.edge_diffs()):
                self.get_tree_min_common_ancestor_times(tree, edge_diff)
                pbar.update(1)

                
    def write_ibd_pair(self, ibd_start, ibd_end, pair):
        if ibd_start[pair] != 0:
            assert ibd_end[pair] != 0
            ind1, ind2 = pair
            
            ## Here we undo the offset we used when building segments
            start_idx = ibd_start[pair] - 1
            end_idx = ibd_end[pair] - 1
            
            start = self.bps[start_idx]
            end = self.bps[end_idx]
            record = [ind1, ind2, start, end]
            self.ibd_list.append(record)
            
            ## Reset the pair for building the next segment
            ibd_start[pair] = 0
            ibd_end[pair] = 0
            
            
    def write_ibd_all(self, ibd_start, ibd_end):
        coo = ibd_start.tocoo()
        for ind1, ind2, start in zip(coo.row, coo.col, coo.data):
            pair = (ind1, ind2)
            self.write_ibd_pair(ibd_start, ibd_end, pair)
            
            
    def get_children_of_edges_crossing_max_time(self, locus=None):
        children = []
        for edge in self.ts.edges():
            if locus is not None and (edge.left > locus or edge.right <= locus):
                continue
            parent_time = self.node_times[edge.parent]
            child_time = self.node_times[edge.child]
            if parent_time > self.max_time and child_time <= self.max_time:
                children.append(edge.child)
                
        return children
    

    def get_initial_clusters(self, first_tree):
        # self.draw_trees(locus=0)
        ancs = self.get_children_of_edges_crossing_max_time(locus=0)
        clusters = defaultdict(set)
        
        ## If no clusters cross max_time, all nodes must be below max_time,
        ## so everyone is IBD
        if len(ancs) == 0:
            clusters[first_tree.get_root()] = set(first_tree.samples())
            return clusters
        
        for a in ancs:
            clusters[a] = set(first_tree.get_leaves(a))

        # print("Initial clusters:")
        # print(clusters)
        
        return clusters
            
            
    def update_clusters(self, tree, diff, clusters):
        changed_samples, oldest_time = self.get_samples_below_edge_diff(tree, diff)
        
        ## First remove changed samples from existing clusters and update
        ## oldest anc if necessary
        ancs = list(clusters.keys())
        num_updates = 0
        for anc in ancs:

            if self.node_times[anc] <= oldest_time:
                ## These ancs may no longer be the closest to max_time
                new_anc = self.get_first_ancestor_below_max_time(tree, anc)
                cluster = clusters.pop(anc)
                clusters[new_anc].update(cluster.difference(changed_samples))
                num_updates += 1
            
        ## Now add to new clusters
        for sample in changed_samples:
            anc = self.get_first_ancestor_below_max_time(tree, sample)
            clusters[anc].add(sample)
            
        ## Sanity check - no clusters should share samples
        # for c1, c2 in combinations(new_clusters.values(), 2):
        #     assert len(c1.intersection(c2)) == 0

        # print(len(ancs), num_updates)

        return clusters
        
                
    # @profile
    def get_ibd(self):        
        n = self.ts.num_samples
        ibd_start = scipy.sparse.lil_matrix((n, n), dtype=int)
        ibd_end = scipy.sparse.lil_matrix((n, n), dtype=int)
        
        ## Initialize with the first tree
        trees = self.ts.trees()
        diffs = self.ts.edge_diffs()
        first_tree = next(trees)
        first_diff = next(diffs)
        clusters = self.get_initial_clusters(first_tree)
        
        ## NOTE 1: We compute the upper-triangular portion of the IBD
        ## matrix only. If a < b, IBD[b, a] == 0 does NOT imply no IBD
        ## shared, only IBD[a, b] is meaningful.
        ##
        ## NOTE 2: Indices are offset +1 so that sparse matrix default fill
        ## of zero implies no IBD, not IBD starting at 0 in genomic coords.
        ## IBD starting at 0 is denoted by index 1.
        with tqdm(total=self.ts.num_trees, desc="Writing IBD pairs") as pbar:
            for i, (tree, diff) in enumerate(zip(trees, diffs)):
#                 print(i)
                for cluster in clusters.values():
#                     print(cluster)
                    ibd_pairs = combinations(sorted(cluster), 2)
                    for pair in ibd_pairs:
                        ## Check if we are starting a new IBD segment
                        ## or continuing one
                        if ibd_end[pair] == i + 1:
                            ibd_end[pair] += 1
                        else:
                            ## If we start a new segment, write the old one first
                            self.write_ibd_pair(ibd_start, ibd_end, pair)
                            ibd_start[pair] = i + 1
                            ibd_end[pair] = i + 2
                clusters = self.update_clusters(tree, diff, clusters)
                pbar.update(1)
        
        ## TODO: Add last tree - integrate into main loop above
            i += 1
            # print(i)
            for cluster in clusters.values():
                # print(cluster)
                ibd_pairs = combinations(sorted(cluster), 2)
                for pair in ibd_pairs:
                    ## Check if we are starting a new IBD segment
                    ## or continuing one
                    if ibd_end[pair] == i + 1:
                        ibd_end[pair] += 1
                    else:
                        ## If we start a new segment, write the old one first
                        self.write_ibd_pair(ibd_start, ibd_end, pair)
                        ibd_start[pair] = i + 1
                        ibd_end[pair] = i + 2
            pbar.update(1)

        ## Write out all remaining segments, which reached the end of the
        ## simulated region
        self.write_ibd_all(ibd_start, ibd_end)
        

def get_tree_spanning_locus(locus):
    for t in ts.trees():
        if t.get_interval()[0] > locus:
            return tree

        
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


def plot_ibd_df(df, ca_times=None, min_length=1e6, ax=None, max_ca_time=5):
    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 4))

    df = df[df['len'] >= min_length]

    df = df.assign(count=np.ones(df.shape[0]))
    max_ca_time = 0
    for ind in set(df['ind1'].values):
        ind_df = df[df['ind1'] == ind][['ind2', 'len', 'count']]
        df2 = ind_df.groupby('ind2').sum()
        df2 = df2.reset_index()
        
        df2 = df2.assign(colours=np.zeros(df2.shape[0]))
        if ca_times is not None:
            get_ca_time = lambda x: ca_times[sorted([ind, x])[0], sorted([ind, x])[1]]
            df2 = df2.assign(
                    tmrca=df2['ind2'].apply(get_ca_time), axis=1
                    )

            if df2['tmrca'].max() > max_ca_time:
                max_ca_time = df2['tmrca'].max()

        ## TODO: Should assert that all ca_times > 0
        df2 = df2[df2['tmrca'] <= max_ca_time]
        num_segments = df2['count'].values
        total_IBD = df2['len'].values
        colours = df2['tmrca'].values
    
        ax.scatter(total_IBD, num_segments, c=colours, s=5,
                cmap='viridis', vmin=0, vmax=max_ca_time)

    ax.set_ylabel('Number of IBD segments')
    ax.set_xlabel('')
    if ax.is_last_row():
        ax.set_xlabel('Total IBD')
    ax.set_xscale('log')
    ax.set_xlim((0.7 * min_length, 4e9))

    return ax


def paper_plot():
    ts_file = '/home/dnelson/project/wf_coalescent/results/' +\
            'IBD/Ne500_samples500_WG_ts_dtwf.h5'
    ts = msprime.load(ts_file)
    max_ibd_time_gens = 10

    ## Set up plot axes
    fig, ax_arr = plt.subplots(2, 1, figsize=(7, 7), sharex=True, sharey=True)

    ## Get DTWF common ancestor times for IBD segments
    print("Loading DTWF")
    T_dtwf = TSRelatives(max_ibd_time_gens, ts)
    T_dtwf.get_all_min_common_ancestor_times()
    # T_dtwf.get_ibd()

    ## Load DTWF IBD and plot
    ibd_file = '/home/dnelson/project/wf_coalescent/results/' +\
            'IBD/Ne500_samples500_WG_dtwf.npz'
    loaded = np.load(ibd_file)
    ibd_array = loaded['ibd_array']
    ibd_df = ibd_list_to_df(ibd_array)
    plot_ibd_df(ibd_df, T_dtwf.ca_times, ax=ax_arr[0])
    ax_arr[0].set_title('DTWF')

    ## Get Hudson common ancestor times for IBD segments
    print("Loading Hudson")
    ts_file = '/home/dnelson/project/wf_coalescent/results/' +\
            'IBD/Ne500_samples500_WG_ts_hudson.h5'
    ts = msprime.load(ts_file)
    T_hudson = TSRelatives(max_ibd_time_gens, ts)
    T_hudson.get_all_min_common_ancestor_times()

    ## Load Hudson IBD and plot
    ibd_file = '/home/dnelson/project/wf_coalescent/results/' +\
            'IBD/Ne500_samples500_WG_hudson.npz'
    loaded = np.load(ibd_file)
    ibd_array = loaded['ibd_array']
    ibd_df = ibd_list_to_df(ibd_array)
    plot_ibd_df(ibd_df, T_hudson.ca_times, ax=ax_arr[1])
    ax_arr[1].set_title('Hudson')

    ## Plot expected IBD cluster means
    max_theory_ca_time = 5
    K = 22
    length_in_morgans = ts.get_sequence_length() / 1e8
    plot_expected_ibd(max_theory_ca_time, length_in_morgans, ax=ax_arr[0], K=K)
    plot_expected_ibd(max_theory_ca_time, length_in_morgans, ax=ax_arr[1], K=K)

    ## Legend only on upper plot
    ax_arr[0].legend()

    ## Jump through a few hoops to set colourbars
    sm = plt.cm.ScalarMappable(cmap='viridis',
            norm=plt.Normalize(vmin=0, vmax=max_ibd_time_gens))
    sm._A = []

    # for ax in ax_arr:
    # cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
    # cbar = fig.colorbar(sm, cax=cbar_ax)
    cbar = fig.colorbar(sm, aspect=50, ax=ax_arr.ravel().tolist())
    cbar.ax.get_yaxis().labelpad = 5
    cbar.ax.set_ylabel('TMRCA', rotation=90)

    return T_dtwf, T_hudson, fig, ax_arr


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python ts_ibd.py paper_plot_outfile")
        sys.exit()

    outfile = os.path.expanduser(sys.argv[1])

    T_dtwf, T_hudson, fig, ax_arr = paper_plot()

    fig.savefig(outfile)

    import IPython; IPython.embed()
    
