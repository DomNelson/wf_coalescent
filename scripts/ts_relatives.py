import sys, os
sys.path.append(
        os.path.expanduser('/home/dnelson/project/msprime/')
        )
sys.path.append(
        os.path.expanduser('/home/dnelson/project/msprime/' +\
                'lib/subprojects/git-submodules/tskit/python/')
        )
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


def get_positions_rates(chrom_lengths, rho):
    """
    Takes a list of chromosome lengths and returns lists of positions and
    rates to pass to msprime.RecombinationMap
    """
    positions = []
    rates = []
    total_length = 0
    for length in chrom_lengths:
        positions.extend([int(total_length), int(total_length) + int(length) - 1])
        rates.extend([rho, 0.5])
        total_length += length

    rates[-1] = 0

    return positions, rates


class TSRelatives:
    def __init__(self, max_time, ts=None, ts_file=None):
        if ts is None and ts_file is None:
            print("One of ts or ts_file must be specified")
            raise ValueError
            
        self.max_time = max_time
        self.bps = list(ts.breakpoints())

        
        self.ts = ts
        if self.ts is None:
            self.ts = msprime.load(ts_file)
            
        self.node_times = self.get_node_times()
        self.ca_times = scipy.sparse.lil_matrix((ts.num_samples, ts.num_samples))
        self.ca_last = scipy.sparse.lil_matrix((ts.num_samples, ts.num_samples))
        self.ca_count = scipy.sparse.lil_matrix((ts.num_samples, ts.num_samples))
        self.ibd_list = []
            
            
    def get_node_times(self):
        node_times = {}
        for node in self.ts.nodes():
            node_times[node.id] = node.time
            
        node_times[-1] = np.inf

        return node_times
    
    
    def get_node_diff_candidates(self, diff):
        """
        Returns nodes which may have been added or removed
        in the diff - will be involved in the most edge diffs
        """
        counts = defaultdict(int)
        for edge in diff:
            counts[edge.parent] += 1
            counts[edge.child] += 1
            
        max_count = max(counts.values())            
        out_candidates = [node for node, count in counts.items() if count == max_count]

        return out_candidates
            
            
    def get_first_ancestor_below_max_time(self, tree, node):
        if self.node_times[node] > self.max_time:
            return -1
        
        parent = tree.get_parent(node)
        while parent >= 0:
            parent_time = self.node_times[parent]
            if parent_time > self.max_time:
                break
            node = parent
            parent = tree.get_parent(node)

        assert self.node_times[node] <= self.max_time
        return node


    def get_samples_below_edge_diff(self, tree, edge_diff):
        nodes = set()
        samples = set()
        segment, edges_in, edges_out = edge_diff
        assert type(edges_in) == list

        for edge in edges_in + edges_out:
            oldest_relevant_anc = self.get_first_ancestor_below_max_time(
                    tree,
                    edge.child
                    )
            if oldest_relevant_anc >= 0:
                nodes.add(oldest_relevant_anc)

        for node in nodes:
            samples.update(tree.get_leaves(node))
            
        return list(samples)
    
    
    def get_tree_min_common_ancestor_times(self, tree, edge_diff):
        nodes_in_diff = self.get_samples_below_edge_diff(tree, edge_diff)
        
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
        children = set()
        for edge in self.ts.edges():
            if locus is not None and (edge.left > locus or edge.right <= locus):
                continue
            parent_time = self.node_times[edge.parent]
            child_time = self.node_times[edge.child]
            if parent_time > self.max_time and child_time <= self.max_time:
                # print("Adding", edge)
                children.add(edge.child)
                
        return children
    

    def get_initial_clusters(self, first_tree):
        ancs = self.get_children_of_edges_crossing_max_time(locus=0)
        clusters = defaultdict(set)
        
        ## If no clusters cross max_time, all nodes must be below max_time,
        ## so everyone is IBD
        if len(ancs) == 0:
            clusters[first_tree.get_root()] = set(first_tree.samples())
            return clusters
        
        for a in ancs:
            clusters[a] = set(first_tree.get_leaves(a))

        # print("Initial clusters:", clusters)

        return clusters


    # @profile
    def update_clusters(self, tree, diff, clusters):
        new_clusters = defaultdict(set)
        changed_samples = self.get_samples_below_edge_diff(tree, diff)
        # print("Changed samples:", len(changed_samples))
        # print("Num clusters:", len(clusters))
        
        ## First remove changed samples from existing clusters and update
        ## oldest anc if necessary
        for anc, cluster in clusters.items():
            if len(cluster) == 0:
                continue
            new_anc = self.get_first_ancestor_below_max_time(tree, anc)
            new_clusters[new_anc].update(cluster.difference(changed_samples))
            
        ## Now add to new clusters
        for sample in changed_samples:
            anc = self.get_first_ancestor_below_max_time(tree, sample)
            new_clusters[anc].add(sample)
            
        ## Sanity check - no clusters should share samples
        # for c1, c2 in combinations(new_clusters.values(), 2):
        #     assert len(c1.intersection(c2)) == 0

        return new_clusters
        
                
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
                # if i == 1000:
                #     print("Num clusters:", len(clusters))
                #     sys.exit()
                for cluster in clusters.values():
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
            for cluster in clusters.values():
                ibd_pairs = combinations(sorted(cluster), 2)
                for pair in ibd_pairs:
                    ## Check if we are starting a new IBD segment
                    ## or continuing one
                    if ibd_end[pair] == i + 1:
                        ibd_end[pair] += 1
                    else:
                        ## If we start a new segment, write the
                        ## old one first
                        self.write_ibd_pair(ibd_start, ibd_end, pair)
                        ibd_start[pair] = i + 1
                        ibd_end[pair] = i + 2
            pbar.update(1)

        ## Write out all remaining segments, which reached the end of the
        ## simulated region
        self.write_ibd_all(ibd_start, ibd_end)


def plot_ibd(ibd_list, ca_times=None, min_length=0, out_file=None):
    print("Plotting!")
    fig, ax = plt.subplots()
    if out_file is None:
        out_file = os.path.expanduser('~/temp/ibd_plot.png')
	
    cols = ["ind1", "ind2", "start", "end"]

    ibd_df = pd.DataFrame(ibd_list, columns=cols)
    ibd_df['len'] = ibd_df['end'] - ibd_df['start']
    ibd_df = ibd_df[ibd_df['len'] >= min_length]
    
    ibd_df['count'] = 1
    inds = set(ibd_df['ind1'].values)
    for ind in inds:
        # print(ind)
        # print(len(inds))
        ind_df = ibd_df[ibd_df['ind1'] == ind]
        ind_pairwise_ibd_df = ind_df.groupby('ind2').sum()
        ind_pairwise_ibd_df = ind_pairwise_ibd_df.reset_index()

        total_IBD = ind_pairwise_ibd_df['len'].values
        num_segments = ind_pairwise_ibd_df['count'].values
        ind2 = ind_pairwise_ibd_df['ind2'].values
        sizes = np.ones(total_IBD.shape[0]) * 8
        # print(len(sizes))
        # print(total_IBD.shape)

        colours = None
        cmap = None
        if ca_times is not None:
            colours = []
            for i in range(len(ind2)):
                x, y = sorted([ind, ind2[i]])
                assert(x <= y)
                colours.append(ca_times[x, y])
                if colours[-1] == 0:
                    print(x, y)

            colours = np.array(colours)
            cmap = 'viridis'

        # print("Plotting")
        ax.scatter(total_IBD, num_segments, s=sizes, c=colours, cmap=cmap)
        # print("Break!")
        # break

    ax.set_ylabel('Number of IBD segments')
    ax.set_xlabel('Total IBD')
    ax.set_xscale('log')
    print("Saving...")
    fig.savefig(out_file)
    print("Done!")


def ibd_list_to_df(ibd_list, ca_times=None):
    cols = ["ind1", "ind2", "start", "end"]

    df = pd.DataFrame(np.array(tsr.ibd_list), columns=cols)
    df['len'] = df['end'] - df['start']

    
    if ca_times is not None:
        df['tmrca'] = df[['ind1', 'ind2']].apply(
                lambda x: tsr.ca_times[x['ind1'], x['ind2']], axis=1
                )
    
    return df
    

def simulate():
    rho = 1e-8
    chrom_lengths = [247249719, 242951149, 199501827, 191273063, 180857866,
            170899992, 158821424, 146274826, 140273252, 135374737, 134452384,
            132349534, 114142980, 106368585, 100338915, 88827254, 78774742,
            76117153, 63811651, 62435964, 46944323, 49691432]
    num_loci = chrom_lengths[-1] + 1

    positions, rates = get_positions_rates(chrom_lengths, rho)
    recombination_map = msprime.RecombinationMap(
            positions, rates, num_loci=num_loci
    )
    population_configuration = msprime.PopulationConfiguration(
	sample_size=50,
	initial_size=50
    )

    ts = msprime.simulate(
            population_configurations=[population_configuration],
            model='dtwf',
            recombination_map=recombination_map
            )

    return ts


def load_triple_merger():
    fname = 'test/test_data/small_triple_merger.h5'
    print("A good max_time is 16")

    return msprime.load(fname)


def draw_trees(ts):
    for t in ts.trees():
        print(t.draw(format='unicode'))


def run_tsr(ts, max_time=10):
    tsr = TSRelatives(max_time, ts)
    tsr.get_all_min_common_ancestor_times()
    tsr.get_ibd()

    return tsr

# ts_file = '/home/dnelson/project/wf_coalescent/results/IBD/Ne500_samples500_WG_ts_dtwf.h5'
# ibd_file = '/home/dnelson/project/wf_coalescent/results/IBD/Ne500_samples500_WG_dtwf.npz'
# loaded = np.load(ibd_file)
# ibd_array = loaded['ibd_array']
# ts = msprime.load(ts_file)

if __name__ == "__main__":
    ts = simulate()

    tsr = TSRelatives(5, ts)
    tsr.get_ibd()
    tsr.get_all_min_common_ancestor_times()

    try:
        plot_ibd(tsr.ibd_list, tsr.ca_times, min_length=1e6)
    except:
        print("Error!")
        import IPython; IPython.embed()
    import IPython; IPython.embed()


