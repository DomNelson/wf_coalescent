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
            
            
    def get_node_diff(self, diff_out, diff_in):
        out_candidates = self.get_node_diff_candidates(diff_out)        
        in_candidates = self.get_node_diff_candidates(diff_in)
        
        if len(in_candidates) > 1:
            out_inds = []
            for edge_out in diff_out:
                out_inds.extend([edge_out.parent, edge_out.child])
            in_candidates = [c for c in in_candidates if c not in out_inds]
        
        assert(len(in_candidates) == 1)
        in_node = in_candidates[0]
        
        if len(out_candidates) > 1:
            in_parents = [edge.parent for edge in diff_in]
            out_candidates = [c for c in out_candidates if c not in in_parents]
        
        assert(len(out_candidates) == 1)
        out_node = out_candidates[0]
        
        return out_node, in_node
    
    
#     def get_new_lineage(self, diff_out, diff_in):
#         """
#         Returns new lineage as tuple (start_time, end_time). Ambiguous start times
#         average over possible values.
#         """
#         out_node, in_node = self.get_node_diff(diff_out, diff_in)
            
    def get_first_ancestor_below_max_time(self, tree, node):
        if self.node_times[node] > self.max_time:
            print(-1)
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
            
        print(tree.draw(format='unicode'))
        print("Updating:", samples)
            
        return list(samples)
    
    
    def get_tree_min_common_ancestor_times(self, tree, edge_diff):
        nodes_in_diff = self.get_samples_below_edge_diff(tree, edge_diff)
        
        for a, b in combinations(nodes_in_diff, 2):
            if a == b:
                continue
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
            if locus is not None and (edge.left <= locus and edge.right > locus):
                continue
            parent_time = self.node_times[edge.parent]
            child_time = self.node_times[edge.child]
            if parent_time > self.max_time and child_time <= self.max_time:
                children.append(edge.child)
                
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

        print("Initial clusters:", clusters)
        
        return clusters
            
            
    def update_clusters(self, tree, diff, clusters):
        new_clusters = defaultdict(set)
        changed_samples = self.get_samples_below_edge_diff(tree, diff)
#         segment, diff_in, diff_out = diff
#         out_node, in_node = self.get_node_diff(diff_out, diff_in)
        
        ## First remove changed samples from existing clusters and update
        ## oldest anc if necessary
        for anc, cluster in clusters.items():
            new_anc = self.get_first_ancestor_below_max_time(tree, anc)
            new_clusters[new_anc].update(cluster.difference(changed_samples))
            
        ## Now add to new clusters
        for sample in changed_samples:
            anc = self.get_first_ancestor_below_max_time(tree, sample)
            new_clusters[anc].add(sample)
            
        ## Sanity check - no clusters should share samples
        for c1, c2 in combinations(new_clusters.values(), 2):
            assert len(c1.intersection(c2)) == 0

        print(clusters)
        print(new_clusters)
            
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
            print(i)
            for cluster in clusters.values():
                print(cluster)
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
    recmap = msprime.RecombinationMap(
	positions=[0, 10],
	rates=[0.01, 0],
	num_loci = 10
    )
    population_configuration = msprime.PopulationConfiguration(
	sample_size=5,
	initial_size=10
    )

    ts = msprime.simulate(
            population_configurations=[population_configuration],
            model='dtwf',
            recombination_map=recmap
            )

    # ts_file = '/Users/dnelson/temp/ts_dump.h5'
    # ts.dump(ts_file)
    # ts = msprime.load(ts_file)

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


