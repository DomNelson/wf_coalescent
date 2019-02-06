import msprime
import numpy as np
import scipy.sparse
from collections import defaultdict
import bisect
from itertools import combinations
import dask
import sparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
from tqdm import tqdm
import sys, os
from sortedcontainers import sortedlist
import gc
import weakref


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


class SortedMap:
    def __init__(self):
        self.keys = sortedlist.SortedList()
        self.vals = sortedlist.SortedList()
        
    def __setitem__(self, key, val):
        idx = self.keys.bisect_left(key)
        if idx == len(self.keys):
            self.keys.add(key)
            self.vals.add(val)
            return
        
        if self.keys[idx] != key:
            self.keys.add(idx, key)
                
        self.vals.add(idx, val)
    
    def __getitem__(self, key, pop=False):
        idx = bisect.bisect_left(self.keys, key)
        try:
            if self.keys[idx] != key:
                raise KeyError
        except IndexError:
            raise KeyError
            
        val = self.vals[idx]
        
        if pop:
            self.keys.pop(idx)
            self.vals.pop(idx)
            
        return val
    
    def items(self):
        for i in range(len(self.keys)):
            yield self.keys[i], self.vals[i]
    
    def __contains__(self, key):
        idx = bisect.bisect_left(self.keys, key)
        return len(self.keys) > 0 and idx < len(self.keys) and self.keys[idx] == key
    
    def pop(self, key):
        return self.__getitem__(key, pop=True)


class EdgeList:
    def __init__(self, ts, max_time):
        self.ts = ts
        self.max_time = max_time
        self.bps = list(ts.breakpoints())
        self.edges = list(ts.edges())
        self.regions = [{} for _ in range(len(self.bps) - 1)]
        self.node_times = self.get_node_times()
        self.ibd_list = []

        ## This is just a container for the clusters pointed to in self.regions
        self.clusters = [[] for _ in range(len(self.regions))]
        
        ## Flag to start using gc.collect()
        self.collect = False
            
        
    def get_regions(self, edge):
        ## Finds the indices of ts breakpoints which span the region
        ## of the edge
        first = bisect.bisect_left(self.bps, edge.left)
        last = first
        while edge.right > self.bps[last]:
            last += 1
            
        return first, last
    
    
    def check_edges(self, edges):
        ## Edges should all have the same parent and region
        parents = [e.parent for e in edges]
        assert len(set(parents)) == 1

    
    def add_edges(self, edges):
        remove_parent_from_regions = set()
        self.check_edges(edges)
        
        for edge in edges:
            first, last = self.get_regions(edge)
                
            parent = edge.parent
            if self.node_times[parent] > self.max_time:
                return
            
            ## Start using gc.collect() once we're below self.max_time
            if self.collect is False:
                gc.collect()
            self.collect = True

            for i in range(first, last):
                remove_parent_from_regions.add(i)
                
                ## If this edge is a descendant of an existing IBD cluster,
                ## we add the child to the cluster. Otherwise we start a
                ## new cluster.
                if parent in self.regions[i]:
                    parent_cluster = self.regions[i][parent]
                else:
                    parent_cluster = []
                    self.regions[i][parent] = parent_cluster
                    self.clusters[i].append(parent_cluster)
                    
                self.regions[i][edge.child] = parent_cluster

        ## Once we have added all children to the clusters, we remove the parents
        for i in remove_parent_from_regions:
            try:
                self.regions[i].pop(parent)
            except KeyError:
                pass

            
    def get_node_times(self):
        node_times = {}
        for node in self.ts.nodes():
            node_times[node.id] = node.time

        return node_times
    
    
    def get_next_edges(self):
        ## Pulls out next group of edges which share the same parent. Groups
        ## them in a dict by region.
        next_parent = self.edges[-1].parent
        next_edges = []
        while len(self.edges) > 0 and self.edges[-1].parent == next_parent:
            edge = self.edges.pop()
            next_edges.append(edge)
                        
        return next_edges
    
    
    def _build_clusters(self):
        batch_size = 1000
        current_batch = 0
        with tqdm(total=len(self.edges), desc="Reading edges") as pbar:
            while len(self.edges) > 0:
                start_edges = len(self.edges)
                edges_dict = self.get_next_edges()
                self.add_edges(edges_dict)
                
                current_batch += start_edges - len(self.edges)
                if current_batch > batch_size and self.collect is True:
                    gc.collect()
                    pbar.update(current_batch)
                    current_batch = 0
        
        with tqdm(total=len(self.regions), desc="Building clusters") as pbar:
            for r in self.regions:
                for sample, cluster in r.items():
                    cluster.append(sample)
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


    def write_ibd_to_file(self, out_file):
            out_array = np.array(self.ibd_list)
            np.savez_compressed(out_file, ibd_array=out_array)
        
                
    def get_ibd(self):
        self._build_clusters()
        
        n = ts.num_samples
        ibd_start = scipy.sparse.lil_matrix((n, n), dtype=int)
        ibd_end = scipy.sparse.lil_matrix((n, n), dtype=int)
        
        ## NOTE 1: We compute the upper-triangular portion of the IBD
        ## matrix only. If a < b, IBD[b, a] == 0 does NOT imply no IBD
        ## shared, only IBD[a, b] is meaningful.
        ##
        ## NOTE 2: Indices are offset +1 so that sparse matrix default fill
        ## of zero implies no IBD, not IBD starting at 0 in genomic coords.
        ## IBD starting at 0 is denoted by index 1.
        with tqdm(total=len(self.clusters), desc="Writing IBD pairs") as pbar:
            for i in range(len(self.clusters)):
                for cluster in self.clusters[i]:
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
        
        
def plot_ibd(ibd_list, min_length=0, out_file=None):
    if out_file is None:
        out_file = os.path.expanduser('~/temp/ibd_plot.png')
	
    cols = ["ind1", "ind2", "start", "end"]

    ibd_df = pd.DataFrame(ibd_list, columns=cols)
    ibd_df['len'] = ibd_df['end'] - ibd_df['start']
    ibd_df = ibd_df[ibd_df['len'] >= min_length]
    
    ibd_df['count'] = 1
    for ind in set(ibd_df['ind1'].values):
        ind_df = ibd_df[ibd_df['ind1'] == ind]
        ind_pairwise_ibd_df = ind_df.groupby('ind2').sum()

        total_IBD = ind_pairwise_ibd_df['len'].values
        num_segments = ind_pairwise_ibd_df['count'].values

        plt.plot(total_IBD, num_segments, 'o')

    plt.ylabel('Number of IBD segments')
    plt.xlabel('Total IBD')
    plt.xscale('log')
    plt.savefig(out_file)


if __name__ == "__main__":
    rho = 1e-8
    chrom_lengths = [247249719, 242951149, 199501827, 191273063, 180857866,
            170899992, 158821424, 146274826, 140273252, 135374737, 134452384,
            132349534, 114142980, 106368585, 100338915, 88827254, 78774742,
            76117153, 63811651, 62435964, 46944323, 49691432]
    num_loci = chrom_lengths[-1] + 1
    # chrom_lengths = chrom_lengths[:3]

    positions, rates = get_positions_rates(chrom_lengths, rho)
    recombination_map = msprime.RecombinationMap(
            positions, rates, num_loci=num_loci
    )
    population_configuration = msprime.PopulationConfiguration(
	sample_size=100,
	initial_size=500
    )

    model = 'dtwf'
    outfile = '/gs/scratch/dnelson/project/wf_coalescent/Ne500_samples100_WG_' + model # .npz extension is added by numpy
    ts = msprime.simulate(population_configurations=[population_configuration], model=model, recombination_map=recombination_map)

    e = EdgeList(ts, 10)
    e.get_ibd()
    e.write_ibd_to_file(out_file=outfile)
