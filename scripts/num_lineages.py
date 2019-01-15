import sys, os
sys.path.append(os.path.expanduser('~/project/msprime'))
import msprime
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import bisect
from collections import Counter, defaultdict, deque


class SegmentQueue:
    def __init__(self, ts):
        self.ts = ts
        self.edge_diffs = ts.edge_diffs()
        self.trees = ts.trees()
        self.lineage_times = None
        self.segment_queue = None
        
        self.node_times = self.get_node_times()
        
        
    def count_lineages(self):
        self._init_segment_queue()
        
        times = [0]
        num_lineages = [len(self.segment_queue)]
        
        while len(self.segment_queue) > 0:
            try:
                next_add_time, _ = self.lineage_times[-1]
            except IndexError:
                next_add_time = np.inf
            next_remove_time = -self.segment_queue[-1]

            ## Doing these separately means we do both if times are equal
            if next_add_time <= next_remove_time:
                next_time = next_add_time
                while len(self.lineage_times) > 0 and next_add_time == self.lineage_times[-1][0]:
                    _, end_time = self.lineage_times.pop()
                    bisect.insort(self.segment_queue, -end_time)
            
            if next_remove_time <= next_add_time:
                next_time = next_remove_time
                while len(self.segment_queue) > 0 and next_remove_time == -self.segment_queue[-1]:
                    self.segment_queue.pop()

            
            assert(next_time > times[-1])
            times.append(next_time)
            num_lineages.append(len(self.segment_queue))
            
        return times, num_lineages        
        
        
    def _init_segment_queue(self):
        assert(self.segment_queue is None)
        self._fill_lineage_times()
        
        segment_queue = []
        while self.lineage_times[-1][0] == 0:
            _, end_time = self.lineage_times.pop()
            ## Add negative value so order is reversed natively
            segment_queue.append(-end_time)
            
        ## Sort queue by end time for popping off in order
        self.segment_queue = sorted(segment_queue)
        

    def get_node_times(self):
        node_times = {}
        for node in self.ts.nodes():
            node_times[node.id] = node.time

        return node_times
    
    
    def _init_lineage_times(self):
        tree = next(self.trees)
        diff = next(self.edge_diffs)
        
        diff_in = diff[2]
        
        lineage_times = []
        for edge in diff_in:
            end_time = self.node_times[edge.parent]
            start_time = self.node_times[edge.child]
            lineage_times.append((start_time, end_time))
            
        return lineage_times
            
            
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
            in_inds = []
            for edge_in in diff_in:
                in_inds.extend([edge_in.parent, edge_in.child])
            # in_parents = [edge.parent for edge in diff_in]
            out_candidates = [c for c in out_candidates if c not in in_inds]
        
        assert(len(out_candidates) == 1)
        out_node = out_candidates[0]
        
        return out_node, in_node
    
    
    def get_new_lineage(self, diff_out, diff_in):
        """
        Returns new lineage as tuple (start_time, end_time). Ambiguous start times
        average over possible values.
        """
        out_node, in_node = self.get_node_diff(diff_out, diff_in)
        
        ## Child of new lineage is child of removed node
        out_children = [edge.child for edge in diff_out if edge.parent == out_node]
        
        ## Child of new lineage is child of new node
        end_time = self.node_times[in_node]
        new_children = [edge.child for edge in diff_in 
                        if (edge.child in out_children and edge.parent == in_node)]
        
        if len(new_children) > 1:
            ## Ambiguity only occurs when we reattach to the same lineage
            assert(len(diff_in) <= 3)
        
        ## We don't know when the lineage split, so we take the mean time of the new lineage,
        ## or average those mean values if there are more than one.
        mean_times = [np.mean([end_time, self.node_times[c]]) for c in new_children]
        start_time = np.mean(mean_times)
#         start_time = np.min([self.node_times[c] for c in new_children])
        
        return (start_time, end_time)
            
            
    def _fill_lineage_times(self):
        assert(self.lineage_times is None)
        lineage_times = self._init_lineage_times()
        
        for tree, diff in zip(self.trees, self.edge_diffs):
            diff_out = diff[1]
            diff_in = diff[2]
            start_time, end_time = self.get_new_lineage(diff_out, diff_in)
            
            lineage_times.append((start_time, end_time))
            
        ## Sort by start time for adding to queue
        key = lambda x: x[0]
        self.lineage_times = sorted(lineage_times, key=key, reverse=True)

for model in ['hudson', 'dtwf']:
    print("Simulating", model)
    Ne = 500
    sample_size = Ne

    if model == 'dtwf':
        ## Wright Fisher model has diploid individuals
        Ne = int(Ne / 2)

    ts = msprime.simulate(
            sample_size=sample_size,
            Ne=Ne,
            length=1e8,
            recombination_rate=1e-8,
            model=model
            )
    S = SegmentQueue(ts)
    times, num_lineages = S.count_lineages()
    plt.plot(times, num_lineages)

plt.ylabel("Number of Lineages")
plt.xlabel("Generation")
# plt.xlim((0, 100))
plt.savefig('../results/test_plot.png')


