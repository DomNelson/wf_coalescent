import numpy as np
import matplotlib.pyplot as plt
import seaborn
import msprime
from collections import defaultdict


class TreeSequenceAncestry:
    def __init__(self, length, admixture_time=5, admixture_prop=0.3):
        self.ts = None
        self.length = length
        self.admixture_time = admixture_time
        self.admixture_prop = admixture_prop
        self.split_time = self.admixture_time + 100
    
    
    def simulate(self):
        """
        Simulates a recent admixture event between two diverged populations
        """
        population_configurations = [
            msprime.PopulationConfiguration(sample_size=300, initial_size=300),
            msprime.PopulationConfiguration(sample_size=0, initial_size=1)
        ]

        demographic_events = [
            msprime.MassMigration(
                time=self.admixture_time,
                source=0,
                dest=1,
                proportion=self.admixture_prop),
            msprime.PopulationParametersChange(
                time=self.admixture_time+1,
                initial_size=1,
                growth_rate=0,
                population_id=0),
            msprime.MassMigration(
                time=self.split_time,
                source=1,
                dest=0,
                proportion=1.0)
        ]

        migration_matrix = None

        self.ts = msprime.simulate(
            length=self.length,
            recombination_rate=1e-8,
            mutation_rate=2e-8,
            population_configurations=population_configurations,
            migration_matrix=migration_matrix,
            demographic_events=demographic_events
        )
        
        
    def get_local_ancestry(self, tree):
        """
        Returns the ancestry of each sample over the region
        spanned by the given tree
        """
        ## We have specified bottlenecks in each population so that
        ## they should each have a single founder
        pop_founders = tree.children(tree.root)
        pops = [tree.population(f) for f in pop_founders]
        assert len(set(pops)) == 2
        assert np.max([tree.time(x) for x in pop_founders]) < self.split_time
        
        
        local_ancestry = np.zeros(tree.num_samples()) - 1
        for founder, pop in zip(pop_founders, pops):
            sample_indices = np.array(list(tree.get_leaves(founder)), dtype=int)
            assert set(local_ancestry[sample_indices]) == set([-1])
            local_ancestry[sample_indices] = pop
            
        assert np.min(local_ancestry) == 0
        
        return local_ancestry
        
        
    def get_pop_ancestry(self):
        """
        Returns an array of sample ancestry tracts, stored as
        corresponding ancestral population and tract length
        """
        assert self.ts is not None
        
        n_samp = self.ts.num_samples
        n_trees = self.ts.num_trees
        
        ancestry = np.zeros((n_samp, n_trees)) - 1
        segment_lengths = np.zeros(n_trees)
        
        for i, tree in enumerate(self.ts.trees()):
            local_ancestry = self.get_local_ancestry(tree)
            ancestry[:, i] = local_ancestry
            segment_lengths[i] = tree.length
            
        return ancestry, segment_lengths
    

def get_ancestry_mean_var(ancestry, segment_lengths):
    ancestry_props = defaultdict(list)
    length_props = segment_lengths / np.sum(segment_lengths)

    for ind_ancestry in ancestry:
        pops = sorted(set(ind_ancestry))
        for pop in pops:
            pop_segments = np.where(ind_ancestry == pop)[0]
            pop_length = np.sum(length_props[pop_segments])
            ancestry_props[pop].append(pop_length)

    for k, v in ancestry_props.items():
        mean = np.mean(v)
        var = np.var(v)
        
    return mean, var


def get_mean_var_decay(length):
    means = []
    variances = []
    for admixture_time in range(1, 20):
        print(admixture_time)
        TSA = TreeSequenceAncestry(length=length, admixture_time=admixture_time)
        TSA.simulate()

        ancestry, lengths = TSA.get_pop_ancestry()
        mean, var = get_ancestry_mean_var(ancestry, lengths)
        means.append(mean)
        variances.append(var)

    return means, variances
