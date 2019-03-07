import sys, os
sys.path.append('/home/dnelson/project/msprime/')
sys.path.append('/home/dnelson/project/msprime/lib/subprojects/git-submodules/tskit/python/')
import msprime
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
from collections import Counter
from tqdm import tqdm


def get_positions_rates(num_chroms=22, rho=1e-8):
    """
    Takes a list of chromosome lengths and returns lists of positions and
    rates to pass to msprime.RecombinationMap
    """
    ## TODO: Check these lengths!
    chrom_lengths: typing.List[float] = [
            277693825, 263349606, 224483368, 212778391,
            203765845, 192951739, 186795932, 170176519,
            168073935, 178947388, 159485425, 172777271,
            126940475, 116331251, 12554709, 13489110,
            12929210, 11897848, 10779606, 10792434,
            61526812, 72706815]

    positions = []
    rates = []
    total_length = 0
    for length in chrom_lengths[:num_chroms]:
        positions.extend([int(total_length), int(total_length) + int(length) - 1])
        rates.extend([rho, 0.5])
        total_length += length

    rates[-1] = 0

    return positions, rates


def get_recombination_map(num_chroms):
    positions, rates = get_positions_rates(num_chroms)

    ## HACK: Discretization hack to avoid overflow for WG
    num_loci = int(positions[-1] / 100)

    recombination_map = msprime.RecombinationMap(
            positions, rates, num_loci
            )

    return recombination_map


def simulate(model=None, num_replicates=1, num_chroms=22):
    if model is None:
        model = 'hudson'

    pc = msprime.PopulationConfiguration(
            initial_size=100,
            sample_size=100,
            growth_rate=0
            )

    rm = get_recombination_map(num_chroms=num_chroms)

    replicates = msprime.simulate(
            population_configurations = [pc],
            num_replicates=num_replicates,
            model=model,
            recombination_map=rm,
            mutation_rate=2e-8
            )

    return replicates


def plot_sfs_list(sfs_list, models, max_count=None, outfile=None):
    if outfile is None:
        outfile = os.path.expanduser('~/temp/test_sfs_plot.png')

    fig, ax = plt.subplots()

    for sfs_replicates, model in zip(sfs_list, models):
        sfs_df = pd.DataFrame(sfs_replicates).fillna(0)
        mean_sfs = sfs_df.mean(axis=0).values

        if max_count is not None:
            mean_sfs = mean_sfs[:max_count + 1]
        x_values = range(1, len(mean_sfs) + 1)
        ax.plot(x_values, mean_sfs, label=model)

    ax.set_xlabel('Frequency')
    ax.set_ylabel('Count')
    # ax.set_xscale('log')
    # ax.set_yscale('log')
    ax.legend()
    fig.savefig(outfile)


def get_sfs(ts):
    counts = []

    for tree in ts.trees():
        for site in tree.sites():
            assert len(site.mutations) == 1
            mutation = site.mutations[0]
            count = tree.num_samples(mutation.node)
            counts.append(count)

    sfs_dict = Counter(counts)
    sfs = [sfs_dict[n] for n in sorted(sfs_dict.keys())]
    print("Singletons:", sfs[0])
    print("Doubletons:", sfs[1])

    return sfs


def simulate_and_plot_sfs(max_count, num_replicates, num_chroms=1):
    sfs_list = []
    models = ['hudson', 'dtwf']
    for model in models:
        print("Simulating", model)
        replicates = simulate(model=model,
                num_replicates=num_replicates,
                num_chroms=num_chroms
                )
        sfs_replicates = [get_sfs(ts) for ts in replicates]
        sfs_list.append(sfs_replicates)

    plot_sfs_list(sfs_list, models, max_count=max_count)


if __name__ == "__main__":
    simulate_and_plot_sfs(max_count=3, num_replicates=20, num_chroms=2)

