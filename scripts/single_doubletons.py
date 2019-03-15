import sys, os
sys.path.append('/home/dnelson/project/msprime/')
sys.path.append('/home/dnelson/project/msprime/lib/subprojects/git-submodules/tskit/python/')
import msprime
import numpy as np
import matplotlib
matplotlib.use('Agg')
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
    # chrom_lengths: typing.List[float] = [1e8]
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


def get_recombination_map(num_chroms, rho=1e-8):
    positions, rates = get_positions_rates(num_chroms, rho)

    ## HACK: Discretization hack to avoid overflow for WG
    num_loci = int(positions[-1] / 100)

    recombination_map = msprime.RecombinationMap(
            positions, rates, num_loci
            )

    return recombination_map


def simulate(model=None, num_replicates=1, num_chroms=22, rho=1e-8):
    if model is None:
        model = 'hudson'

    pc = msprime.PopulationConfiguration(
            initial_size=10000,
            sample_size=10000,
            growth_rate=0
            )

    rm = get_recombination_map(num_chroms=num_chroms, rho=rho)

    replicates = msprime.simulate(
            population_configurations = [pc],
            num_replicates=num_replicates,
            model=model,
            recombination_map=rm,
            # mutation_rate=1,
            mutation_rate=2e-8
            )

    return replicates


def plot_singletons_doubletons(sfs_df, outfile=None):
    if outfile is None:
        outfile = os.path.expanduser('~/temp/test_singletons_doubletons_plot.png')

    max_count = 2

    fig, ax = plt.subplots()

    ## Pull models out of multi-index dataframe
    models, _iterations = sfs_df.index.levels

    i = 0
    width = 0.35
    for model in models:
        df = sfs_df.loc[model]

        mean_sfs = df.mean(axis=0).values
        mean_sfs = mean_sfs[:max_count]
        sem = df.sem(axis=0).values
        sem = sem[:max_count]

        x_vals = np.arange(1, len(mean_sfs) + 1) + i * width
        ax.errorbar(x_vals, mean_sfs,
                        yerr=2*sem,
                        fmt='.',
                        capsize=2,
                        label=model,
                        )
        i += 1

    xTickMarks = ['Singletons', 'Doubletons']
    ax.set_xticks(x_vals - width / 2)
    xtickNames = ax.set_xticklabels(xTickMarks)
    plt.setp(xtickNames, fontsize=10)

    ax.set_ylabel('Count')
    ax.legend()
    fig.savefig(outfile)


def deprecated_plot_sfs(sfs_dict, models, max_count=None, outfile=None):
    if outfile is None:
        outfile = os.path.expanduser('~/temp/test_sfs_plot.png')

    fig, ax = plt.subplots()

    i = 0
    width = 0.2
    for sfs_replicates, model in zip(sfs_dict, models):
        sfs_df = pd.DataFrame(sfs_replicates).fillna(0)
        mean_sfs = sfs_df.mean(axis=0).values
        sem = sfs_df.sem(axis=0).values

        if max_count is not None:
            mean_sfs = mean_sfs[:max_count]
            sem = sem[:max_count]
        x_vals = np.arange(1, len(mean_sfs) + 1) + i * width
        ax.errorbar(x_vals, mean_sfs,
                        yerr=2*sem,
                        fmt='.',
                        capsize=2,
                        label=model,
                        )
        i += 1

    ax.set_xlabel('Frequency')
    ax.set_ylabel('Count')
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

    sfs_counter = Counter(counts)
    sfs = [sfs_counter[n] for n in range(1, max(sfs_counter.keys()) + 1)]
    print("Singletons:", sfs[0], "Doubletons:", sfs[1])

    return sfs


def write_sfs(sfs_df, outfile=None):
    if outfile is None:
        outfile = os.path.expanduser('~/temp/test_write_sfs.txt')

    sfs_df.to_csv(outfile)


def simulate_and_plot_sfs(max_count, num_replicates, num_chroms=1, rho=1e-8):
    dfs = []
    models = ['hudson', 'dtwf']
    for model in models:
        print("Simulating", model)
        replicates = simulate(model=model,
                num_replicates=num_replicates,
                num_chroms=num_chroms,
                rho=rho,
                )
        sfs_replicates = pd.DataFrame([get_sfs(ts) for ts in replicates])
        dfs.append(sfs_replicates)

    sfs_df = pd.concat(dfs, keys=models)
    sfs_df = sfs_df.fillna(0)

    write_sfs(sfs_df, outfile=None)
    plot_singletons_doubletons(sfs_df)


if __name__ == "__main__":
    simulate_and_plot_sfs(max_count=2, num_replicates=2, num_chroms=1, rho=1e-8)

