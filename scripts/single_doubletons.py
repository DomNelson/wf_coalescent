import sys, os
sys.path.append('/home/dnelson/project/msprime/')
sys.path.append('/home/dnelson/project/msprime/lib/subprojects/git-submodules/tskit/python/')
import msprime
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from collections import Counter
from tqdm import tqdm


def get_positions_rates(num_chroms=22, rho=1e-8):
    """
    Takes a list of chromosome lengths and returns lists of positions and
    rates to pass to msprime.RecombinationMap
    """
    ## TODO: Check these lengths!
    chrom_lengths: typing.List[float] = [1e6]
    # chrom_lengths: typing.List[float] = [
    #         277693825, 263349606, 224483368, 212778391,
    #         203765845, 192951739, 186795932, 170176519,
    #         168073935, 178947388, 159485425, 172777271,
    #         126940475, 116331251, 12554709, 13489110,
    #         12929210, 11897848, 10779606, 10792434,
    #         61526812, 72706815]

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


def run_hybrid_simulations(sim_kwargs, hybrid_wf_gens):
    """
    TODO: Use this for dtwf sims
    """
    hybrid_kwargs = sim_kwargs.copy()
    hybrid_kwargs['model'] = 'dtwf'
    dtwf_ts = msprime.simulate(
                __tmp_max_time=hybrid_wf_gens,
                **hybrid_kwargs,
                )

    hybrid_kwargs['model'] = 'hudson'
    Ne = hybrid_kwargs['population_configurations'][0].initial_size
    hybrid_kwargs['population_configurations'] = None
    hybrid_ts = msprime.simulate(
            from_ts=dtwf_ts,
            Ne=Ne,
            **hybrid_kwargs,
            )

    return hybrid_ts


def simulate(model=None, num_replicates=1, num_chroms=22, rho=1e-8):
    if model is None:
        model = 'hudson'

    pc = msprime.PopulationConfiguration(
            initial_size=100,
            sample_size=50,
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


def plot_singletons_doubletons(sfs_df, outfile=None, column_rename_dict=None):
    if outfile is None:
        outfile = os.path.expanduser('~/temp/test_sfs_plot.png')

    if column_rename_dict is not None:
        sfs_df = sfs_df.replace(column_rename_dict)

    max_count = 2
    samp = sfs_df[sfs_df['Frequency'] <= max_count]

    fig, ax = plt.subplots(figsize=(5, 3.5))
    sns.violinplot(x='Frequency', y='Count', hue='model', data=samp, ax=ax)

    plt.xticks([0, 1], ['Singletons', 'Doubletons'])
    ax.set_xlabel('')

    fig.savefig(outfile)


def plot_sfs(sfs_df, max_count=None, outfile=None, column_rename_dict=None):
    if outfile is None:
        outfile = os.path.expanduser('~/temp/test2_sfs_plot.png')

    if column_rename_dict is not None:
        sfs_df = sfs_df.replace(column_rename_dict)

    samp = sfs_df[sfs_df['Frequency'] <= max_count]

    fig, ax = plt.subplots(figsize=(5, 3.5))
    sns.violinplot(x='Frequency', y='Count', hue='model', data=samp, ax=ax)
    ax.set_xlabel('Frequency')

    ## Could also set custom xticks here
    # plt.xticks([0, 1, 2], ['Singletons', 'Doubletons', 'Tripletons'])
    # ax.set_xlabel('')

    fig.savefig(outfile)


def get_sfs(ts, model):
    counts = []

    for tree in ts.trees():
        for site in tree.sites():
            assert len(site.mutations) == 1
            mutation = site.mutations[0]
            count = tree.num_samples(mutation.node)
            counts.append(count)

    sfs_counter = Counter(counts)
    sfs = [sfs_counter[n] for n in range(1, max(sfs_counter.keys()) + 1)]
    try:
        print("Singletons:", sfs[0], "Doubletons:", sfs[1])
    except:
        ## Only happens if (in very small sims) we get only singletons
        pass

    df = pd.DataFrame()
    df['Frequency'] = range(1, max(sfs_counter.keys()) + 1)
    df['Count'] = sfs
    df['model'] = model

    return df


def plot_sfs_violin(sfs_df=None, sfs_file=None, plot_outfile=None, max_freq=2):
    if sfs_df is None:
        assert sfs_file is not None
        sfs_df = pd.read_csv(os.path.expanduser(sfs_file), index_col=0)

    if plot_outfile is None:
        plot_outfile = os.path.expanduser('~/temp/test_violin.pdf')

    samp = df[df['Frequency'] <= max_freq]
    sns.violinplot(x='Frequency', y='Count', hue='model', data=samp)


def write_sfs(sfs_df, outfile=None):
    if outfile is None:
        outfile = os.path.expanduser('~/temp/test_sfs.txt')

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
        sfs_replicates = [get_sfs(ts, model) for ts in replicates]
        dfs.extend(sfs_replicates)

    sfs_df = pd.concat(dfs)
    sfs_df = sfs_df.fillna(0)

    write_sfs(sfs_df, outfile=None)
    plot_singletons_doubletons(sfs_df)


if __name__ == "__main__":
    simulate_and_plot_sfs(max_count=2, num_replicates=100,
            num_chroms=1, rho=1e-8)

