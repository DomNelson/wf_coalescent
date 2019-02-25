import sys, os
sys.path.append('../../msprime')
sys.path.append('../../msprime/lib/subprojects/git-submodules/tskit/python')
import msprime
import collections
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import argparse
import pandas as pd
from tqdm import tqdm


def haploid_gravel_ancestry_variance(T, m, L, K, N):
    """
    T - admixture time
    m - migration rate
    L - total length
    K - number of chromosomes
    N - effective population size
    """
    A = m * (1 - m) / 2 ** (T - 1)
    B_num = 2 * m * (1 - m) * (1 - 1 / (2 * N)) ** (T - 1)
    B_denom = 1 * (K + 1 * (T - 2) * L)
    B = B_num / B_denom
    
    return A + B


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


def get_whole_genome_positions_rates_loci(args, rho=1e-8):
    all_lengths_morgans = [2.77693825, 2.633496065, 2.24483368, 2.12778391, 
            2.03765845, 1.929517394, 1.867959329, 1.701765192, 1.68073935, 
            1.789473882, 1.594854258, 1.72777271, 1.26940475, 1.16331251, 
            1.2554709, 1.348911043, 1.29292106, 1.18978483, 1.077960694, 
            1.079243479, 0.61526812, 0.72706815]
    all_lengths = [x * 1e8 for x in all_lengths_morgans]

    chrom_lengths = all_lengths[:args.num_chroms]

    positions, rates = get_positions_rates(chrom_lengths, rho)
    num_loci = positions[-1]

    ## HACK: This is to avoid bad edge intervals (right <= left) when
    ## simulating whole genomes. Possible overflow / floating point precision
    ## issue
    if args.discretize_hack:
        num_loci = int(num_loci / 100)
        print("Scaling to", num_loci, "loci")

    return positions, rates, num_loci


def get_ind_tracts(ts, max_time):
    ind_tracts = collections.defaultdict(list)
    
    key = lambda x: x.left
    migrations = sorted(ts.migrations(), key=key)
    trees = ts.trees()
    t = next(trees)
    
    for migration in migrations:
        if migration.time > max_time:
            continue
            
        node = migration.node
        length = migration.right - migration.left

        while t.interval[1] <= migration.left:
            t = next(trees)

        assert t.interval[0] <= migration.left and t.interval[1] > migration.left

        samples = t.get_leaves(node)
        for s in samples:
            ind_tracts[s].append(length)

    return ind_tracts


def get_ancestry_props(replicates, max_time, num_replicates):
    ancestry_props = []

    with tqdm(total=num_replicates, desc=str(max_time)) as pbar:
        for ts in replicates:
            ind_tracts = get_ind_tracts(ts, max_time)
            total_length = ts.get_sequence_length()

            replicate_props = []

            samples = iter(ts.samples())
            for sample in samples:
                sample = next(samples)
                prop = sum(ind_tracts[sample])
                prop = prop / total_length
                replicate_props.append(prop)

            ancestry_props.append(replicate_props)
            pbar.update(1)
        
    return ancestry_props


def simulate(args, recombination_map, admixture_time):

    population_configurations = [
            msprime.PopulationConfiguration(
                    sample_size=args.sample_size,
                    initial_size=args.Ne
            ),
            msprime.PopulationConfiguration(
                    sample_size=0,
                    initial_size=args.Ne
            )
    ]

    demographic_events = [
        msprime.MassMigration(
            time=admixture_time,
            source=0,
            dest=1,
            proportion=args.admixture_prop
        ),
        msprime.MigrationRateChange(time=admixture_time, rate=0),
        msprime.MassMigration(
            time=admixture_time + 1,
            source=1,
            dest=0,
            proportion=1
        ),
        msprime.PopulationParametersChange(
            time=admixture_time + 2,
            initial_size=1.0,
            population_id=0
        )
    ]

    replicates = msprime.simulate(
            recombination_map=recombination_map,
            demographic_events=demographic_events,
            population_configurations=population_configurations,
            model=args.model,
            record_migrations=True,
            num_replicates=args.replicates
    )

    return replicates


def get_output_suffix(args, admixture_time):
    suffix = 'ancestry_variance'
    suffix += 'Ne' + '_' + str(args.Ne) + '_'
    suffix += 'model' + '_' + str(args.model) + '_'
    suffix += 'admix_time' + '_' + str(admixture_time) + '_'
    suffix += 'admix_prop' + '_' + str(args.admixture_prop) + '_'
    suffix += 'nchroms' + '_' + str(args.num_chroms)

    return suffix


def set_paper_params(args):
    print("*" * 79)
    print("Reproducing figure from paper - ignoring all args" +\
            " except --out_dir")
    print("*" * 79)
    # admixture_times = [x for x in range(1, 5)]
    admixture_times = [x for x in range(1, 20)]
    admixture_times += [x for x in range(20, 50, 5)]
    admixture_times += [x for x in range(50, 100, 10)]
    admixture_times += [x for x in range(100, 200, 10)]
    admixture_times += [x for x in range(200, 500, 25)]

    paper_args = argparse.Namespace(
            Ne=80,
            sample_size=80,
            model=None,
            admixture_prop=0.3,
            num_chroms=22,
            replicates=3,
            out_dir=args.out_dir,
            discretize_hack=True,
            plot=True,
            )

    return paper_args, admixture_times


def get_output_prefix(out_dir, admixture_times):
    basedir, ext = os.path.splitext(args.out_dir)
    suffix = get_output_suffix(args, admixture_time=admixture_times[-1])

    return os.path.join(basedir, suffix)


def get_simulation_variance(args, admixture_times, recombination_map):
    nrows = len(admixture_times)
    ncols = args.replicates
    variance_array = np.zeros([nrows, ncols])
    prefix = get_output_prefix(args.out_dir, admixture_times)

    var_df = pd.DataFrame(index=admixture_times)
    CIs = []
    for model in ['dtwf', 'hudson']:
        args.model = model
        for i, t in enumerate(admixture_times):
            replicates = simulate(args, recombination_map, admixture_time=t)
            props_replicates = get_ancestry_props(
                    replicates,
                    max_time=t,
                    num_replicates=args.replicates)

            for j, props in enumerate(props_replicates):
                variance_array[i, j] = np.var(props) / (1 - 1 / len(props))

        model_df = pd.DataFrame(
                variance_array,
                columns=range(args.replicates),
                index=admixture_times)
        average_variance = model_df.mean(axis=1)
        var_df[model] = average_variance
        var_df.to_csv(prefix + '.txt')

        ## Long-winded but other methods don't seem to allow specifying
        ## percentiles
        # model_CI = model_df.transpose().describe([0.025, 0.975]).transpose()
        # CIs.append(model_CI[['2.5%', '97.5%']].transpose().values)
        model_CI = model_df.transpose().describe([0.165, 0.835]).transpose()
        CIs.append(model_CI[['16.5%', '83.5%']].transpose().values)

    errs = np.stack(CIs, axis=1)

    return var_df, errs


def main(args):
    if args.paper_params:
        args, admixture_times = set_paper_params(args)
    else:
        admixture_range = [int(x.strip()) for x in args.admixture_range.split(',')]
        admixture_times = range(*admixture_range)

    if args.model is None:
        models = ['dtwf', 'hudson']
    else:
        models = [x.strip() for x in args.model.split(',')]

    positions, rates, num_loci = get_whole_genome_positions_rates_loci(args)
    rec_map = msprime.RecombinationMap(
            positions, rates, num_loci=num_loci)

    var_df, errs = get_simulation_variance(args, admixture_times, rec_map)

    prefix = get_output_prefix(args.out_dir, admixture_times)
    length_in_morgans = positions[-1] / 1e8
    haploid_gravel_variance = [
            haploid_gravel_ancestry_variance(
                    T,
                    args.admixture_prop,
                    length_in_morgans,
                    args.num_chroms,
                    args.Ne
            ) for T in admixture_times]
    print("Comparing vs tracts with", length_in_morgans, "Morgans")
    var_df['Expected (haploid)'] = haploid_gravel_variance
    var_df.to_csv(prefix + '.txt')

    if args.plot:
        sns.set_palette("muted", 8)
        plot_file = prefix + '.png'

        fig, ax = plt.subplots()

        sim_var_df = var_df.drop(columns='Expected (haploid)')
        sim_var_df.plot(ax=ax, yerr=errs, capsize=2, fmt='o', legend=False)
        var_df['Expected (haploid)'].plot(ax=ax, legend=False)

        ax.set_xscale('log')
        ax.set_yscale('log')
        fig.legend()
        fig.savefig(plot_file)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--paper_params", action='store_true')
    parser.add_argument("--Ne", type=int, default=80)
    parser.add_argument("--sample_size", type=int, default=80)
    parser.add_argument('--model', default=None)
    parser.add_argument('--admixture_range', default="1,10")
    parser.add_argument('--admixture_prop', type=float, default=0.3)
    parser.add_argument('--num_chroms', type=int, default=22)
    parser.add_argument('--replicates', type=int, default=1)
    parser.add_argument('--plot', action='store_true')
    parser.add_argument('--discretize_hack', action='store_true')
    parser.add_argument('--out_dir', default=os.path.expanduser('~/temp/'))

    args = parser.parse_args()
    main(args)
