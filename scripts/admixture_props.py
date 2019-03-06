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


def diploid_gravel_ancestry_variance(T, m, L, K, N):
    """
    T - admixture time
    m - migration rate
    L - total length
    K - number of chromosomes
    N - effective population size
    """
    A = m * (1 - m) / 2 ** (T - 1)
    B_num = 2 * m * (1 - m) * (1 - 1 / (2 * N)) ** (T - 1)
    B_denom = 2 * (K + 1 * (T - 2) * L)
    B = B_num / B_denom

    return A + B


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
            # while True:
            #     try:
            #         sample = next(samples)
            #         sample_copy = next(samples)
            #     except StopIteration:
            #         break
                # len1 = sum(ind_tracts[sample])
                # len2 = sum(ind_tracts[sample_copy])
                # prop = (len1 + len2) / (2 * total_length)
            for sample in samples:
                sample = next(samples)
                prop = sum(ind_tracts[sample])
                prop = prop / total_length
                replicate_props.append(prop)

            ancestry_props.append(replicate_props)
            pbar.update(1)
        
    return ancestry_props


def simulate(args, model, recombination_map, admixture_time):

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
            model=model,
            record_migrations=True,
            num_replicates=args.replicates
    )

    return replicates


def get_output_suffix(args, admixture_time):
    suffix = 'ancestry_variance'
    suffix += 'Ne' + '_' + str(args.Ne) + '_'
    suffix += 'admix_time' + '_' + str(admixture_time) + '_'
    suffix += 'admix_prop' + '_' + str(args.admixture_prop) + '_'
    suffix += 'nchroms' + '_' + str(args.num_chroms)

    return suffix


def set_paper_params(args):
    print("*" * 79)
    print("Reproducing figure from paper - ignoring all args" +\
            " except --out_dir")
    print("*" * 79)
    # admixture_times = [x for x in range(100, 130, 10)]
    admixture_times = [x for x in range(1, 5)]
    admixture_times = [x for x in range(1, 20)]
    admixture_times += [x for x in range(20, 50, 5)]
    admixture_times += [x for x in range(50, 100, 10)]
    admixture_times += [x for x in range(100, 200, 10)]
    admixture_times += [x for x in range(200, 525, 25)]

    paper_args = argparse.Namespace(
            Ne=80,
            sample_size=80,
            dtwf_file=args.dtwf_file,
            hudson_file=args.hudson_file,
            admixture_prop=0.3,
            num_chroms=22,
            replicates=1,
            CI_width=0.95,
            out_dir=args.out_dir,
            discretize_hack=True,
            plot=True,
            )

    return paper_args, admixture_times


def get_output_prefix(out_dir, admixture_times):
    basedir, ext = os.path.splitext(args.out_dir)
    suffix = get_output_suffix(args, admixture_time=admixture_times[-1])

    return os.path.join(basedir, suffix)


def get_simulation_variance(args, model, admixture_times, rec_map):
    nrows = len(admixture_times)
    ncols = args.replicates
    variance_array = np.zeros([nrows, ncols])
    prefix = get_output_prefix(args.out_dir, admixture_times)

    for i, t in enumerate(admixture_times):
        replicates = simulate(args, model,  rec_map, admixture_time=t)
        props_replicates = get_ancestry_props(
                replicates,
                max_time=t,
                num_replicates=args.replicates)

        for j, props in enumerate(props_replicates):
            variance_array[i, j] = np.var(props)# / (1 - 1 / len(props))

    model_df = pd.DataFrame(
            variance_array,
            columns=range(args.replicates),
            index=admixture_times)

    return model_df


def get_mean_variance_with_CIs(model_df, CI_interval=0.95):
    assert 0 < CI_interval <= 1

    low_CI = 0.5 - (CI_interval / 2)
    high_CI = 0.5 + (CI_interval / 2)

    ## Long-winded but other methods don't seem to allow specifying
    ## percentiles
    percentiles = [low_CI, high_CI]
    CI_df = model_df.transpose().describe(percentiles).transpose()

    ## This is awful. Pandas is inconsistent in how it handles trailing
    ## zeros when rounding percentiles and converting to strings for
    ## column names. This seems to handle all cases.
    low_CI_column = str(np.round(low_CI * 100, 1)).strip('0.') + '%'
    high_CI_column = str(np.round(high_CI * 100, 1)).strip('0.') + '%'
    try:
        CI_df = CI_df[['mean', low_CI_column, high_CI_column]]
    except KeyError:
        low_CI_column = str(np.round(low_CI * 100, 1)) + '%'
        CI_df = CI_df[['mean', low_CI_column, high_CI_column]]

    return CI_df


def format_CIs_for_plot(CI_df):
    errs = CI_df.drop(columns='mean').transpose().values
    errs = errs.reshape(1, 2, -1)
    errs[:, 0, :] = CI_df['mean'].transpose().values - errs[:, 0, :]
    errs[:, 1, :] = errs[:, 1, :] - CI_df['mean'].transpose().values

    return errs


def main(args):
    if args.paper_params:
        args, admixture_times = set_paper_params(args)
    else:
        admixture_range = [int(x.strip()) for x in args.admixture_range.split(',')]
        admixture_times = range(*admixture_range)

    prefix = get_output_prefix(args.out_dir, admixture_times)

    positions, rates, num_loci = get_whole_genome_positions_rates_loci(args)
    rec_map = msprime.RecombinationMap(
            positions, rates, num_loci=num_loci)

    dfs = {}

    ## DTWF variance
    if args.dtwf_file is not None:
        df = pd.read_csv(args.dtwf_file, index_col=0)
        dtwf_times = list(df.index)
    else:
        df = get_simulation_variance(args, 'dtwf', admixture_times, rec_map)
        df.to_csv(prefix + '_replicates_dtwf.txt')

    CI_df = get_mean_variance_with_CIs(df, CI_interval=args.CI_width)
    errs = format_CIs_for_plot(CI_df)
    dfs['dtwf'] = [CI_df, errs]


    ## Hudson variance
    if args.hudson_file is not None:
        df = pd.read_csv(args.hudson_file, index_col=0)
        hudson_times = list(df.index)

        ## To help compare to the proper theory curve
        ## TODO: Check other args as well
        assert(hudson_times == dtwf_times)
        admixture_times = hudson_times
    else:
        pass
    #     df = get_simulation_variance(args, 'hudson', admixture_times, rec_map)
    #     df.to_csv(prefix + '_replicates_hudson.txt')
    #
    CI_df = get_mean_variance_with_CIs(df, CI_interval=args.CI_width)
    errs = format_CIs_for_plot(CI_df)
    dfs['hudson'] = [CI_df, errs]

    ## Tracts expected variance
    length_in_morgans = positions[-1] / 1e8
    diploid_gravel_variance = [
            diploid_gravel_ancestry_variance(
                    T,
                    args.admixture_prop,
                    length_in_morgans,
                    args.num_chroms,
                    args.Ne
            ) for T in admixture_times]
    haploid_gravel_variance = [
            haploid_gravel_ancestry_variance(
                    T,
                    args.admixture_prop,
                    length_in_morgans,
                    args.num_chroms,
                    args.Ne
            ) for T in admixture_times]
    print("Comparing vs tracts with", length_in_morgans, "Morgans")
    expected_df = pd.DataFrame(index=admixture_times)
    expected_df['Expected (haploid)'] = haploid_gravel_variance
    # expected_df['Expected (diploid)'] = diploid_gravel_variance

    if args.plot:
        sns.set_palette("muted", 8)
        plot_file = prefix + '.png'

        fig, ax = plt.subplots()

        for model, (CI_df, errs) in dfs.items():
            CI_df['mean'].plot(ax=ax, yerr=errs, capsize=2, fmt='.', legend=False,
                    label=model)
        expected_df.plot(ax=ax, legend=False)

        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.legend()
        print("Plotting to", plot_file)
        fig.savefig(plot_file)

        import IPython; IPython.embed()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--paper_params", action='store_true')
    parser.add_argument("--Ne", type=int, default=80)
    parser.add_argument("--sample_size", type=int, default=80)
    parser.add_argument('--dtwf_file', default=None)
    parser.add_argument('--hudson_file', default=None)
    parser.add_argument('--admixture_range', default="1,10")
    parser.add_argument('--admixture_prop', type=float, default=0.3)
    parser.add_argument('--num_chroms', type=int, default=22)
    parser.add_argument('--replicates', type=int, default=1)
    parser.add_argument('--CI_width', type=float, default=0.66)
    parser.add_argument('--plot', action='store_true')
    parser.add_argument('--discretize_hack', action='store_true')
    parser.add_argument('--out_dir', default=os.path.expanduser('~/temp/'))

    args = parser.parse_args()
    main(args)
