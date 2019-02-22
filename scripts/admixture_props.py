import sys, os
sys.path.append('../../msprime')
sys.path.append('../../msprime/lib/subprojects/git-submodules/tskit/python')
import msprime
import collections
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import argparse
import pandas as pd
from tqdm import tqdm


def gravel_ancestry_variance(T, m, L, K, N):
    """
    T - admixture time
    m - migration rate
    L - total length
    K - number of chromosomes
    N - effective population size
    """
    A = m * (1 - m) / 2 ** (T - 1)
    B_num = 2 * m * (1 - m) * (1 - 1 / (2 * N)) ** (T - 1)
    B_denom = 2 * K + 2 * (T - 2) * L
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


def get_ind_tracts(ts, max_time):
    ind_tracts = collections.defaultdict(list)
    
    key = lambda x: x.left
    migrations = sorted(ts.migrations(), key=key)
    trees = ts.trees()
    
    for migration in migrations:
        if migration.time > max_time:
            continue
            
        stored = False
        node = migration.node
        length = migration.right - migration.left

        while not stored:
            t = next(trees)
            try:
                samples = t.get_leaves(node)
            except ValueError:
                continue
                
            for s in samples:
                ind_tracts[s].append(length)
            stored = True

    # import IPython; IPython.embed()
                
    return ind_tracts


def get_ancestry_props(replicates, max_time, num_replicates):
    ancestry_props = []

    with tqdm(total=num_replicates, desc=str(max_time)) as pbar:
        for ts in replicates:
            # from IPython import embed; embed()
            # sys.exit()
            ind_tracts = get_ind_tracts(ts, max_time)
            total_length = ts.get_sequence_length()

            replicate_props = []
            samples = iter(ts.samples())
            for sample in samples:
                sample_copy = next(samples)
                ## Convert to Morgans
                prop = sum(ind_tracts[sample]) / total_length
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

    migration_matrix = None
    # [
    #     [0, args.admixture_prop],
    #     [0,       0],
    # ]

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
            migration_matrix=migration_matrix,
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


def main(args):

    if args.paper_params:
        print("*" * 79)
        print("Reproducing figure from paper - ignoring all other args" +\
                " except --model and --out_dir")
        print("*" * 79)
        admixture_times = [x for x in [1, 2]]
        # admixture_times = [x for x in range(1, 20)]
        # admixture_times += [x for x in range(20, 50, 4)]
        # admixture_times += [x for x in range(50, 100, 10)]
        # admixture_times += [x for x in range(100, 200, 10)]
        # admixture_times += [x for x in range(200, 500, 25)]

        args = argparse.Namespace(
                Ne=100,
                sample_size=100,
                model=args.model,
                admixture_prop=0.3,
                num_chroms=22,
                replicates=50,
                out_dir=args.out_dir,
                discretize_hack=True,
                plot=True,
                )
    else:
        admixture_range = [int(x.strip()) for x in args.admixture_range.split(',')]
        admixture_times = range(*admixture_range)


    rho = 1e-8
    all_lengths = [247249719, 242951149, 199501827, 191273063, 180857866,
            170899992, 158821424, 146274826, 140273252, 135374737, 134452384,
            132349534, 114142980, 106368585, 100338915, 88827254, 78774742,
            76117153, 63811651, 62435964, 46944323, 49691432]
    chrom_lengths = all_lengths[:args.num_chroms]

    positions, rates = get_positions_rates(chrom_lengths, rho)
    num_loci = positions[-1]

    ## HACK: This is to avoid bad edge intervals (right <= left) when
    ## simulating whole genomes. Possible overflow / floating point precision
    ## issue
    if args.discretize_hack:
        num_loci = int(num_loci / 100)
        print("Scaling to", num_loci, "loci")

    recombination_map = msprime.RecombinationMap(
            positions, rates, num_loci=num_loci)
    
    nrows = len(admixture_times)
    ncols = args.replicates
    variance_array = np.zeros([nrows, ncols])
    bins = np.arange(0, 1.1, 0.1)

    df = pd.DataFrame()
    for i, t in enumerate(admixture_times):
        replicates = simulate(args, recombination_map, admixture_time=t)
        props_replicates = get_ancestry_props(
                replicates,
                max_time=t,
                num_replicates=args.replicates)

        for j, props in enumerate(props_replicates):
            variance_array[i, j] = np.var(props)

    df = pd.DataFrame(
            variance_array,
            columns=range(args.replicates),
            index=admixture_times)
    average_variance = df.mean(axis=1)

    basedir, ext = os.path.splitext(args.out_dir)
    suffix = get_output_suffix(args, admixture_time=t)
    average_variance.to_csv(os.path.join(basedir, suffix + '.txt'))

    if args.plot:

        length_in_morgans = positions[-1] / 1e8
        gravel_variance = [
                gravel_ancestry_variance(
                        T,
                        args.admixture_prop,
                        length_in_morgans,
                        args.num_chroms,
                        args.Ne
                ) for T in admixture_times]
        print("Comparing vs tracts with", length_in_morgans, "Morgans")

        plot_file = os.path.join(basedir, suffix + '.png')

        fig, ax = plt.subplots()
        ax.plot(admixture_times, gravel_variance, label='Expected')
        average_variance.plot(ax=ax, label=str(args.model))
        ax.set_xscale('log')
        ax.set_yscale('log')
        fig.legend()
        fig.savefig(plot_file)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--paper_params", action='store_true')
    parser.add_argument("--Ne", type=int, default=80)
    parser.add_argument("--sample_size", type=int, default=80)
    parser.add_argument('--model', default='Hudson')
    parser.add_argument('--admixture_range', default="1,10")
    parser.add_argument('--admixture_prop', type=float, default=0.3)
    parser.add_argument('--num_chroms', type=int, default=22)
    parser.add_argument('--replicates', type=int, default=1)
    parser.add_argument('--plot', action='store_true')
    parser.add_argument('--discretize_hack', action='store_true')
    parser.add_argument('--out_dir', default=os.path.expanduser('~/temp/'))

    args = parser.parse_args()
    main(args)
