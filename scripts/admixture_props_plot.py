import sys, os
sys.path.append('../../msprime')
import msprime
import collections
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import argparse


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
                
    return ind_tracts


def get_ancestry_props(replicates, max_time):
    ancestry_props = []

    for ts in replicates:
        ind_tracts = get_ind_tracts(ts, max_time)
        total_length = ts.get_sequence_length()

        replicate_props = []
        for sample in ts.samples():
            prop = sum(ind_tracts[sample]) / total_length
            replicate_props.append(prop)

        ancestry_props.append(replicate_props)
        
    return ancestry_props


def simulate(args):
    rho = 1e-8
    all_lengths = [247249719, 242951149, 199501827, 191273063, 180857866,
            170899992, 158821424, 146274826, 140273252, 135374737, 134452384,
            132349534, 114142980, 106368585, 100338915, 88827254, 78774742,
            76117153, 63811651, 62435964, 46944323, 49691432]
    chrom_lengths = all_lengths[:args.num_chroms]
    num_loci = chrom_lengths[-1] + 1

    positions, rates = get_positions_rates(chrom_lengths, rho)
    recombination_map = msprime.RecombinationMap(
            positions, rates, num_loci=num_loci)

    population_configurations = [
            msprime.PopulationConfiguration(
                    sample_size=args.Ne,
                    initial_size=args.Ne
            ),
            msprime.PopulationConfiguration(
                    sample_size=0,
                    initial_size=args.Ne
            )
    ]

    demographic_events = [
        msprime.MassMigration(
            time=args.admixture_time,
            source=0,
            dest=1,
            proportion=args.admixture_prop
        ),
        msprime.MassMigration(
            time=args.admixture_time + 1,
            source=1,
            dest=0,
            proportion=1
        ),
        msprime.PopulationParametersChange(
            time=args.admixture_time + 2,
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


def get_output_suffix(args):
    suffix = 'admixture_props_'
    suffix += 'Ne' + '_' + str(args.Ne) + '_'
    suffix += 'model' + '_' + str(args.model) + '_'
    suffix += 'admix_time' + '_' + str(args.admixture_time) + '_'
    suffix += 'admix_prop' + '_' + str(args.admixture_prop) + '_'
    suffix += 'nchroms' + '_' + str(args.num_chroms)

    return suffix


def main(args):
    replicates = simulate(args)
    props_replicates = get_ancestry_props(replicates, args.admixture_time)

    bins = np.arange(0, 1.1, 0.1)
    for props in props_replicates:
        print("Ancestry mean:", np.mean(props), "variance:", np.var(props))
        plt.hist(props, alpha=0.3, bins=bins)

    basename, ext = os.path.splitext(args.out_dir)
    suffix = get_output_suffix(args)
    plot_file = basename + suffix + '.png'
    plt.savefig(plot_file)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--Ne", type=int, default=1000)
    parser.add_argument('--model', default='Hudson')
    parser.add_argument('--admixture_time', type=float, default=1)
    parser.add_argument('--admixture_prop', type=float, default=0.2)
    parser.add_argument('--num_chroms', type=int, default=1)
    parser.add_argument('--replicates', type=int, default=1)
    parser.add_argument('--out_dir', default='admixture_props.txt')

    args = parser.parse_args()
    main(args)
