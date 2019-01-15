import sys, os
import numpy as np
sys.path.append(os.path.expanduser('~/project/msprime'))
import msprime
import subprocess
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


def main(args):
    Ne = args.Ne
    model = args.model
    sample_size = Ne
    rho = 1e-7

    chrom_lengths = [247249719, 242951149, 199501827, 191273063, 180857866,
            170899992, 158821424, 146274826, 140273252, 135374737, 134452384,
            132349534, 114142980, 106368585, 100338915, 88827254, 78774742,
            76117153, 63811651, 62435964, 46944323, 49691432]
    num_loci = chrom_lengths[-1] + 1

    # chrom_lengths = [1e8] * 10

    positions, rates = get_positions_rates(chrom_lengths, rho)
    recombination_map = msprime.RecombinationMap(
            positions, rates, num_loci=num_loci
    )
    for p, r in zip(positions, rates):
        print(p, r)

    for i in range(args.niter):
        seed = np.random.randint(1e9)
        print("Simulating", model, "with random seed", seed)
        try:
            ts = msprime.simulate(Ne=Ne,
                    recombination_map=recombination_map,
                    sample_size=sample_size,
                    model=model,
                    random_seed=seed
                    )
        except Exception as e:
            print("Error!")
            with open('multiple_chroms_err.txt', 'a') as f:
                f.write("seed: " + str(seed) + '\n')
                f.write(str(args) + '\n')
                f.write(str(e) + '\n')
                f.write('=' * 60 + '\n')



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--Ne', default=50, type=int)
    parser.add_argument('--model', default='Hudson')
    parser.add_argument('--niter', default=50, type=int)

    args = parser.parse_args()
    main(args)


