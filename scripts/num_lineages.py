import sys, os
import numpy as np
import subprocess
import argparse
from datetime import datetime

msprime_dir = os.path.expanduser('../../msprime')
sys.path.insert(1, msprime_dir)
sys.path.insert(1, msprime_dir + '/lib/subprojects/git-submodules/tskit/python')
import msprime

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


def srun_job(args):
    out_dir = os.path.expanduser(args.out_dir)
    timestamp = datetime.strftime(datetime.now(), '%Y-%m-%d_%H-%M-%S')

    slurm_cmd = 'srun'
    job_cmd = 'python num_lineages.py'
    out_file = os.path.join(out_dir, timestamp)
    for arg, value in args.__dict__.items():
        if arg in ['out_dir']:
            continue

        if arg in ['srun', 'mem', 'time']:
            slurm_cmd += ' --' + str(arg) + ' ' + str(value)
        else:
            out_file += '_' + str(arg) + str(value)

        cmd += ' --' + str(arg) + ' ' + str(value)

    ## TODO: Need to submit to srun with proper args, which haven't been
    ## included yet
    out_file += '.gz'
    cmd = slurm_cmd + ' ' + job_cmd
    cmd += ' | gzip > ' + out_file
    cmd += ' &'

    print(cmd)
    sys.exit()


def simulate_num_lineages(args):
    Ne = args.Ne
    model = args.model
    sample_size = args.sample_size
    rho = 1e-8

    chrom_lengths = [247249719, 242951149, 199501827, 191273063, 180857866,
            170899992, 158821424, 146274826, 140273252, 135374737, 134452384,
            132349534, 114142980, 106368585, 100338915, 88827254, 78774742,
            76117153, 63811651, 62435964, 46944323, 49691432]
    chrom_lengths = chrom_lengths[:args.nchroms]
    # chrom_lengths = [1e8] * 10
    # print(sum(chrom_lengths) / 1e6)
    # sys.exit()

    positions, rates = get_positions_rates(chrom_lengths, rho)
    # num_loci = positions[-1] - 1
    num_loci = int(positions[-1] / 100) ## Avoids overflow for whole-genome sims
    recombination_map = msprime.RecombinationMap(
            positions, rates, num_loci=num_loci
    )

    seed = np.random.randint(1e9)
    ts = msprime.simulate(Ne=Ne,
            recombination_map=recombination_map,
            sample_size=sample_size,
            model=model,
            random_seed=seed
            )


def main(args):
    if args.srun is not None:
        srun_job(args)
    else:
        simulate_num_lineages(args)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--Ne', default=50, type=int)
    parser.add_argument('--sample_size', default=50, type=int)
    parser.add_argument('--model', default='Hudson')
    parser.add_argument('--nchroms', default=22, type=int)
    parser.add_argument('--srun', action="store_true")
    parser.add_argument('--mem', default="7000")
    parser.add_argument('--time', default="3:00:00")
    parser.add_argument('--out_dir', default='~/temp/')

    args = parser.parse_args()

    main(args)


