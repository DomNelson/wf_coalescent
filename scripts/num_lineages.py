import sys, os
import numpy as np
import subprocess
import argparse
from datetime import datetime

# Kind of a dirty hack but it works for now
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


def run_sbatch_job(args):
    job_script = 'job_script.sh'
    out_dir = args.out_dir
    timestamp = datetime.strftime(datetime.now(), '%Y-%m-%d_%H-%M-%S')

    header = "#!/bin/bash\n#SBATCH --mem={}mb\n#SBATCH --time={}\n\n".format(
            args.mem, args.time)

    cmd = 'python num_lineages.py'
    out_file = os.path.join(out_dir, 'num_lineages')
    for arg, value in args.__dict__.items():
        if value is None:
            continue

        if arg not in ['srun', 'out_dir', 'mem', 'time']:
            out_file += '_' + str(arg) + str(value)
            cmd += ' --' + str(arg) + ' ' + str(value)

    out_file += '_' + timestamp + '.txt'
    cmd += ' | python mask_gens.py > ' + out_file

    with open(job_script, 'w') as f:
        f.write(header)
        f.write(cmd)
        f.write('\n')

    print("Submitting:")
    print('\n', cmd, '\n')
    subprocess.run('sbatch job_script.sh', shell=True)


def simulate_num_lineages(args):
    Ne = args.Ne
    model = args.model
    sample_size = args.sample_size
    growth_rate = args.growth_rate
    rho = 1e-8

    # chrom_lengths = [247249719, 242951149, 199501827, 191273063, 180857866,
    #         170899992, 158821424, 146274826, 140273252, 135374737, 134452384,
    #         132349534, 114142980, 106368585, 100338915, 88827254, 78774742,
    #         76117153, 63811651, 62435964, 46944323, 49691432]
    chrom_lengths_morgans = [2.77693825, 2.633496065, 2.24483368, 
		2.12778391, 2.03765845, 1.929517394, 
		1.867959329, 1.701765192, 1.68073935, 
		1.789473882, 1.594854258, 1.72777271, 
		1.26940475, 1.16331251, 1.2554709, 
		1.348911043, 1.29292106, 1.18978483, 
		1.077960694, 1.079243479, 0.61526812, 
		0.72706815]

    chrom_lengths = [l * 1e8 for l in chrom_lengths_morgans]
    chrom_lengths = chrom_lengths[:args.nchroms]
    # print(sum(chrom_lengths))

    # chrom_lengths = [1e8] * 10
    # print(sum(chrom_lengths) / 1e6)
    # sys.exit()

    positions, rates = get_positions_rates(chrom_lengths, rho)
    # num_loci = positions[-1] - 1
    num_loci = int(positions[-1] / 100) ## Avoids overflow for whole-genome sims
    recombination_map = msprime.RecombinationMap(
            positions, rates, num_loci=num_loci
    )

    demographic_events = None
    if args.min_pop_size is not None:
        time_to_min_pop_size = (np.log(Ne) - np.log(args.min_pop_size)) / growth_rate
        if time_to_pop_size_100 > 0:
            demographic_events = [
                    msprime.PopulationParametersChange(
                        time=time_to_pop_size_100,
                        growth_rate=0
                        )
                    ]

    population_configurations = [
            msprime.PopulationConfiguration(
                sample_size=sample_size,
                initial_size=Ne,
                growth_rate=growth_rate
                )
            ]

    seed = np.random.randint(1e9)
    ts = msprime.simulate(
            population_configurations=population_configurations,
            demographic_events=demographic_events,
            recombination_map=recombination_map,
            model=model,
            random_seed=seed
        )


def main(args):
    if args.srun is True:
        run_sbatch_job(args)
    else:
        simulate_num_lineages(args)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--Ne', default=50, type=int)
    parser.add_argument('--min_pop_size', type=int)
    parser.add_argument('--sample_size', default=50, type=int)
    parser.add_argument('--growth_rate', default=0.005, type=float)
    parser.add_argument('--model', default='Hudson')
    parser.add_argument('--nchroms', default=22, type=int)
    parser.add_argument('--srun', action="store_true")
    parser.add_argument('--mem', default="7000")
    parser.add_argument('--time', default="3:00:00")
    parser.add_argument('--out_dir', default='~/temp/')

    args = parser.parse_args()

    if args.min_pop_size is not None:
        assert (args.min_pop_size <= args.Ne)

    main(args)
