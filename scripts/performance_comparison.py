import sys, os
sys.path.append('../../msprime')
sys.path.append('../../msprime/lib/subprojects/git-submodules/tskit/python')
import msprime
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import argparse
import pandas as pd
import time
from tqdm import tqdm


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


def get_chrom_lengths(num_chroms, rho=1e-8):
    all_lengths_morgans = [2.77693825, 2.633496065, 2.24483368, 2.12778391, 
            2.03765845, 1.929517394, 1.867959329, 1.701765192, 1.68073935, 
            1.789473882, 1.594854258, 1.72777271, 1.26940475, 1.16331251, 
            1.2554709, 1.348911043, 1.29292106, 1.18978483, 1.077960694, 
            1.079243479, 0.61526812, 0.72706815]

    all_lengths = [x / rho for x in all_lengths_morgans]
    chrom_lengths = all_lengths[:num_chroms]

    return chrom_lengths


def get_recombination_map(num_chroms, rho=1e-8):
    chrom_lengths = get_chrom_lengths(num_chroms, rho)
    positions, rates = get_positions_rates(chrom_lengths, rho)

    ## HACK: Discretization hack to avoid overflow for WG
    num_loci = int(positions[-1] / 10)

    recombination_map = msprime.RecombinationMap(
            positions, rates, num_loci
            )

    return recombination_map


def get_population_configuration(Ne, sample_size):
    population_configurations = []
    p = msprime.PopulationConfiguration(Ne, sample_size, growth_rate=0)
    population_configurations.append(p)

    return population_configurations


def get_sim_kwargs():
    Ne = 100
    sample_size = 100

    population_configuration = get_population_configuration(Ne, sample_size)
    sim_kwargs = {
            'population_configurations': population_configuration,
            }

    return sim_kwargs


def time_simulation(sim_kwargs, model, replicates=10):
    sim_kwargs['model'] = model

    times = []
    for i in tqdm(range(replicates)):
        start_time = time.time()
        ts = msprime.simulate(**sim_kwargs)
        times.append(time.time() - start_time)

    return times


def simulation_times(sim_kwargs, model, max_chroms=22):
    times = []
    for num_chroms in range(1, max_chroms + 1):
        recombination_map = get_recombination_map(num_chroms)
        sim_kwargs['recombination_map'] = recombination_map

        print("Model:", model + ',', "first", num_chroms, "chromosomes:")
        sim_times = time_simulation(sim_kwargs, model)
        times.append(sim_times)

    return times


def main():
    sim_kwargs = get_sim_kwargs()

    hudson_times = simulation_times(sim_kwargs, model='hudson', max_chroms=2)
    dtwf_times = simulation_times(sim_kwargs, model='dtwf', max_chroms=2)


if __name__ == "__main__":
    main()
