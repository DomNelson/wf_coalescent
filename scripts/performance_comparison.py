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
import attr
import typing


@attr.s(auto_attribs=True)
class PerformanceComparison:
    Ne: int
    sample_size: int
    rho: float = 1e-8
    max_chroms: int = 22
    replicates: int = 10
    chrom_lengths: typing.List[float] = [277693825, 2633496065, 224483368,
            212778391, 203765845, 1929517394, 1867959329, 1701765192,
            168073935, 1789473882, 1594854258, 172777271, 126940475, 116331251,
            12554709, 1348911043, 129292106, 118978483, 1077960694, 1079243479,
            61526812, 72706815]

    ## TODO: Could allow passing of extra kwargs for simulations, to allow more
    ## flexible demographic events, for example

    simulated_lengths: typing.List[float] = attr.ib(init=False, default=attr.Factory(list))
    simulated_num_chroms: typing.List[int] = attr.ib(init=False, default=attr.Factory(list))
    simulation_times: dict = attr.ib(init=False, default=attr.Factory(dict))


    def get_positions_rates(self, num_chroms=22):
        """
        Takes a list of chromosome lengths and returns lists of positions and
        rates to pass to msprime.RecombinationMap
        """
        positions = []
        rates = []
        total_length = 0
        for length in self.chrom_lengths[:num_chroms]:
            positions.extend([int(total_length), int(total_length) + int(length) - 1])
            rates.extend([self.rho, 0.5])
            total_length += length

        rates[-1] = 0

        return positions, rates


    def get_recombination_map(self, num_chroms):
        positions, rates = self.get_positions_rates(num_chroms)
        print(positions, rates)

        ## HACK: Discretization hack to avoid overflow for WG
        num_loci = int(positions[-1] / 10)

        recombination_map = msprime.RecombinationMap(
                positions, rates, num_loci
                )

        return recombination_map


    def get_population_configuration(self):
        population_configurations = []
        p = msprime.PopulationConfiguration(self.Ne, self.sample_size)
        population_configurations.append(p)

        return population_configurations


    def get_sim_kwargs(self):
        population_configuration = self.get_population_configuration()
        sim_kwargs = {
                'population_configurations': population_configuration,
                }

        return sim_kwargs


    def time_simulation(self, model, num_chroms):
        recombination_map = self.get_recombination_map(num_chroms)
        sim_kwargs = self.get_sim_kwargs()
        sim_kwargs['recombination_map'] = recombination_map
        sim_kwargs['model'] = model

        times = []
        for i in tqdm(range(self.replicates)):
            start_time = time.time()
            ts = msprime.simulate(**sim_kwargs)
            times.append(time.time() - start_time)

        return times


    def update_simulated_lengths(self, num_chroms):
        assert len(self.simulated_num_chroms) == len(self.simulated_lengths)

        new_length = sum(self.chrom_lengths[:num_chroms])
        self.simulated_num_chroms.append(num_chroms)
        self.simulated_lengths.append(new_length)


    def store_simulation_times(self, model):
        times = []
        for num_chroms in range(1, self.max_chroms + 1):
            print("Model:", model + ',', "first", num_chroms, "chromosomes:")
            sim_times = self.time_simulation(model, num_chroms)
            times.append(sim_times)

            self.update_simulated_lengths(num_chroms)

        times_array = np.array(times)
        self.simulation_times[model] = times_array


    def save(self, outfile):
        np.savez(outfile,
                num_chroms=self.simulated_num_chroms,
                lengths=self.simulated_lengths,
                **self.simulation_times
                )


    def format_for_plot():
        """
        Returns simulation results in the same format as they would be if
        loaded from saved results using np.load(outfile)
        """
        plot_dict = {}
        plot_dict['num_chroms'] = self.simulated_num_chroms
        plot_dict['lengths'] = self.simulated_lengths
        plot_dict.update(self.simulation_times)

        return plot_dict


def plot_times(**kwargs):
    for key, value in kwargs.items():
        print(key, value)


def main():
    outfile = os.path.expanduser('~/temp/times.npz')
    P = PerformanceComparison(
            Ne=100,
            sample_size=100,
            max_chroms=1,
            replicates=10,
            )

    P.store_simulation_times(model='hudson')
    P.store_simulation_times(model='dtwf')

    if outfile is not None:
        P.save(outfile)

    plot_times(dtwf_times=dtwf_times, hudson_times=hudson_times)


if __name__ == "__main__":
    main()
