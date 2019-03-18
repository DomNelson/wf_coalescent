import sys, os
sys.path.append('../../msprime')
sys.path.append('../../msprime/lib/subprojects/git-submodules/tskit/python')
import msprime
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from collections import defaultdict
import argparse
import pandas as pd
import time
from tqdm import tqdm
import attr
import typing
import subprocess
# from cycler import cycler
# mpl.rcParams['axes.prop_cycle'] = cycler(color='bgrcmyk')


@attr.s(auto_attribs=True)
class PerformanceComparison:
    Ne: int
    sample_size: int
    rho: float = 1e-8
    max_chroms: int = 22
    replicates: int = 10
    hybrid_wf_gens: typing.List[int] = attr.ib(
            default=attr.Factory(list)
            )

    ## Workaround to simplify multiple hybrid simulations
    _temp_hybrid_wf_gens: int = attr.ib(init=False, default=-1)

    ## Lengths from GRCh38.p12
    chrom_lengths: typing.List[float] = [
	    248956422, 242193529, 198295559, 190214555, 181538259,
	    170805979, 159345973, 145138636, 138394717, 133797422,
	    135086622, 133275309, 114364328, 107043718, 101991189,
	    90338345, 83257441, 80373285, 58617616, 64444167,
	    46709983, 50818468]

    ## TODO: Could allow passing of extra kwargs for simulations, to allow more
    ## flexible demographic events, for example

    ## Complicated way of defining uninitialized empty collections
    simulated_lengths: typing.List[float] = attr.ib(
            init=False, default=attr.Factory(list)
            )
    simulated_num_chroms: typing.List[int] = attr.ib(
            init=False, default=attr.Factory(list)
            )
    simulation_times: dict = attr.ib(init=False, default=defaultdict(list))


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
        num_loci = int(positions[-1] / 100)

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


    def run_hybrid_simulations(self, sim_kwargs):
        assert len(self.hybrid_wf_gens) > 0
        assert self._temp_hybrid_wf_gens > 0

        hybrid_kwargs = sim_kwargs.copy()
        hybrid_kwargs['model'] = 'dtwf'
        dtwf_ts = msprime.simulate(
                    __tmp_max_time=self._temp_hybrid_wf_gens,
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

        # import IPython; IPython.embed()

        return hybrid_ts


    def time_simulation(self, model, num_chroms):
        recombination_map = self.get_recombination_map(num_chroms)
        sim_kwargs = self.get_sim_kwargs()
        sim_kwargs['recombination_map'] = recombination_map
        sim_kwargs['model'] = model

        times = []
        for i in tqdm(range(self.replicates)):
            ## This isn't the cleanest timing but we don't need super high
            ## precision here
            start_time = time.time()

            if model.lower() == 'hybrid':
                ts = self.run_hybrid_simulations(sim_kwargs)
            else:
                ts = msprime.simulate(**sim_kwargs)

            times.append(time.time() - start_time)

        return np.array(times)


    def update_simulated_lengths(self, num_chroms):
        assert len(self.simulated_num_chroms) == len(self.simulated_lengths)

        new_length = sum(self.chrom_lengths[:num_chroms])
        self.simulated_num_chroms.append(num_chroms)
        self.simulated_lengths.append(new_length)


    def initialize(self, models):
        ## These lists should be empty
        assert len(self.simulated_lengths) == 0
        assert len(self.simulated_num_chroms) == 0

        ## This dict should be empty
        assert len(self.simulation_times) == 0


    def store_simulation_times(self, models):
        self.initialize(models)

        for num_chroms in range(1, self.max_chroms + 1):
            for model in models:
                ## Not the cleanest but it'll do
                if model.lower() == 'hybrid':
                    for t in self.hybrid_wf_gens:
                        self._temp_hybrid_wf_gens = t
                        print("Model:", model + ',',
                                "with", t, "WF generations,",
                                "first", num_chroms, "chromosomes:")

                        sim_times = self.time_simulation(model, num_chroms)
                        label = 'hybrid_' + str(t)
                        self.simulation_times[label].append(sim_times)

                        ## Reset the temp variable
                        self._temp_hybrid_wf_gens = -1
                else:
                    print("Model:", model + ',', "first",
                            num_chroms, "chromosomes:")
                    sim_times = self.time_simulation(model, num_chroms)
                    self.simulation_times[model].append(sim_times)

            self.update_simulated_lengths(num_chroms)

        ##TODO: Could rewrite these as filling in pre-allocated array


    def save(self, outfile):
        np.savez(outfile,
                num_chroms=self.simulated_num_chroms,
                lengths=self.simulated_lengths,
                **self.simulation_times
                )


    def format_for_plot(self):
        """
        Returns simulation results in the same format as they would be if
        loaded from saved results using np.load(outfile)
        """
        plot_dict = {}
        plot_dict['num_chroms'] = self.simulated_num_chroms
        plot_dict['lengths'] = self.simulated_lengths
        plot_dict.update(self.simulation_times)

        return plot_dict


def plot_times(plotfile=None, **sim_times):
    lengths = sim_times.pop('lengths')
    num_chroms = sim_times.pop('num_chroms')

    fig, ax = plt.subplots()

    models = sorted(sim_times.keys())
    for model in models:
        times_list = sim_times[model]
        mean_times = [np.mean(times) for times in times_list]
        ax.plot(lengths, mean_times, label=model)

    cur_x_low, cur_x_high = ax.get_xlim()
    ax.set_xlim([0, cur_x_high])
    ax.set_xlabel('Simulated length (base pairs)')
    ax.set_ylabel('Time (s)')
    ax.legend()
    fig.savefig(plotfile)


def submit_qsub_script(args):
    tmp_file = '.qsub.sh'

    header = '#!/bin/bash\n'
    header += '#PBS -l walltime=' + str(args.walltime) + '\n'
    header += '#PBS -l nodes=1:ppn=' + str(args.nprocs) + '\n'
    header += '\n'

    cmd = 'python ' + os.path.realpath(__file__) + ' '
    for arg, value in args._get_kwargs():
        ## For regular kwargs
        if value is not None and type(value) is not bool:
            cmd += '--' + arg + ' ' + str(value) + ' '

        ## For action = 'store_true'
        if value is True:
            cmd += '--' + arg + ' '

    cmd += '\n'

    with open(tmp_file, 'w') as f:
        f.write(header)
        f.write(cmd)

    ret = subprocess.run(['qsub', tmp_file], check=True)

    import IPython; IPython.embed()


def main(args):
    if args.qsub:
        assert args.walltime is not None
        assert args.nprocs is not None
        submit_qsub_script(args)
        print("Job submitted!")
        sys.exit()


    outfile = os.path.expanduser('~/temp/times_hybrid2.npz')
    if args.outfile is not None:
        outfile = os.path.expanduser(args.outfile)

    plotfile  = os.path.expanduser('~/temp/times_hybrid2_plot.png')
    if args.plotfile is not None:
        plotfile  = os.path.expanduser(args.plotfile)

    models = [x.lower() for x in args.models.strip().split(',')]
    if 'hybrid' in models:
        assert args.hybrid_wf_gens is not None

    # loaded = np.load(outfile)
    # plot_times(plotfile, **loaded)
    short_chroms = [1e6] * 22

    hybrid_wf_gens = None
    if args.hybrid_wf_gens is not None:
        hybrid_wf_gens = [int(x) for x in args.hybrid_wf_gens.strip().split(',')]

    P = PerformanceComparison(
            # chrom_lengths=short_chroms,
            Ne=args.Ne,
            sample_size=args.sample_size,
            max_chroms=args.max_chroms,
            replicates=args.replicates,
            hybrid_wf_gens=hybrid_wf_gens,
            )

    P.store_simulation_times(models=models)

    if outfile is not None:
        P.save(outfile)

    pd = P.format_for_plot()

    if plotfile is not None:
        plot_times(plotfile, **pd)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--paper_params", action='store_true')
    parser.add_argument("--Ne", type=int, default=100)
    parser.add_argument("--sample_size", type=int, default=100)
    parser.add_argument('--max_chroms', type=int, default=5)
    parser.add_argument('--replicates', type=int, default=1)
    parser.add_argument('--plotfile', default=None)
    parser.add_argument('--discretize_hack', action='store_true')
    parser.add_argument('--hybrid_wf_gens', default=None)
    parser.add_argument('--models', required=True)
    parser.add_argument('--outfile', required=True)
    parser.add_argument('--qsub', action='store_true')
    parser.add_argument('--walltime')
    parser.add_argument('--nprocs')

    args = parser.parse_args()
    main(args)
