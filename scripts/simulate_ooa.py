import sys, os
import math
from datetime import datetime

# Kind of a dirty hack but it works for now
msprime_dir = os.path.expanduser('../../msprime')
sys.path.insert(1, msprime_dir)
sys.path.insert(1, msprime_dir + '/lib/subprojects/git-submodules/tskit/python')
import msprime


def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


def current_time():
    return(' [' + datetime.strftime(datetime.now(), '%Y-%m-%d %H:%M:%S') + ']')


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


def get_whole_genome_positions_rates_loci(num_chroms=22, rho=1e-8):
    all_lengths_morgans = [2.77693825, 2.633496065, 2.24483368, 2.12778391, 
            2.03765845, 1.929517394, 1.867959329, 1.701765192, 1.68073935, 
            1.789473882, 1.594854258, 1.72777271, 1.26940475, 1.16331251, 
            1.2554709, 1.348911043, 1.29292106, 1.18978483, 1.077960694, 
            1.079243479, 0.61526812, 0.72706815]
    all_lengths = [x * 1e8 for x in all_lengths_morgans]

    chrom_lengths = all_lengths[:num_chroms]

    positions, rates = get_positions_rates(chrom_lengths, rho)
    num_loci = positions[-1]

    ## HACK: This is to avoid bad edge intervals (right <= left) when
    ## simulating whole genomes. Possible overflow / floating point precision
    ## issue
    num_loci = int(num_loci / 100)
    eprint("Scaling to", num_loci, "loci")

    return positions, rates, num_loci


def out_of_africa_alicia(nhaps):
    """
    Specify the demographic model used in these simulations (Gravel et al, 2011 PNAS)
    """
    # First we set out the maximum likelihood values of the various parameters
    # given in Gravel et al, 2011 Table 2.
    N_A = 7300
    N_B = 1861
    N_AF = 14474
    N_EU0 = 1032
    N_AS0 = 554
    # Times are provided in years, so we convert into generations.
    generation_time = 25
    T_AF = 148e3 / generation_time
    T_B = 51e3 / generation_time
    T_EU_AS = 23e3 / generation_time
    # We need to work out the starting (diploid) population sizes based on
    # the growth rates provided for these two populations
    r_EU = 0.0038
    r_AS = 0.0048
    N_EU = N_EU0 / math.exp(-r_EU * T_EU_AS)
    N_AS = N_AS0 / math.exp(-r_AS * T_EU_AS)
    # Migration rates during the various epochs.
    m_AF_B = 15e-5
    m_AF_EU = 2.5e-5
    m_AF_AS = 0.78e-5
    m_EU_AS = 3.11e-5
    
    # Population IDs correspond to their indexes in the population
    # configuration array. Therefore, we have 0=YRI, 1=CEU and 2=CHB
    # initially.
    population_configurations = [
        msprime.PopulationConfiguration(
            sample_size=nhaps[0], initial_size=N_AF),
        msprime.PopulationConfiguration(
            sample_size=nhaps[1], initial_size=N_EU, growth_rate=r_EU),
        msprime.PopulationConfiguration(
            sample_size=nhaps[2], initial_size=N_AS, growth_rate=r_AS)
    ]
    migration_matrix = [
        [      0, m_AF_EU, m_AF_AS],
        [m_AF_EU,       0, m_EU_AS],
        [m_AF_AS, m_EU_AS,       0],
    ]
    demographic_events = [
        # CEU and CHB merge into B with rate changes at T_EU_AS
        msprime.MassMigration(
            time=T_EU_AS, source=2, destination=1, proportion=1.0),
        msprime.MigrationRateChange(time=T_EU_AS, rate=0),
        msprime.MigrationRateChange(
            time=T_EU_AS, rate=m_AF_B, matrix_index=(0, 1)),
        msprime.MigrationRateChange(
            time=T_EU_AS, rate=m_AF_B, matrix_index=(1, 0)),
        msprime.PopulationParametersChange(
            time=T_EU_AS, initial_size=N_B, growth_rate=0, population_id=1),
        # Population B merges into YRI at T_B
        msprime.MassMigration(
            time=T_B, source=1, destination=0, proportion=1.0),
        # Size changes to N_A at T_AF
        msprime.PopulationParametersChange(
            time=T_AF, initial_size=N_A, population_id=0)
    ]
    # Use the demography debugger to print out the demographic history
    # that we have just described.
    dp = msprime.DemographyDebugger(
        Ne=N_A,
        population_configurations=population_configurations,
        migration_matrix=migration_matrix,
        demographic_events=demographic_events)
    dp.print_history()
    
    return(population_configurations, migration_matrix, demographic_events)

def out_of_africa_msprime(nhaps):
    """
    Specify the demographic model used in these simulations (msprime.readthedocs.io)
    """
    # First we set out the maximum likelihood values of the various parameters
    # given in Table 1.
    N_A = 7300
    N_B = 2100
    N_AF = 12300
    N_EU0 = 1000
    N_AS0 = 510
    # Times are provided in years, so we convert into generations.
    generation_time = 25
    T_AF = 220e3 / generation_time
    T_B = 140e3 / generation_time
    T_EU_AS = 21.2e3 / generation_time
    # We need to work out the starting (diploid) population sizes based on
    # the growth rates provided for these two populations
    r_EU = 0.004
    r_AS = 0.0055
    N_EU = N_EU0 / math.exp(-r_EU * T_EU_AS)
    N_AS = N_AS0 / math.exp(-r_AS * T_EU_AS)
    # Migration rates during the various epochs.
    m_AF_B = 25e-5
    m_AF_EU = 3e-5
    m_AF_AS = 1.9e-5
    m_EU_AS = 9.6e-5
    # Population IDs correspond to their indexes in the population
    # configuration array. Therefore, we have 0=YRI, 1=CEU and 2=CHB
    # initially.
    population_configurations = [
        msprime.PopulationConfiguration(
            sample_size=nhaps[0], initial_size=N_AF),
        msprime.PopulationConfiguration(
            sample_size=nhaps[1], initial_size=N_EU, growth_rate=r_EU),
        msprime.PopulationConfiguration(
            sample_size=nhaps[2], initial_size=N_AS, growth_rate=r_AS)
    ]
    migration_matrix = [
        [      0, m_AF_EU, m_AF_AS],
        [m_AF_EU,       0, m_EU_AS],
        [m_AF_AS, m_EU_AS,       0],
    ]
    demographic_events = [
        # CEU and CHB merge into B with rate changes at T_EU_AS
        msprime.MassMigration(
            time=T_EU_AS, source=2, destination=1, proportion=1.0),
        msprime.MigrationRateChange(time=T_EU_AS, rate=0),
        msprime.MigrationRateChange(
            time=T_EU_AS, rate=m_AF_B, matrix_index=(0, 1)),
        msprime.MigrationRateChange(
            time=T_EU_AS, rate=m_AF_B, matrix_index=(1, 0)),
        msprime.PopulationParametersChange(
            time=T_EU_AS, initial_size=N_B, growth_rate=0, population_id=1),
        # Population B merges into YRI at T_B
        msprime.MassMigration(
            time=T_B, source=1, destination=0, proportion=1.0),
        # Size changes to N_A at T_AF
        msprime.PopulationParametersChange(
            time=T_AF, initial_size=N_A, population_id=0)
    ]
    # Use the demography debugger to print out the demographic history
    # that we have just described.
    dd = msprime.DemographyDebugger(
        population_configurations=population_configurations,
        migration_matrix=migration_matrix,
        demographic_events=demographic_events)
    dd.print_history()

    return(population_configurations, migration_matrix, demographic_events)


def simulate_ooa(population_configurations, migration_matrix, demographic_events, recomb_map, models):
    """
    simulate according to the specified demographic model with recombination
    """
    ret = []
    for model in models:
        eprint('Starting simulations (' + model + ')' + current_time())
        simulation = msprime.simulate(
            population_configurations = population_configurations,
            migration_matrix=migration_matrix,
            demographic_events=demographic_events,
            model=model,
            mutation_rate=2e-8,
            recombination_map=recomb_map
        )
        eprint('Ending simulations' + current_time())

        ret.append(simulation)

    return ret


if __name__ == "__main__":
    # num_chroms = 1
    num_chroms = 22
    # nhaps = [1, 1, 0]
    nhaps = [1000, 1000, 1000]
    positions, rates, num_loci = get_whole_genome_positions_rates_loci(
            num_chroms=num_chroms)
    recomb_map = msprime.RecombinationMap(
            positions, rates, num_loci=num_loci)
    # models = ['dtwf']
    models = ['hudson']
    # models = ['dtwf', 'hudson']

    (pop_config, mig_mat, demog) = out_of_africa_msprime(nhaps)
    simulations = simulate_ooa(pop_config, mig_mat, demog, recomb_map, models)

    for sim, model in zip(simulations, models):
        sim.dump(os.path.expanduser('./ooa_') + model + \
                        '_nhaps_' + '_'.join(map(str, nhaps)) + '.hdf5', True
                        )


