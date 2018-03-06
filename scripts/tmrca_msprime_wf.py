import sys, os
sys.path.append(os.path.expanduser('~/projects/msprime'))
import msprime

try:
    outpath = os.path.expanduser(sys.argv[1])
    N = int(sys.argv[2])
    sample_size = int(sys.argv[3])
    length = int(sys.argv[4])
    n_iter = int(sys.argv[5])
except:
    print('''
    Usage:
    python tmrca_msprime_wf.py OUTPATH N_E SAMPLE_SIZE LEN_IN_BP N_ITER
    ''')
    raise

rho = 1e-8
mu = 1e-8

prefix = outpath + '/'
prefix += 'N' + str(N)
prefix += 'samplesize' + str(sample_size)
prefix += 'rec' + str(rho)
prefix += 'mut' + str(mu)
prefix += 'length' + str(length)
tmrca_file = prefix + '_tmrca.txt'

pop_config = msprime.PopulationConfiguration(
        sample_size=sample_size,
        initial_size=N)

for i in range(n_iter):
    ts = msprime.simulate(
            population_configurations=[pop_config],
            recombination_rate=rho,
            mutation_rate=mu,
            model='dtwf')

    ## There could (?) be more than one tree, but then their TMRCAs would
    ## be correlated, which we don't want
    t = next(ts.trees())
    tmrca = t.time(t.get_root())

    with open(tmrca_file, 'a') as g:
        g.write(str(tmrca) + '\n')

