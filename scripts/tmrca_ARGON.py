import numpy as np
import sys, os
import subprocess
import gzip

import ete3

try:
    jarfile = os.path.expanduser(sys.argv[1])
    outpath = os.path.expanduser(sys.argv[2])
    N = int(sys.argv[3])
    sample_size = int(sys.argv[4])
    base_pairs = float(sys.argv[5])
    n_iter = int(sys.argv[6])
except:
    print('''
    Usage:
    python tmrca.py PATH_TO_ARGON_JARFILE OUTPATH N_E SAMPLE_SIZE LEN_IN_BPS N_ITER
    ''')
    sys.exit()

size = base_pairs / 1e6

rho = 1e-8
mu = 1e-8

assert N >= sample_size

prefix = outpath + '/'
prefix += 'N' + str(N)
prefix += 'samplesize' + str(sample_size)
prefix += 'rec' + str(rho)
prefix += 'mut' + str(mu)
prefix += 'length' + str(base_pairs)
treefile = prefix + '.trees'
tmrca_file = prefix + '_tmrca.txt'

command = 'java -jar ' + jarfile
command += ' -N ' + str(N)
command += ' -pop 1 ' + str(sample_size)
command += ' -rec ' + str(rho)
command += ' -mut ' + str(mu)
command += ' -size ' + str(size)
command += ' -trees -age'
command += ' -out ' + prefix

for i in range(n_iter):
    ret = subprocess.run(command, shell=True, check=True)

    with open(treefile, 'r') as f:
        with open(tmrca_file, 'a') as g:
            for line in f:
                tree = ete3.Tree(line.split()[3], format=1)
                root = tree.get_tree_root()
                tmrca = root.name.split('_')[0]

                g.write(tmrca + '\n')

cleanup = 'rm ' + prefix + '.*'
subprocess.run(cleanup, shell=True, check=True)
