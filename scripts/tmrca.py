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
except:
    print('''
    Usage: python tmrca.py PATH_TO_ARGON_JARFILE OUTPATH N_E SAMPLE_SIZE LEN_IN_BPS
    ''')
    sys.exit()

size = base_pairs / 1e6
tmrca_file = os.path.join(outpath, 'tmrca.txt')

rho = 1e-8
mu = 1e-7

assert N > sample_size

prefix = outpath + '/'
prefix += 'N' + str(N)
prefix += 'rec' + str(rho)
prefix += 'mut' + str(mu)
treefile = prefix + '.trees'

command = 'java -jar ' + jarfile
command += ' -N ' + str(N)
command += ' -pop 1 ' + str(sample_size)
command += ' -rec ' + str(rho)
command += ' -mut ' + str(mu)
command += ' -size ' + str(size)
command += ' -trees -age'
command += ' -out ' + prefix

ret = subprocess.run(command, shell=True, check=True)

with open(treefile, 'r') as f:
    with open(tmrca_file, 'w') as g:
        for line in f:
            tree = ete3.Tree(line.split()[3], format=1)
            root = tree.get_tree_root()
            tmrca = root.name.split('_')[0]

            g.write(tmrca + '\n')
