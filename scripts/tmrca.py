import numpy as np
import sys, os
import subprocess
import gzip

import ete3

jarfile = os.path.expanduser(sys.argv[1])
outpath = os.path.expanduser(sys.argv[2])
tmrca_file = os.path.join(outpath, 'tmrca.txt')

N = 1000
sample_size = 100
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
