import sys, os
sys.path.append(os.path.expanduser('~/project/msprime'))
import msprime

ts_file = os.path.expanduser(sys.argv[1])
vcf_file = os.path.expanduser(sys.argv[2])

ts = msprime.load(ts_file)
print(len(list(ts.trees())), "trees")

print("Writing to", vcf_file)
with open(vcf_file, 'w') as f:
    ts.write_vcf(f, ploidy=2)
