import sys, os
import argparse


def main(args):
    
    chrom_lengths = [
            247249719, 242951149, 199501827, 191273063, 180857866,
            170899992, 158821424, 146274826, 140273252, 135374737,
            134452384, 132349534, 114142980, 106368585, 100338915,
            88827254, 78774742, 76117153, 63811651, 62435964,
            46944323, 49691432
            ]

    out_file = os.path.expanduser(args.vcf_out)
    in_file = os.path.expanduser(args.vcf_in)

    shift = 0
    chrom_idx = 0
    with open(out_file, 'w') as g:
        with open(in_file, 'r') as f:
            for line in f.readlines():
                if line[0] == '#':
                    g.write(line)
                    continue

                split = line.split('\t')
                chrom = int(split[0])
                pos = int(split[1])

                if pos > chrom_lengths[chrom_idx] + shift:
                    shift += chrom_lengths[chrom_idx]
                    chrom_idx += 1

                shifted_chrom = chrom_idx + 1
                shifted_pos = pos - shift
                assert shifted_pos > 0

                split[0] = str(shifted_chrom)
                split[1] = str(shifted_pos)
                shifted_line = '\t'.join(split)
                g.write(shifted_line)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--vcf_in",
            metavar='|',
            required=True,
            help="VCF from ts with chromosomes defined in" +\
                    " recombination map, but all set to 1")
    parser.add_argument("--vcf_out",
            metavar='|',
            required=True,
            help="Path to write VCF with split chromosomes")

    args = parser.parse_args()
    main(args)
