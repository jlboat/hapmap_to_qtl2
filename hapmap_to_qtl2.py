import sys
import argparse
import pandas as pd

def parse_arguments():
    """Parse arguments passed to script"""
    parser = argparse.ArgumentParser(description="This script was designed to convert " +
                                                 "a hapmap file into R qtl2 format.\n\n")

    required_named = parser.add_argument_group('required arguments')

    required_named.add_argument(
        "--hapmap",
        type=str,
        required=True,
        help="The name of the input hapmap file. ",
        action="store")

    parser.add_argument(
        "--output",
        type=str,
        required=False,
        default="Rqtl2.genotypes.csv",
        help="File name for qtl2 output.",
        action="store")

    return parser.parse_args()


def write_qtl2(header, samples, variants, output):
    sys.stderr.write("Number of samples: " + 
                     str(len(variants[0])) + "\n")
    sys.stderr.write("Number of variants: " + 
                     str(len(variants)) + "\n")
    with open(output, 'w') as f:
        f.write(",".join(header) + "\n")
        for x, row in enumerate(pd.DataFrame(variants).transpose().values):
            f.write(samples[x] + "," + ",".join(row) + "\n")
 

def hapmap_to_qtl2(hapmap, output):
    """Given a hapmap file, convert genotypes into qtl2 format"""
    # ref -> A
    # alt -> B
    # N -> -
    variants_header = ["id"]
    all_variants = []

    for x, line in enumerate(hapmap):
        split_line = line.split("\t")
        if x == 0:
            samples = split_line[11:]
        else:
            variants_header.append(split_line[0])
            split_alleles = split_line[1].split('/')
            ref_allele = split_alleles[0]
            try:
                alt_allele = split_alleles[1]
            except IndexError:
                sys.stderr.write("ERROR: Monomorphic sites must be removed.\n")
                sys.exit(1)
            variant_list = []
            for y, variant in enumerate(split_line[11:]):
                if variant == ref_allele:
                    variant_list.append("A")
                elif variant == alt_allele:
                    variant_list.append("B")
                else:
                    variant_list.append("-")
            all_variants.append(variant_list)
    write_qtl2(variants_header, samples, 
                   all_variants, output)


if __name__ == '__main__':
    arguments = parse_arguments()
    with open(arguments.hapmap) as f:
        hapmap = f.read().splitlines()
    hapmap_to_qtl2(hapmap, arguments.output)

