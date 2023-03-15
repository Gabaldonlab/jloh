#!/usr/bin/env python3

"""
###
JLOH - Inferring Loss of Heterozygosity Blocks from Short-read sequencing data

Copyright (C) 2023 Matteo Schiavinato

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
###
"""


import argparse as ap
import sys
import pandas as pd
import pybedtools
from pybedtools import BedTool


# help section
if len(sys.argv) == 1:
    sys.argv.append("--help")

if (sys.argv[1] in ["--help", "-h", "-help", "help", "getopt", "usage"]):
    sys.stderr.write("""

Perform intersection/removal operations with output files

Usage:
jloh intersect --loh-A <LOH_A.tsv> --loh-B <LOH_B.tsv> [options]

[I/O/E]
--loh-A             TSV file produced by "jloh extract"                         [!]
--loh-B             . . . . . . . . . . . . . . . . . .                         [!]

[operations]
--mode              Type of operation to carry on, among:                       [intersection]
                    - "intersection" (default): keep lines in common
                    - "complement": keep lines unique to --loh-A
                                  (may be subject to --min-ovl)
                    - "unique": exclude common lines
--min-ovl           Fraction of block overlap to trigger retention/removal      [0.01]

""")
    sys.exit(0)


p = ap.ArgumentParser()
# input and output
p.add_argument("--loh-A", required=True)
p.add_argument("--loh-B", required=True)
p.add_argument("--mode", choices=["intersection", "complement", "unique"], default="intersection")
p.add_argument("--min-ovl", default=0.01, type=float)
args = p.parse_args()


# functions
def get_intersection(file_A, file_B, args):

    """
    13/07/2022
    Matteo Schiavinato
    """

    btA = BedTool(file_A)
    btB = BedTool(file_B)

    bt_out = btA.intersect(b=btB, header=True, f=args.min_ovl)
    return bt_out


def get_complement(file_A, file_B, args):

    """
    13/07/2022
    Matteo Schiavinato
    """

    btA = BedTool(file_A)
    btB = BedTool(file_B)

    bt_out = btA.intersect(b=btB, v=True, header=True, f=args.min_ovl)
    return bt_out


def get_unique(file_A, file_B, args):

    """
    13/07/2022
    Matteo Schiavinato
    """

    btA = BedTool(file_A)
    btB = BedTool(file_B)

    bt_out_A = btA.intersect(b=btB, v=True, header=True, f=args.min_ovl)
    bt_out_B = btB.intersect(b=btA, v=True, header=True, f=args.min_ovl)

    bt_out = BedTool(list(bt_out_A) + list(bt_out_B))
    return bt_out



# main code
if __name__ == "__main__":

    ###

    file_A = args.loh_A
    file_B = args.loh_B

    ###

    if args.mode == "intersection":

        """
        Intersection:
        This is the default mode
        If nothing is specified, it runs this one
        """

        Tsv_out = get_intersection(file_A, file_B, args)

    ###

    elif args.mode == "complement":

        """
        Complement:
        This mode removes lines found in both files and leaves only lines found
        only in --loh-A
        The comparison is made only on the first three fields (chrom, start, end)
        The comparison takes into account other parameters, such as --min-ovl
        """

        Tsv_out = get_complement(file_A, file_B, args)

    ###

    elif args.mode == "unique":

        """
        Complement:
        This mode removes common lines between the files, and leaves the ones that
        are unique to either one of them (i.e. not just --loh-A, as in "complement")
        """

        Tsv_out = get_unique(file_A, file_B, args)

    ###

    else:
        sys.stderr.write(f"ERROR: unrecognized --mode: {args.mode}\n")
        sys.exit()

    ###
    
    header="#Chrom\tStart\tEnd\tCov\tCov_frac\tCov_ratio_up\tCov_ratio_down\tZygosity\tLength\tAllele\tHomo_snps\tHetero_snps"
    sys.stdout.write(header + "\n")
    for line in list(Tsv_out):
        sys.stdout.write(str(line))
