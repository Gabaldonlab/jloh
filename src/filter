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


# help section
if len(sys.argv) == 1:
    sys.argv.append("--help")

if (sys.argv[1] in ["--help", "-h", "-help", "help", "getopt", "usage"]):
    sys.stderr.write("""

Filter JLOH results according to different criteria

Usage:
jloh filter --loh <LOH.tsv> [options]

[I/O/E]
--loh               TSV file produced by "jloh extract"                         [!]
--bed               Emit intervals also in BED file with this name              [off]

[basic filters]
--length            Min. LOH length                                             [off]
--alleles           Choices: REF, ALT, NA, multiple choices sep by SPACE        [off]
--zygosity          Choices: homo, hemi, NA, multiple choices sep by SPACE      [off]
--coverage          Min. coverage to retain a block, e.g. "10"                  [off]
--region            Chr:start-end, or simply Chr                                [off]

[by SNPs]
--snps              Min. total number of SNPs in block                          [off]
--het-snps          Min. total number of heterozygous SNPs in block             [off]
--homo-snps         Min. total number of homozygous SNPs in block               [off]
--snps-kbp          Min. SNPs/kbp in block                                      [off]
--het-snps-kbp      Min. heterozygous SNPs/kbp in block                         [off]
--homo-snps-kbp     Min. homozygous SNPs/kbp in block                           [off]


""")
    sys.exit(0)


p = ap.ArgumentParser()
# input and output
p.add_argument("--loh", required=True)
p.add_argument("--bed", type=str)
# basic filters
p.add_argument("--length", default=0, type=int)
p.add_argument("--alleles", nargs="*", default="ALL")
p.add_argument("--zygosity", nargs="*", default="ALL")
p.add_argument("--coverage", default=0, type=int)
p.add_argument("--region", default="ALL")
# by SNPs
p.add_argument("--snps", default=0, type=int)
p.add_argument("--het-snps", default=0, type=int)
p.add_argument("--homo-snps", default=0, type=int)
# by density
p.add_argument("--snps-kbp", default=0, type=float)
p.add_argument("--het-snps-kbp", default=0, type=float)
p.add_argument("--homo-snps-kbp", default=0, type=float)
args = p.parse_args()

# functions

def filter_by_length(args, df):

    """
    07/03/2022
    Remove blocks shorter than a certain length
    """

    df = df.loc[ df["Length"] >= args.length , : ]

    return df


def filter_by_allele(args, df):

    """
    07/03/2022
    Keep only blocks of a certain allele
    """

    if "ALL" in args.alleles:
        return df

    else:
        df = df.loc[ df["Allele"].isin(args.alleles), : ]
        return df


def filter_by_zygosity(args, df):

    """
    07/03/2022
    Keep only blocks of a certain zygosity
    """

    if "ALL" in args.zygosity:
        return df

    else:
        df = df.loc[ df["Zygosity"].isin(args.zygosity), : ]
        return df


def filter_by_coverage(args, df):

    """
    07/03/2022
    Keep only blocks with at least a certain coverage
    """

    df.loc[:,"Cov"] = df.loc[:,"Cov"].str.replace("x", "")
    df.loc[:,"Cov"] = pd.to_numeric(df.loc[:,"Cov"])
    df = df.loc[ df["Cov"] >= args.coverage, : ]

    return df


def filter_by_region(args, df):

    """
    07/03/2022
    Keep only blocks found within a certain region
    """

    if args.region != "ALL":

        Region = args.region.split(":")

        if len(Region) == 1:
            chr = str(Region[0])
            start = 0
            end = float("inf")

        elif len(Region) > 1:
            chr = str(Region[0])
            start, end = Region[1].split("-")
            start, end = int(start), int(end)

        else:
            sys.stderr.write(f"ERROR: region {args.region} cannot be parsed. Should be formatted as:\n")
            sys.stderr.write("Chromosome:Start-End, e.g. Chr_XI:12142-1251985\n")
            sys.stderr.write("Or simply with a chromosome name, e.g. Chr_XI\n")
            sys.stderr.write("Exiting... \n")
            sys.exit()

        mask_1 = df["#Chrom"] == chr
        mask_2 = df["Start"] >= start
        mask_3 = df["End"] <= end

        df = df.loc[ mask_1 & mask_2 & mask_3, : ]

        return df

    else:
        return df


def filter_by_snp_count(args, df):

    """
    07/03/2022
    Apply SNP filters based on raw SNP counts, e.g. num of heterozygous SNPs
    """

    # min snps
    df = df.loc[ df["Homo_snps"] + df["Hetero_snps"] >= args.snps , : ]

    # min het snps
    df = df.loc[ df["Hetero_snps"] >= args.het_snps , : ]

    # min homo snps
    df = df.loc[ df["Homo_snps"] >= args.homo_snps , : ]

    return df


def filter_by_snp_density(args, df):

    """
    07/03/2022
    Apply SNP filters based on SNP density, e.g. het SNPs/kbp
    """

    # min snps/kbp
    df = df.loc[ (df["Homo_snps"] + df["Hetero_snps"]) / (df["Length"] / 1000) >= args.snps_kbp , : ]

    # min het snps/kbp
    df = df.loc[ df["Hetero_snps"] / (df["Length"] / 1000) >= args.het_snps_kbp , : ]

    # min homo snps/kbp
    df = df.loc[ df["Homo_snps"] / (df["Length"] / 1000) >= args.homo_snps_kbp , : ]

    return df


def apply_basic_filters(args, df):

    """
    07/03/2022
    Filter blocks according to general-purpose criteria
    """

    df = filter_by_length(args, df)
    df = filter_by_allele(args, df)
    df = filter_by_zygosity(args, df)
    df = filter_by_coverage(args, df)
    df = filter_by_region(args, df)

    return df


def apply_snp_filters(args, df):

    """
    07/03/2022
    Filter blocks according to criteria based on SNPs and their density
    """

    df = filter_by_snp_count(args, df)
    df = filter_by_snp_density(args, df)

    return df


def main(args):

    """
    07/03/2022
    Main function
    """

    df = pd.read_csv(args.loh, comment="@", sep="\t")
    df = df.fillna("NA")
    df = apply_basic_filters(args, df)
    df = apply_snp_filters(args, df)

    df.to_csv(sys.stdout, sep="\t", header=True, index=False)

    if args.bed:
        df_bed = df.loc[:,["#Chrom", "Start", "End"]].copy()
        df_bed["Start"] = df_bed["Start"].apply(lambda x : max(0, x-1))
        df_bed.to_csv(args.bed, sep="\t", header=False, index=False)


if __name__ == "__main__":
    main(args)
