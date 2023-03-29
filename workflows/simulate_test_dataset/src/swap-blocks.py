#!/usr/bin/env python3 

import argparse as ap 
import sys 
from Bio import SeqIO
import pandas as pd 
import numpy as np 
import random 
from time import asctime as at 

ss = sys.exit 

# parser 
p = ap.ArgumentParser()
p.add_argument("--genome-A", required=True, help="Input FASTA genome A")
p.add_argument("--genome-B", required=True, help="Input FASTA genome B")
p.add_argument("--out-A", required=True, help="Output FASTA genome A")
p.add_argument("--out-B", required=True, help="Output FASTA genome A")
p.add_argument("--loh-fraction", type=float, default=0.20, help="Fraction of output genomes to undergo LOH")
p.add_argument("--mean-haplotype-size", default=25000, type=int, help="Minimum size of the generated haplotypes")
p.add_argument("--min-haplotype-length", default=1000, type=int, help="Minimum length of the generated haplotypes")
p.add_argument("--towards", type=str, choices=["A", "B"], default="B", help="Which subgenome is taking over")
p.add_argument("--weight", type=float, default=0.75, help="Extent of takeover (0.5-1.0)")
args = p.parse_args()

# functions 
def get_chrom_sizes(fasta):

    Chrom_sizes = {}
    for record in SeqIO.parse(fasta, "fasta"):
        Chrom_sizes[str(record.id)] = int(len(record.seq))

    return Chrom_sizes


def get_chrom_sequences(fasta):

    Chroms = {}
    for record in SeqIO.parse(fasta, "fasta"):
        Chroms[str(record.id)] = str(record.seq)

    return Chroms


def read_genomes(genome_A, genome_B):

    Chrom_sizes_A = get_chrom_sizes(genome_A)
    Chrom_sizes_B = get_chrom_sizes(genome_B)

    Chrom_seqs_A = get_chrom_sequences(genome_A)
    Chrom_seqs_B = get_chrom_sequences(genome_B)

    INPUT = open(genome_A, "r")
    Chrom_order_A = [line.lstrip(">").rstrip("\b\r\n") for line in INPUT if line[0] == ">"]
    INPUT.close()

    INPUT = open(genome_B, "r")
    Chrom_order_B = [line.lstrip(">").rstrip("\b\r\n") for line in INPUT if line[0] == ">"]
    INPUT.close()

    Data = (
        Chrom_sizes_A, Chrom_seqs_A, Chrom_order_A, 
        Chrom_sizes_B, Chrom_seqs_B, Chrom_order_B
        )

    return Data


def get_chrom_random_intervals(mean_haplotype_size, min_haplotype_len, stdev_factor, chrom_len):

    # initialize variables and containers 
    Intervals = []
    start, end = 0, 0

    # iterate over chrom length 
    while (end < chrom_len):
        start = int(end)
        increment = np.random.normal( loc=float(mean_haplotype_size), 
                                      scale=float(mean_haplotype_size)*stdev_factor)
        
        # if increment is smaller than min_haplotype_len, take min_haplotype_len
        increment = max(min_haplotype_len, increment)

        # if end is above chrom length, pick chrom length 
        end = int(min(chrom_len, end+increment))

        # if end is near chrom length, pick chrom length 
        if (chrom_len - end < min_haplotype_len * 1.01):
            end = chrom_len

        # append interval 
        Intervals.append((start, end))

    return Intervals


def break_into_haplotypes(Data, mean_haplotype_size, min_haplotype_length):

    stdev_factor=0.25

    (
        Chrom_sizes_A, Chrom_seqs_A, Chrom_order_A, 
        Chrom_sizes_B, Chrom_seqs_B, Chrom_order_B
    ) = Data
    
    Chrom_intervals_A = { chrom : get_chrom_random_intervals(mean_haplotype_size, min_haplotype_length, stdev_factor, Chrom_sizes_A[chrom]) for chrom in Chrom_sizes_A.keys() }
    Chrom_intervals_B = { chrom : get_chrom_random_intervals(mean_haplotype_size, min_haplotype_length, stdev_factor, Chrom_sizes_B[chrom]) for chrom in Chrom_sizes_B.keys() }

    Data = (    Chrom_sizes_A, Chrom_seqs_A, Chrom_order_A, Chrom_intervals_A,
                Chrom_sizes_B, Chrom_seqs_B, Chrom_order_B, Chrom_intervals_B )

    return Data


def simulate_loh_haplotypes(Data, loh_fraction, towards, weight):
    
    (   Chrom_sizes_A, Chrom_seqs_A, Chrom_order_A, Chrom_intervals_A,
        Chrom_sizes_B, Chrom_seqs_B, Chrom_order_B, Chrom_intervals_B ) = Data

    Loh_intervals_A = {}
    Loh_intervals_B = {}

    # genome A
    for chrom in Chrom_intervals_A.keys():

        num_items = len(Chrom_intervals_A[chrom])
        loh_items = int(round(num_items * loh_fraction, 0))

        if loh_items > 0:
            idx = random.choices(range(0, num_items), k=loh_items)
            Ints = [Chrom_intervals_A[chrom][i] for i in idx]
            Loh_intervals_A[chrom] = Ints 
        else:
            Loh_intervals_A[chrom] = []
    
    # genome B 
    for chrom in Chrom_intervals_B.keys():

        num_items = len(Chrom_intervals_B[chrom])
        loh_items = int(round(num_items * loh_fraction, 0))

        if loh_items > 0:
            idx = random.choices(range(0, num_items), k=loh_items)
            Ints = [Chrom_intervals_B[chrom][i] for i in idx]
            Loh_intervals_B[chrom] = Ints 
        else:
            Loh_intervals_B[chrom] = []

    Data = (
        Chrom_sizes_A, Chrom_seqs_A, Chrom_order_A, Chrom_intervals_A, Loh_intervals_A, 
        Chrom_sizes_B, Chrom_seqs_B, Chrom_order_B, Chrom_intervals_B, Loh_intervals_B
    )
    
    return Data 


def swap_blocks(Data):
    
    (
        Chrom_sizes_A, Chrom_seqs_A, Chrom_order_A, Chrom_intervals_A, Loh_intervals_A, 
        Chrom_sizes_B, Chrom_seqs_B, Chrom_order_B, Chrom_intervals_B, Loh_intervals_B
    ) = Data

    New_chrom_seqs_A = {}
    New_chrom_seqs_B = {}

    Matching_chrom_A = { Chrom_order_A[i] : Chrom_order_B[i] for i in range(0, len(Chrom_order_A))}
    Matching_chrom_B = { Chrom_order_B[i] : Chrom_order_A[i] for i in range(0, len(Chrom_order_B))}

    # genome A
    for chrom in Loh_intervals_A.keys():
        newseq = Chrom_seqs_A[chrom]

        for interval in Loh_intervals_A[chrom]:
            start = interval[0]
            end = interval[1]
            match = Matching_chrom_A[chrom]
            newseq = newseq[:start] + Chrom_seqs_B[match][start:end] + newseq[end:]
            New_chrom_seqs_A[chrom] = newseq 

    # genome B
    for chrom in Loh_intervals_B.keys():
        newseq = Chrom_seqs_B[chrom]

        for interval in Loh_intervals_B[chrom]:
            start = interval[0]
            end = interval[1]
            match = Matching_chrom_B[chrom]
            newseq = newseq[:start] + Chrom_seqs_A[match][start:end] + newseq[end:]
            New_chrom_seqs_B[chrom] = newseq 

    Data = (
        Chrom_sizes_A, New_chrom_seqs_A, Chrom_order_A, Chrom_intervals_A, Loh_intervals_A, 
        Chrom_sizes_B, New_chrom_seqs_B, Chrom_order_B, Chrom_intervals_B, Loh_intervals_B
    )

    return Data 


def write_to_output(Data, Outs):

    (
        Chrom_sizes_A, New_chrom_seqs_A, Chrom_order_A, Chrom_intervals_A, Loh_intervals_A, 
        Chrom_sizes_B, New_chrom_seqs_B, Chrom_order_B, Chrom_intervals_B, Loh_intervals_B
    ) = Data

    # sequences
    out_A = open(Outs[0], "w")
    out_B = open(Outs[1], "w")

    for chrom in Chrom_order_A:
        seq = New_chrom_seqs_A[chrom]
        out_A.write(f">{chrom}\n{seq}\n")

    for chrom in Chrom_order_B:
        seq = New_chrom_seqs_B[chrom]
        out_B.write(f">{chrom}\n{seq}\n")

    out_A.close()
    out_B.close()


    # loh haplotypes
    out_A = open(f"{Outs[0]}.true_loh_blocks.bed", "w")
    out_B = open(f"{Outs[1]}.true_loh_blocks.bed", "w")

    for chrom in Loh_intervals_A:
        for x in Loh_intervals_A[chrom]:
            start = x[0]
            end = x[1]
            out_A.write(f"{chrom}\t{start}\t{end}\n")

    for chrom in Loh_intervals_B:
        for x in Loh_intervals_B[chrom]:
            start = x[0]
            end = x[1]
            out_B.write(f"{chrom}\t{start}\t{end}\n")

    out_A.close()
    out_B.close()


### MAIN ###

def main(args):
    
    # read genomes 
    sys.stderr.write(f"[{at()}] Reading input genomes ... ")
    Data = read_genomes(args.genome_A, args.genome_B)
    sys.stderr.write(f"DONE\n")

    # break into haplotypes 
    sys.stderr.write(f"[{at()}] Breaking them into haplotypes ... ")
    Data = break_into_haplotypes(Data, args.mean_haplotype_size, args.min_haplotype_length)
    sys.stderr.write(f"DONE\n")

    # simulate loh based on --loh-fraction --towards and --weight
    sys.stderr.write(f"[{at()}] Simulating LOH blocks between genomes ... ")
    Data = simulate_loh_haplotypes(Data, args.loh_fraction, args.towards, args.weight)
    sys.stderr.write(f"DONE\n")
    
    # swap blocks based on assigned identities
    sys.stderr.write(f"[{at()}] Swapping blocks based on assignment ... ")
    Data = swap_blocks(Data)
    sys.stderr.write(f"DONE\n")

    # write to output 
    sys.stderr.write(f"[{at()}] Writing to output ... ")
    write_to_output(Data, [args.out_A, args.out_B])
    sys.stderr.write(f"DONE\n")


if __name__ == "__main__":
    main(args)
    