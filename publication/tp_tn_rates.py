#!/usr/bin/env python3

# modules
import sys
import argparse as ap
import pybedtools
from pybedtools import BedTool

# parser
p = ap.ArgumentParser()
p.add_argument("--generated", required=True)
p.add_argument("--predicted", required=True)
p.add_argument("--genome-file", required=True)
p.add_argument("--ovl-rate", default=0.50, type=float)
p.add_argument("--out")
args = p.parse_args()

# functions
def get_statistics(generated, predicted, ovl_rate, loh_type, genome_file):

    # read files
    bt_gen = BedTool(generated)
    bt_pred = BedTool(predicted)

    # get total positions
    INPUT = open(args.genome_file, "r")
    Lines = [line.rstrip("\b\r\n").split("\t") for line in INPUT]
    INPUT.close()
    tot_pos = sum([int(x[1]) for x in Lines])

    # get all positives (tp + fn)
    Lines = [str(line).rstrip("\b\r\n").split("\t") for line in bt_gen]
    all_true = sum([int(x[2])-int(x[1]) for x in Lines])

    # get true positives
    # a.k.a. positions that ARE predicted inside LOH blocks and DO belong to one
    bt_tp = bt_pred.intersect(b=bt_gen, wo=True)
    tp = [str(line).rstrip("\b\r\n").split("\t") for line in bt_tp]
    tp = sum([int(x[6]) for x in bt_tp])

    # get false negatives
    # a.k.a. positions that ARE NOT predicted inside LOH blocks but DO belong to one
    bt_pred_comp = bt_pred.complement(g=genome_file)
    bt_fn = bt_pred_comp.intersect(b=bt_gen, wo=True)
    fn = [str(line).rstrip("\b\r\n").split("\t") for line in bt_fn]
    fn = sum([int(x[6]) for x in bt_fn])

    # get false positives
    # a.k.a. positions that ARE predicted inside LOH blocks but DO NOT belong to one
    bt_gen_comp = bt_gen.complement(g=genome_file)
    bt_fp = bt_pred.intersect(b=bt_gen_comp, wo=True)
    fp = [str(line).rstrip("\b\r\n").split("\t") for line in bt_fp]
    fp = sum([int(x[6]) for x in bt_fp])

    # get true negatives
    # a.k.a. positions that ARE NOT predicted inside LOH blocks and DO NOT belong to one
    bt_pred_comp = bt_pred.complement(g=genome_file)
    bt_gen_comp = bt_gen.complement(g=genome_file)
    bt_tn = bt_pred_comp.intersect(b=bt_gen_comp, wo=True)
    tn = [str(line).rstrip("\b\r\n").split("\t") for line in bt_tn]
    tn = sum([int(x[6]) for x in bt_tn])

    # form statistics
    try:
        precision = round(tp / (tp+fp), 3)
    except ZeroDivisionError:
        precision = 0

    try:
        recall = round(tp / (tp+fn), 3)
        sensitivity = recall
    except ZeroDivisionError:
        recall = 0
        sensitivity = 0

    try:
        specificity = round(tn / (tn+fp), 3)
    except ZeroDivisionError:
        specificity = 0

    return (loh_type, tp, fp, fn, tn, precision, recall, sensitivity, specificity)


def main(args):

    # predicted blocks
    OUT = open(args.out, "w")
    out_pred = [str(x) for x in get_statistics(args.generated, args.predicted, args.ovl_rate, "predicted", args.genome_file)]
    get_statistics(args.generated, args.predicted, args.ovl_rate, "predicted", args.genome_file)
    OUT.write("TYPE\tTP\tFP\tFN\tTN\tPRECISION\tRECALL\tSENSITIVITY\tSPECIFICITY\n")
    OUT.write("\t".join(out_pred) + "\n")
    OUT.close()


if __name__ == "__main__":
    main(args)
