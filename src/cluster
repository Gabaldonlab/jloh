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
import multiprocessing as mp 
import scipy
import scipy.cluster 
from time import asctime as at


ss = sys.exit 


# help section
if len(sys.argv) == 1:
    sys.argv.append("--help")

if (sys.argv[1] in ["--help", "-h", "-help", "help", "getopt", "usage"]):
    sys.stderr.write("""

Cluster different runs by overlap 

Usage:
jloh cluster [options] --loh <LOH_A.tsv> ... <LOH_N.tsv> 

[I/O/E]
--loh               TSV files produced by "jloh extract"                        [!]
--out-prefix        Prefix to use for output files                              [jloh_clust_out]

[parameters]
--max-dist          Maximum distance (float, 0-1) between elements in cluster   [0.1]
--threads           Number of parallel threads                                  [4]

""")
    sys.exit(0)

# parser 
def create_parser():
    p = ap.ArgumentParser()
    p.add_argument("--loh", nargs="+", required=True)
    p.add_argument("--out-prefix", default="jloh_clust_out", type=str)
    p.add_argument("--threads", default=4, type=int)
    p.add_argument("--max-dist", default=0.1, type=float)
    args = p.parse_args()

    return args 

# functions
def dump_queue(q):

    """
    16/03/2022
    """

    out = []
    while not q.empty():
        x = q.get()
        out.append(x)

    return out


def read_input_files(Lohs):

    """
    31/05/2023
    """

    Bts = [(infile, BedTool([line for line in open(infile)])) for infile in Lohs]
    return Bts


def get_jaccard_distance(bt1, infile_1, bt2, infile_2, queue):

    """
    31/05/2023
    """

    jaccard = bt1.jaccard(b=bt2)["jaccard"]
    queue.put((infile_1, infile_2, jaccard))


def get_distances(Bts, args):

    """
    31/05/2023
    """

    # calculate tot distances 
    tot_dists = len(Bts) * len(Bts)

    # create pool and queue for multiprocessing 
    pool = mp.Pool(processes=args.threads)
    queue = mp.Manager().Queue()

    # for each file in the list of BedTool objects 
    # calculate distance from every other file including itself 
    counter = 0
    for i in range(0,len(Bts)):
        infile_1 = Bts[i][0]
        bt1 = Bts[i][1]
        for k in range(0,len(Bts)):
            counter += 1
            infile_2 = Bts[k][0]
            bt2 = Bts[k][1]
            pool.apply_async(get_jaccard_distance, args=(bt1, infile_1, bt2, infile_2, queue,))
            sys.stderr.write(f"[{at()}] Analysing pairwise distance {counter} of {tot_dists}\r")

    pool.close()
    pool.join()
    sys.stderr.write("\n")

    out = dump_queue(queue)
    return out 


def get_distance_matrix(Distances):

    """
    31/05/2023
    """

    Dist_dict = {}
    
    # create dictionary of dictionaries which will be used by 
    # pandas to create a matrix 
    for x in Distances:
        if x[0] in Dist_dict.keys():
            Dist_dict[x[0]][x[1]] = x[2]
        else:
            Dist_dict[x[0]] = {x[1] : x[2]}
    
    # create pandas dataframe
    df = pd.DataFrame(Dist_dict)

    # sort columns as row names so that the matrix is symmetrical
    df = df.loc[: , df.index]

    # invert values to have diagonal 0 
    # jaccard index measures proximity but to cluster we want distance 
    df = abs(df -1) 
    
    return df 


def get_clusters(df_dist, args):

    """
    31/05/2023
    """

    # get condensed matrix 
    df_dist_condensed = scipy.spatial.distance.squareform(df_dist)

    # obtain linkage between objects in matrix
    clust = scipy.cluster.hierarchy.linkage(df_dist_condensed, method='single', metric='euclidean')

    # get belonging of items to clusters 
    clust = scipy.cluster.hierarchy.fcluster(clust, args.max_dist, criterion="distance")

    # subdivide objects by cluster belonging 
    df_clust = pd.DataFrame({"Sample":df_dist.index.tolist(), "Cluster":clust})
    df_clust = df_clust.sort_values(by="Cluster", ascending=True)

    return df_clust 


def main(args):

    """
    31/05/2023
    """

    # read input files into BedTool objects 
    sys.stderr.write(f"[{at()}] Reading input files\n")
    Bts = read_input_files(args.loh)
    sys.stderr.write(f"[{at()}] Read {len(Bts)} files\n")

    # estimate distances using a jaccard index 
    sys.stderr.write(f"[{at()}] Calculating pairwise distances\n")
    Distances = get_distances(Bts, args)
    sys.stderr.write(f"[{at()}] Calculated {len(Distances)} pairwise distances\n")

    # convert into a distance matrix 
    sys.stderr.write(f"[{at()}] Getting a distance matrix\n")
    df_dist = get_distance_matrix(Distances)

    # cluster samples based on distance matrix 
    sys.stderr.write(f"[{at()}] Inferring clusters from distance matrix\n")
    df_clust = get_clusters(df_dist, args)

    # write to output 
    sys.stderr.write(f"[{at()}] Writing to output\n")
    df_dist.round(3).reset_index(drop=False).to_csv(f"{args.out_prefix}.dist.tsv", index=False, header=True, sep="\t")
    df_clust.to_csv(f"{args.out_prefix}.clust.tsv", index=False, header=True, sep="\t")

    sys.stderr.write(f"[{at()}] Done\n")

# main 
if __name__ == "__main__":
    args = create_parser()
    main(args)