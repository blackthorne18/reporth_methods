import pickle
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import collections
import warnings
from copy import deepcopy
import os

FILEPATH = "./temp/cluster_output_Sep20_875/"
CLUSTERFILENAME = FILEPATH + "/clusters_Sep20.txt"
PATHFILENAME = FILEPATH + "/path_making_Sep20.txt"
OUTPUTFILENAME = "./temp/clusters_extragenic_map.csv"
IGNOREDGENOMES = ['chlPCL1606', "chlO6"]
READMESTR = """
PROGRAM COMPLETE
The output file is a csv file where each line contains the following information:
<cluster_number>,<genome_name>,<left_flanking_region>,<right_flanking_region>

For example:
0, chlK2, left_100_200, right_300_400
Meaning:
Corresponding to cluster 0, the extragenic space in chlK2 is between 200bp and 300bp where the left flanking region is found between 100bp and 200bp and the right flanking region between 300bp and 400bp. The bp measurement is based on the original genome file provided.

Sometimes multiple flanking regions can exist, for example:
0, chlK2, left_100_200, right_300_400, right_800_900
Meaning:
This could arise in a situation where one of the REPINs in a cluster was matched only due to the left flanking sequence. So the right flanking sequence of this REPIN is different from the rest. Consider this situation:
Cluster 1:
(genome 1) A-REPIN1-B
(genome 2) A-REPIN2-B
(genome 3) A-REPIN3-B
(genome 4) A-REPIN4-B
(genome 5) A-REPIN5-X
(genome 6) no repin
For such a cluster, if flanking region 'X' is found in those genomes not containing a REPIN (genome 6), it will be shown along with flanking region 'B'. 'X' will not be shown as a flanking region for genomes 1-4 where a REPIN is present.
"""


def read_input():
    global clusters, pathmk, lhs_hits, rhs_hits, store_nearby_reps, flank_pairwise_dists, clusters, clus_cols, alltypes, genomenames, repbytype
    lhs_hits = pickle.load(open(f"{FILEPATH}/lhs_hits.p", "rb"))
    rhs_hits = pickle.load(open(f"{FILEPATH}/rhs_hits.p", "rb"))
    store_nearby_reps = pickle.load(
        open(f"{FILEPATH}/store_nearby_reps.p", "rb"))
    flank_pairwise_dists = pickle.load(
        open(f"{FILEPATH}/flank_pairwise_dists.p", "rb"))

    with open(CLUSTERFILENAME) as f:
        f_read = [i.split(" ") for i in f.read().split("\n")]
        for x in range(len(f_read)):
            for y in range(len(f_read[x])):
                try:
                    f_read[x][y] = int(f_read[x][y])
                except Exception:
                    pass
        clusters = {}
        clus_cols = {}
        for f in f_read:
            if len(f) <= 1:
                continue
            if f[0] not in clusters:
                clusters[f[0]] = []
            if f[1] in IGNOREDGENOMES:
                continue
            clusters[f[0]].append(f"{f[1]}_{f[2]}_{f[3]}")
            clus_cols[f"{f[1]}_{f[2]}_{f[3]}"] = f[4]

    pathmk = []
    with open(PATHFILENAME) as f:
        fread = f.read().split("\n")
        for i in range(len(fread)):
            if ">" in fread[i]:
                bob = fread[i][1:].split(" ")
                pathmk.append({'clus': None, 'path': None, 0: []})
                pathmk[-1]['clus'] = [int(x) for x in bob[0].split("_")]
                pathmk[-1]['path'] = bob[-1]
            else:
                pathmk[-1][0].append(fread[i].replace(" ", "_"))

    alltypes = list(set(clus_cols.values()))
    genomenames = list(set([v.split("_")[0] for v in clus_cols]))
    repbytype = {t: [] for t in alltypes}
    for k, v in clus_cols.items():
        repbytype[v].append(k)


def map_egs():
    # os.system("clear")
    # print("............Running extragenic_map.py")
    egs_map = {}
    for key, val in clusters.items():
        # print(f"Cluster {key}")
        egs_map[key] = {}
        # Genomes containing a REPIN and thus the EGS
        gen_with_repins = [v.split("_")[0] for v in val]
        # Genomes not containing a REPIN, unsure about EGS
        gen_wo_repins = list(set(genomenames) - set(gen_with_repins))
        gen_wo_repins = {v:{'left':{}, 'right':{}} for v in gen_wo_repins}
        for rep in val:
            # Add an entry of the repin name against the respective genome
            gen_name = rep.split("_")[0]
            rep_pos = rep.split("_")[1:]
            rep_pos  = [int(rep_pos[0]), int(rep_pos[1])]
            left_flank_rep = "left_{}_{}".format(rep_pos[0] - 1250, rep_pos[0] - 250)
            right_flank_rep = "right_{}_{}".format(rep_pos[1] + 250, rep_pos[1] + 1250)
            egs_map[key][gen_name] = [f"{left_flank_rep},{right_flank_rep}"]
            # Does the LHS of this REPIN exist in the genomes without a REPIN
            lhits = [x for x in lhs_hits[rep] if x[1] in gen_wo_repins]
            for item in lhits:
                gen_wo_repins[item[1]]['left'][f"{item[4]}_{item[5]}"] = gen_name
            # Does the RHS of this REPIN exist in the genomes without a REPIN
            rhits = [x for x in rhs_hits[rep] if x[1] in gen_wo_repins]
            for item in rhits:
                gen_wo_repins[item[1]]['right'][f"{item[4]}_{item[5]}"] = gen_name
        # If REPINs have left egs A and B. The gen_wo_repins[genome1]['left'] will have entries of multiple genomes each for A and B. How to pick one or chose one?
        for k2, v2 in gen_wo_repins.items():
            left = list(v2['left'].keys())
            right = list(v2['right'].keys())
            notuniq_left = []
            notuniq_right = []
            for i in range(len(left)):
                if i in notuniq_left:
                    continue
                for j in range(i + 1, len(left)):
                    if j in notuniq_left:
                        continue
                    a1,a2 = left[i].split("_")
                    b1,b2 = left[j].split("_")
                    a1 = int(a1)
                    a2 = int(a2)
                    b1 = int(b1)
                    b2 = int(b2)
                    t1 = a1 <= max(b1,b2) and a1 >= min(b1,b2)
                    t2 = a2 <= max(b1,b2) and a2 >= min(b1,b2)
                    if t1 or t2:
                        # i and j are the same
                        notuniq_left.append(j)
            left = [left[x] for x in range(len(left)) if x not in notuniq_left]
            
            for i in range(len(right)):
                if i in notuniq_right:
                    continue
                for j in range(i + 1, len(right)):
                    if j in notuniq_right:
                        continue
                    a1,a2 = right[i].split("_")
                    b1,b2 = right[j].split("_")
                    a1 = int(a1)
                    a2 = int(a2)
                    b1 = int(b1)
                    b2 = int(b2)
                    t1 = a1 <= max(b1,b2) and a1 >= min(b1,b2)
                    t2 = a2 <= max(b1,b2) and a2 >= min(b1,b2)
                    if t1 or t2:
                        # i and j are the same
                        notuniq_right.append(j)
            right = [right[x] for x in range(len(right)) if x not in notuniq_right]
            
            left = [f"left_{x}" for x in left]
            right = [f"right_{x}" for x in right]
            egs_map[key][k2] = left + right
        
    
    printstring = []
    for clusnum, val in egs_map.items():
        for gen, flank in val.items():
            printflank = ",".join(flank)
            printstring.append(f"{clusnum},{gen},{printflank}")
    with open(OUTPUTFILENAME, "w+") as f:
        f.write("\n".join(printstring))


def main():
    read_input()
    map_egs()
    print(READMESTR)

if __name__ == "__main__":
    main()
