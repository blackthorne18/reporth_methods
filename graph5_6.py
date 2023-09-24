import pickle
import numpy as np
import matplotlib.pyplot as plt
from prajwaltools import timethis

"""
Graph 5:
    (x axis) Plot graph  : for each cluster
    (y axis) showing     : How many REPINs are included because both flanking sequences match, and how many for LHS, RHS
    Example              : In a cluster of six REPINs, say A,B,C are grouped because both flanks match, then D because LHS
                           matches and E and F because RHSs match. Then the stacked bar graph would have 3, 1, 2 stacked
    Note                 : Also show as percentage or one graph for each type of stacking
    Data Needed          : post BLAST clustering stage
"""


@timethis
def read_input():
    global clusters, pathmk
    with open("./output/cluster_output_Sep20_875/clusters_Sep20.txt") as f:
        f_read = [i.split(" ") for i in f.read().split("\n")]
        for x in range(len(f_read)):
            for y in range(len(f_read[x])):
                try:
                    f_read[x][y] = int(f_read[x][y])
                except Exception:
                    pass
        clusters = {}
        for f in f_read:
            if len(f) <= 1:
                continue
            if f[0] not in clusters:
                clusters[f[0]] = []
            clusters[f[0]].append(f"{f[1]}_{f[2]}_{f[3]}")

    pathmk = []
    with open("./output/cluster_output_Sep20_875/path_making_Sep20.txt") as f:
        fread = f.read().split("\n")
        for i in range(len(fread)):
            if ">" in fread[i]:
                bob = fread[i][1:].split(" ")
                pathmk.append({'clus': None, 'path': None, 0: []})
                pathmk[-1]['clus'] = [int(x) for x in bob[0].split("_")]
                pathmk[-1]['path'] = bob[1]
            else:
                pathmk[-1][0].append(fread[i].replace(" ", "_"))


@timethis
def graph5():
    for item in pathmk:
        pass


@timethis
def graph6():
    pass


def main():
    read_input()
    graph5()
    graph6()


if __name__ == "__main__":
    main()
