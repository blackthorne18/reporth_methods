import pickle
import numpy as np
import matplotlib.pyplot as plt
from prajwaltools import timethis

"""
Graph 3:
    (x axis) Plot graph  : for each cluster
    (y axis) showing     : How many REPINs are included because both flanking sequences match, and how many for LHS, RHS
    Example              : In a cluster of six REPINs, say A,B,C are grouped because both flanks match, then D because LHS
                           matches and E and F because RHSs match. Then the stacked bar graph would have 3, 1, 2 stacked
    Note                 : Also show as percentage or one graph for each type of stacking
    Data Needed          : post BLAST clustering stage
"""


@timethis
def read_input():
    global clusters, pathmk, lhs_hits, rhs_hits, store_nearby_reps, flank_pairwise_dists, clusters
    lhs_hits = pickle.load(open("./temp/lhs_hits.p", "rb"))
    rhs_hits = pickle.load(open("./temp/rhs_hits.p", "rb"))
    store_nearby_reps = pickle.load(open("./temp/store_nearby_reps.p", "rb"))
    flank_pairwise_dists = pickle.load(
        open("./temp/flank_pairwise_dists.p", "rb"))

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
                pathmk[-1]['path'] = bob[-1]
            else:
                pathmk[-1][0].append(fread[i].replace(" ", "_"))


@timethis
def graph1b():
    """
    Note: lhs and rhs are structured as:
    lhs = {
        repin: [repinname, genome, pident, length, start, end]
    }
    """
    lhs, rhs = {}, {}
    lg, rg = {}, {}
    for key, val in lhs_hits.items():
        lhs[key] = {}
        for item in val:
            lhs[key][item[1]] = item[2:]
    for key, val in rhs_hits.items():
        rhs[key] = {}
        for item in val:
            rhs[key][item[1]] = item[2:]

    for key, clus in clusters.items():
        if len(clus) < 2:
            continue
        lg[key] = []
        rg[key] = []
        for x in range(len(clus)):
            for y in range(x + 1, len(clus)):
                gx = clus[x].split("_")[0]
                gy = clus[y].split("_")[0]
                try:
                    lg[key].append(lhs[clus[x]][gy][0])
                except Exception:
                    # lg[key].append(0)
                    pass
                try:
                    lg[key].append(lhs[clus[y]][gx][0])
                except Exception:
                    # lg[key].append(0)
                    pass
                try:
                    rg[key].append(rhs[clus[x]][gy][0])
                except Exception:
                    # rg[key].append(0)
                    pass
                try:
                    rg[key].append(rhs[clus[y]][gx][0])
                except Exception:
                    # rg[key].append(0)
                    pass

    for key in lg:
        lg[key] = np.mean(lg[key])
        rg[key] = np.mean(rg[key])

    # for key in lg:
    #     if lg[key] < 50:
    #         print(key, len(clusters[key]))
    # exit()
    plt.subplot(121)
    plt.hist(lg.values(), bins=10)
    plt.xticks(range(90, 102, 2))
    plt.ylabel("Number of Clusters")
    plt.xlabel("Average similarity")
    plt.title("1.A. Left Flanking Sequnce")
    plt.subplot(122)
    plt.hist(rg.values(), bins=10)
    plt.xticks(range(90, 102, 2))
    plt.ylabel("Number of Clusters")
    plt.xlabel("Average similarity")
    plt.title("1.B. Right Flanking Sequnce")
    plt.suptitle(
        "Average Similarity of Flanking Sequences Within A Cluster")
    plt.show()


@timethis
def graph4():
    lhs, rhs = {}, {}
    mpara_L = {}
    mpara_R = {}
    for key, val in lhs_hits.items():
        manygen = {}
        for item in val:
            if item[1] not in manygen:
                manygen[item[1]] = []
            manygen[item[1]].append(item[2:])
        mpara_L[key] = {}
        for k2, v2 in manygen.items():
            if len(v2) > 1:
                mpara_L[key][k2] = v2

    for key, val in rhs_hits.items():
        manygen = {}
        for item in val:
            if item[1] not in manygen:
                manygen[item[1]] = []
            manygen[item[1]].append(item[2:])
        mpara_R[key] = {}
        for k2, v2 in manygen.items():
            if len(v2) > 1:
                mpara_R[key][k2] = v2

    for key, val in mpara_L.items():
        if len(val) > 1:
            lhs[key] = list(val.keys())
    for key, val in mpara_R.items():
        if len(val) > 1:
            rhs[key] = list(val.keys())

    clus_graph = {'L': {}, 'R': {}}
    gen_graph = {'L': {}, 'R': {}}
    for key, val in clusters.items():
        if len(val) < 2:
            continue
        clus_graph['L'][key] = 0
        clus_graph['R'][key] = 0
        for rep in val:
            gen = rep.split("_")[0]
            if gen not in gen_graph['L']:
                gen_graph['L'][gen] = 0
            if gen not in gen_graph['R']:
                gen_graph['R'][gen] = 0
            if rep in lhs:
                clus_graph['L'][key] += 1
                gen_graph['L'][gen] += 1
            if rep in rhs:
                clus_graph['R'][key] += 1
                gen_graph['R'][gen] += 1

    # ------------------------------------------------------------
    # Graph by cluster
    clus_graph['L'] = {key: val for key,
                       val in clus_graph['L'].items() if val > 0}
    clus_graph['R'] = {key: val for key,
                       val in clus_graph['R'].items() if val > 0}

    cgticks_L = [str(x) for x in list(clus_graph['L'].keys())]
    cgticks_R = [str(x) for x in list(clus_graph['R'].keys())]

    clussize_L = [len(clusters[key]) for key in clus_graph['L']]
    clussize_R = [len(clusters[key]) for key in clus_graph['R']]

    plt.subplot(121)
    plt.bar(range(len(cgticks_L)), clus_graph['L'].values())
    plt.scatter(range(len(cgticks_L)), clussize_L, color='red')
    plt.xticks(range(len(cgticks_L)), cgticks_L, rotation=60)
    plt.yticks(range(30))
    plt.ylabel("Number of genomes")
    plt.xlabel("Cluster Number")
    plt.title("1.A. Left Flanking Sequnce")
    plt.subplot(122)
    plt.bar(range(len(cgticks_R)), clus_graph['R'].values())
    plt.scatter(range(len(cgticks_R)), clussize_R, color='red')
    plt.xticks(range(len(cgticks_R)), cgticks_R, rotation=60)
    plt.yticks(range(30))
    plt.ylabel("Number of genomes")
    plt.xlabel("Cluster Number")
    plt.title("1.B. Right Flanking Sequnce")
    plt.suptitle(
        "Number of genomes in a cluster that have a potential paralog of a flanking sequence")
    plt.show()
    # ------------------------------------------------------------

    # ------------------------------------------------------------
    # Graph by Genome

    # genL = {k: v for k, v in gen_graph['L'].items() if v > 0}
    # genR = {k: v for k, v in gen_graph['R'].items() if v > 0}

    # plt.subplot(121)
    # plt.bar(range(len(genL)), list(genL.values()), align='center')
    # plt.xticks(range(len(genL)), list(genL.keys()), rotation=90)
    # plt.yticks(range(5))
    # plt.ylabel("Number of potential paralogs")
    # plt.xlabel("Genome")
    # plt.title("1.A. Left Flanking Sequnce")
    # plt.subplot(122)
    # plt.bar(range(len(genR)), list(genR.values()), align='center')
    # plt.xticks(range(len(genR)), list(genR.keys()), rotation=90)
    # plt.yticks(range(5))
    # plt.ylabel("Number of potential paralogs")
    # plt.xlabel("Genome")
    # plt.title("1.B. Right Flanking Sequnce")
    # plt.suptitle("Number of potential paralogs coming from different genomes")
    # plt.tight_layout()
    # plt.show()
    # ------------------------------------------------------------


@timethis
def graph3():
    phylo_gens = ['ChPhzTR38', 'ChPhzS24', 'ChPhzTR18', 'ChPhzTR39', 'PA23', 'ChPhzS23', '66', 'Lzh-T5', 'O6', 'ChPhzTR36', '6698', 'C50', 'P2', '2210', '19603', 'CW2', 'Q16', '464', '449', 'M12', 'StFRB508',
                  'JD37', 'M71', 'K27', 'Pb-St2', 'B25', 'SLPH10', 'ChPhzS135', 'DTR133', 'ToZa7', 'ChPhzTR44', 'PCL1607', 'PCL1391', 'ChPhzS140', '17809', '17411', 'ZJU60', '21509', '189', '17415', '50083', 'TAMOak81']

    reversekey = {}
    for item in pathmk:
        for rep in item[0]:
            reversekey[rep] = item['path']

    graphs = {'both': {}, 'left': {}, 'right': {}}
    bygen = {'both': {}, 'left': {}, 'right': {}}
    for key, val in clusters.items():
        if len(val) < 2:
            continue
        graphs['both'][key] = 0
        graphs['left'][key] = 0
        graphs['right'][key] = 0
        for rep in val:
            graphs[reversekey[rep]][key] += 1
            gen = rep.split("_")[0]
            if gen not in bygen[reversekey[rep]]:
                bygen[reversekey[rep]][gen] = 0
            bygen[reversekey[rep]][gen] += 1

    bygen['both'] = {newkey: bygen['both']["chl" + newkey]
                     for newkey in phylo_gens}
    bygen['left'] = {newkey: bygen['left']["chl" + newkey]
                     for newkey in phylo_gens}
    bygen['right'] = {newkey: bygen['right']["chl" + newkey]
                      for newkey in phylo_gens}

    gsum = {k: sum(v.values()) for k, v in graphs.items()}
    gbb = np.array(list(bygen['both'].values()))
    gbl = np.array(list(bygen['left'].values()))
    gbr = np.array(list(bygen['right'].values()))
    lbb = range(len(bygen['both']))
    gennames = [key for key in bygen['both']]

    plt.subplot(1, 3, 1)
    plt.bar(range(len(gsum)), gsum.values())
    plt.xticks(range(len(gsum)), ["Both", "Left", "Right"])
    plt.ylabel("Number of REPINs")
    # plt.title("1.A. Left Flanking Sequnce")
    plt.subplot(1, 3, (2, 3))
    plt.bar(lbb, gbb, label='Both')
    plt.bar(lbb, gbl, bottom=gbb, label='Left')
    plt.bar(lbb, gbr, bottom=gbb + gbl, label='Right')
    plt.xticks(lbb, gennames, rotation=90)
    plt.ylabel("Number of REPINs")
    plt.legend()
    # plt.title("1.B. Right Flanking Sequnce")
    plt.suptitle(
        "REPINs clustered based on both, or one of the flanking sequences")
    plt.tight_layout()
    plt.show()


def main():
    read_input()
    # graph1b()
    # graph3()
    graph4()


if __name__ == "__main__":
    main()
