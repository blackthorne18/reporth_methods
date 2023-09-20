import pickle
import numpy as np
import matplotlib.pyplot as plt
from prajwaltools import timethis


@timethis
def read_input():
    global lhs_hits, rhs_hits, store_nearby_reps, flank_pairwise_dists, clusters
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


@timethis
def graph1():
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
def graph2():
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
    # clus_graph['L'] = {key: val for key,
    #                    val in clus_graph['L'].items() if val > 0}
    # clus_graph['R'] = {key: val for key,
    #                    val in clus_graph['R'].items() if val > 0}

    # cgticks_L = [str(x) for x in list(clus_graph['L'].keys())]
    # cgticks_R = [str(x) for x in list(clus_graph['R'].keys())]

    # clussize_L = [len(clusters[key]) for key in clus_graph['L']]
    # clussize_R = [len(clusters[key]) for key in clus_graph['R']]

    # plt.subplot(121)
    # plt.bar(range(len(cgticks_L)), clus_graph['L'].values())
    # plt.scatter(range(len(cgticks_L)), clussize_L, color='red')
    # plt.xticks(range(len(cgticks_L)), cgticks_L, rotation=60)
    # plt.yticks(range(30))
    # plt.ylabel("Number of genomes")
    # plt.xlabel("Cluster Number")
    # plt.title("1.A. Left Flanking Sequnce")
    # plt.subplot(122)
    # plt.bar(range(len(cgticks_R)), clus_graph['R'].values())
    # plt.scatter(range(len(cgticks_R)), clussize_R, color='red')
    # plt.xticks(range(len(cgticks_R)), cgticks_R, rotation=60)
    # plt.yticks(range(30))
    # plt.ylabel("Number of genomes")
    # plt.xlabel("Cluster Number")
    # plt.title("1.B. Right Flanking Sequnce")
    # plt.suptitle(
    #     "Number of genomes in a cluster that have a potential paralog of a flanking sequence")
    # plt.show()
    # ------------------------------------------------------------

    # ------------------------------------------------------------
    # Graph by Genome

    genL = {k: v for k, v in gen_graph['L'].items() if v > 0}
    genR = {k: v for k, v in gen_graph['R'].items() if v > 0}

    plt.subplot(121)
    plt.bar(range(len(genL)), list(genL.values()), align='center')
    plt.xticks(range(len(genL)), list(genL.keys()), rotation=90)
    plt.yticks(range(5))
    plt.ylabel("Number of potential paralogs")
    plt.xlabel("Genome")
    plt.title("1.A. Left Flanking Sequnce")
    plt.subplot(122)
    plt.bar(range(len(genR)), list(genR.values()), align='center')
    plt.xticks(range(len(genR)), list(genR.keys()), rotation=90)
    plt.yticks(range(5))
    plt.ylabel("Number of potential paralogs")
    plt.xlabel("Genome")
    plt.title("1.B. Right Flanking Sequnce")
    plt.suptitle("Number of potential paralogs coming from different genomes")
    plt.tight_layout()
    plt.show()
    # ------------------------------------------------------------


def main():
    read_input()
    # graph1()
    graph2()


if __name__ == "__main__":
    main()
