import pickle
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import collections
import warnings
from copy import deepcopy

FILEPATH = "./temp/cluster_output_Sep20_875/"
CLUSTERFILENAME = FILEPATH + "/clusters_Sep20.txt"
PATHFILENAME = FILEPATH + "/path_making_Sep20.txt"
IGNOREDGENOMES = ['chlPCL1606']
FIGLOCATION = './temp/images/'
plt.rcParams.update({'font.size': 20})
LABELFONTSIZE = 24
FIGUREDIMENSION = (12, 8)
# Only for figure 7 use (14, 8)



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


def graph2():
    for keep_type in ['all', 'type0', 'type1', 'type2']:
        copy_clusters = deepcopy(clusters)
        if keep_type != 'all':
            for key in copy_clusters:
                copy_clusters[key] = [v for v in copy_clusters[key] if v in repbytype[keep_type]]
            _types = []
            for key, val in copy_clusters.items():
                for v in val:
                    _types.append(clus_cols[v])

        cluslen = {}
        for key, val in copy_clusters.items():
            gens = [v.split("_")[0] for v in val]
            cluslen[key] = len(list(set(gens)))

        colorguide = {
            'type0': '#4682b4', #Blue
            'type1': '#8b0000', #Red
            'type2': '#008b00', #Green
            'all': 'black'
        }[keep_type]
        yax = cluslen.values()
        print(keep_type, sum(yax))
        ybins = range(1, max(yax) + 1)
        ybins = range(1, 42)
        ybin_alternate = [str(x) if x % 2 == 0 or x == 1 else "" for x in ybins]
        fig, ax = plt.subplots(figsize=FIGUREDIMENSION)
        plt.hist(yax, bins=ybins, color=colorguide)
        plt.xticks(ybins, ybin_alternate)
        plt.ylim(0, 300)
        plt.xlabel("Number of genomes present in a cluster", fontsize=LABELFONTSIZE)
        plt.ylabel("Number of clusters", fontsize=LABELFONTSIZE)
        plt.title(keep_type)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        # plt.show()
        plt.savefig(FIGLOCATION + f'graph2_{keep_type}.pdf', dpi=500, format='pdf')
        plt.figure().clear()
        plt.close()
        plt.cla()
        plt.clf()


def graph2b():
    # The number of times each extragenic space occurs in P. chlororaphis for REPINs that occur in only one genome.
    oneset = {}
    reversekey = {}
    EGS_THRESHOLD = 700

    for key, val in clusters.items():
        gen = [v.split("_")[0] for v in val]
        if len(set(gen)) > 1:
            continue
        # Only looking at clusters of size 1
        if len(val) == 0:
            continue
        oneset[key] = val[0]
        reversekey[val[0]] = key

    egs_exists_check = {key: {gen: {'left': None, 'right': None}
                              for gen in genomenames} for key in oneset}

    for key, val in lhs_hits.items():
        if key not in oneset.values():
            continue
        clusnum = reversekey[key]
        for item in val:
            if item[1] in IGNOREDGENOMES:
                continue
            egs_exists_check[clusnum][item[1]]['left'] = [item[4], item[5]]
    for key, val in rhs_hits.items():
        if key not in oneset.values():
            continue
        clusnum = reversekey[key]
        for item in val:
            if item[1] in IGNOREDGENOMES:
                continue
            egs_exists_check[clusnum][item[1]]['right'] = [item[4], item[5]]

    egs_exists = {key: [0] for key in oneset}
    for key, val in egs_exists_check.items():
        confirm = []
        # print(key)
        for gen, item in val.items():
            # print(gen, item)
            if item['left'] is None or item['right'] is None:
                continue
            d1 = abs(item['left'][0] - item['right'][0])
            d2 = abs(item['left'][1] - item['right'][0])
            d3 = abs(item['left'][0] - item['right'][1])
            d4 = abs(item['left'][1] - item['right'][1])
            diff = min(d1, d2, d3, d4)
            if diff <= EGS_THRESHOLD:
                confirm.append(gen)
        egs_exists[key] = list(set(confirm))
        # exit()

    egs_exists = {k: len(v) for k, v in egs_exists.items()}
    yax = list(egs_exists.values())
    ybins = range(1, max(yax) + 1)
    ybin_alternate = [str(x) if x % 2 == 0 or x == 1 else "" for x in ybins]

    fig, ax = plt.subplots(figsize=FIGUREDIMENSION)
    plt.hist(yax, bins=ybins, color='black')
    plt.xticks(ybins, ybin_alternate)
    plt.yticks(range(0, 100, 10))
    plt.xlabel(
        "Number of genomes that contain the given extragenic space", fontsize=LABELFONTSIZE)
    plt.ylabel("Number of clusters", fontsize=LABELFONTSIZE)
    # plt.title("The number of times each extragenic space occurs in P. chlororaphis for REPINs that occur in only one genome")
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    # plt.show()
    plt.tight_layout()
    plt.savefig(FIGLOCATION + 'graph2_b.pdf', dpi=500, format='pdf')
    plt.figure().clear()
    plt.close()
    plt.cla()
    plt.clf()


def graph4():
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
    # I expect to see RuntimeWarnings in this block
    # This is because for clusters of size 1, we will be taking mean
    # of an empty slice
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        for key in lg:
            lg[key] = np.mean(lg[key])
            rg[key] = np.mean(rg[key])

    bins = range(90, 101)
    either = list(lg.values()) + list(rg.values())
    either = [x for x in either if str(x) != 'nan']
    # How many of these are above 97% and how many between 90-91%
    n96 = round(len([x for x in either if x >= 96]) / len(either), 3)
    n90 = round(len([x for x in either if x < 92]) / len(either), 3)
    print(f'Graph4\nProportion of hits >=96%: {n96}\nProportion of hits <92%: {n90}')
    print(f'AUC', len(either))

    fig, ax = plt.subplots(figsize=FIGUREDIMENSION)
    plt.hist(either, bins=bins, density=True, color='black')
    # plt.hist(lg.values(), bins=bins, alpha=0.5, label='L')
    # plt.hist(rg.values(), bins=bins, alpha=0.5, label='R')
    plt.xticks(range(90, 102, 2))
    # plt.hist([list(lg.values()), list(rg.values())], bins=range(
    #     90, 101), label=['L', 'R'])
    # plt.xticks(range(90, 101))
    # plt.legend(loc='upper right')
    plt.ylabel("Proportion of clusters", fontsize=LABELFONTSIZE)
    plt.xlabel("Average similarity", fontsize=LABELFONTSIZE)
    # plt.title("Average Similarity of Flanking Sequences Within A Cluster")
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    # plt.show()
    fig.savefig(FIGLOCATION + 'graph4.pdf', dpi=500, format='pdf')
    plt.figure().clear()
    plt.close()
    plt.cla()
    plt.clf()


def graph5():
    phylo_gens = [
            'ChPhzTR18', 'ChPhzS24', 'ChPhzTR38', 'ChPhzTR39', 'PA23', 'ChPhzS23', '66', 'O6', 'Lzh-T5', 'ChPhzTR36',
            '6698', 'C50', 'P2', '19603', '2210', 'Q16', 'StFRB508', 'JD37', 'CW2', '449', '464', 'M12', 'K27', 'M71',
            'Pb-St2', '189', '17415', '50083', 'TAMOak81', 'B25', 'ChPhzS140', 'PCL1391', 'PCL1607', 'ToZa7', 'ChPhzS135',
            'DTR133', 'SLPH10', 'ChPhzTR44', 'ZJU60', '17809', '17411', '21509'
    ]

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
    gbe = gbl + gbr
    lbb = range(len(bygen['both']))
    gennames = [key for key in bygen['both']]
    gsum['either'] = gsum['left'] + gsum['right']
    del gsum['left']
    del gsum['right']

    plt.figure(figsize=(20, 12))
    plt.subplot(1, 3, 1)
    # This subplot in numbers
    nboth = gsum['both'] / (gsum['both'] + gsum['either'])
    neither = 1 - nboth
    print(f'Graph5\nProportion of both:{nboth}\nProportion of either:{neither}')
    plt.bar(range(len(gsum)), gsum.values(), color=['black', 'gray'])
    plt.xticks(range(len(gsum)), ["Both", "Either"])
    plt.ylabel("Number of REPINs", fontsize = LABELFONTSIZE)
    # plt.title("1.A. Left Flanking Sequnce")
    plt.subplot(1, 3, (2, 3))
    plt.bar(lbb, gbb, label='Both', color='black')
    plt.bar(lbb, gbe, bottom=gbb, label='Either', color='gray')
    plt.xticks(lbb, gennames, rotation=90)
    plt.ylabel("Number of REPINs", fontsize = LABELFONTSIZE)
    plt.legend()
    # plt.title("1.B. Right Flanking Sequnce")
    # plt.suptitle(
    #     "REPINs clustered based on both, or one of the flanking sequences")
    plt.tight_layout()
    # plt.show()
    plt.savefig(FIGLOCATION + 'graph5.pdf', dpi=500, format='pdf')
    plt.figure().clear()
    plt.close()
    plt.cla()
    plt.clf()


def graph6():
    reversekey = {}
    for item in pathmk:
        for rep in item[0]:
            reversekey[rep] = item['path']

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

    clus_graph = {}
    for key, val in clusters.items():
        if len(val) < 2:
            continue
        clus_graph[key] = 0
        for rep in val:
            _match = ''
            if rep in lhs:
                _match = 'left'
            if rep in rhs:
                if _match == 'left':
                    _match = 'both'
                else:
                    _match = 'right'

            # If this REPIN was merged because of the '_match' flanking sequence i.e. left or right
            # and that same flanking sequence is likely present twice, then this is false positive
            # Or if it was matched because of both flanking sequences and one of them is a paralog
            # then this is false positive
            if _match != '':
                if reversekey[rep] == _match or reversekey[rep] == 'both':
                    clus_graph[key] += 1
                    print(f'Cluster:{key}\nMerged based on: {reversekey[rep]}\tPotential Paralog: {_match}')
    # ------------------------------------------------------------
    # Graph by cluster
    clus_graph = {key: val for key, val in clus_graph.items() if val > 0}
    clus_graph = dict(collections.OrderedDict(sorted(clus_graph.items())))

    cgticks = [str(x) for x in list(clus_graph.keys())]
    clussize = {key: len(clusters[key]) for key in clus_graph}
    fig, ax = plt.subplots(figsize=FIGUREDIMENSION)

    # Sorting dataset by cluster size
    clussize = dict(sorted(clussize.items(), key=lambda item: item[1]))
    cgboth = sorted(clus_graph.items(), key=lambda pair: list(clussize.keys()).index(pair[0]))
    cgboth = {x[0]:x[1] for x in cgboth}

    plt.bar(range(len(cgticks)), clussize.values(), color='gray', label='Cluster size')
    plt.bar(range(len(cgticks)), cgboth.values(), color='black', label='Potential paralogs')
    plt.xticks(range(len(cgticks)), cgboth.keys(), rotation=60)
    yticks = [str(x) if x%2==0 else '' for x in range(30)]
    plt.yticks(range(len(yticks)), yticks)
    plt.ylabel("Number of genomes", fontsize=LABELFONTSIZE)
    plt.xlabel("Cluster number", fontsize=LABELFONTSIZE)
    plt.legend()
    plt.tight_layout()
    # plt.title("Number of genomes in a cluster that have a potential paralog of a flanking sequence")
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    # plt.show()
    plt.savefig(FIGLOCATION + 'graph6.pdf', dpi=500, format='pdf')
    plt.figure().clear()
    plt.close()
    plt.cla()
    plt.clf()


def graph7():
    # Histogram of distances between REPINs when there are two or more from the same genome
    # Capture all instances, highlight the ones that are above 1000bp
    HIGH_THRESH=1000
    wgroups = {}
    for key, val in clusters.items():
        samegens = {}
        wgroups[key] = []
        for rep in val:
            gen = rep.split('_')[0]
            if gen not in samegens:
                samegens[gen] = []
            samegens[gen].append(rep)
        samegens = [v for k,v in samegens.items() if len(v) > 1]
        for item in samegens:
            for x, r1 in enumerate(item):
                for y, r2 in enumerate(item):
                    if x == y:
                        continue
                    da = r1.split('_')
                    da = [int(da[1]), int(da[2])]
                    db = r2.split('_')
                    db = [int(db[1]), int(db[2])]
                    diffab = min(abs(da[1] - db[0]), abs(db[1] - da[0]))
                    wgroups[key].append([diffab, r1, r2])
    # Capturing the problematic ones first and all the ones to print
    probs = {}
    wgvals = []
    for key, itemsx in wgroups.items():
        if len(itemsx) == 0:
            continue
        for item in itemsx:
            wgvals.append(item[0])
            if item[0] > HIGH_THRESH:
                if key not in probs:
                    probs[key] = []
                probs[key].append(item)
    # print(probs)
    # Plotting the graph
    fig, ax = plt.subplots(figsize=FIGUREDIMENSION)
    plt.hist(wgvals, color='black')
    plt.ylabel('Number of occurrences',fontsize=LABELFONTSIZE)
    plt.xlabel('Distance between REPINs from the same genome within the same cluster',fontsize=LABELFONTSIZE)
    # plt.legend()
    plt.tight_layout()
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    # plt.show()
    plt.savefig(FIGLOCATION + 'graph7.pdf', dpi=500, format='pdf')
    plt.figure().clear()
    plt.close()
    plt.cla()
    plt.clf()


def main():
    # Standard Reading of Inputs
    read_input()

    # Primary Functions
    # graph2()
    # graph2b()
    # graph4()
    # graph5()
    # graph6()
    graph7()

if __name__ == "__main__":
    main()
