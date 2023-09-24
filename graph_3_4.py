import pickle
import numpy as np
import matplotlib.pyplot as plt
from prajwaltools import timethis


@timethis
def read_input():
    global clusters, flank_pairwise_dists, nearby_reps, lhs_hits, rhs_hits
    lhs_hits = pickle.load(open("./temp/lhs_hits.p", "rb"))
    rhs_hits = pickle.load(open("./temp/rhs_hits.p", "rb"))
    nearby_reps = pickle.load(open("./temp/store_nearby_reps.p", "rb"))
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
def graph3():
    lefty, righty = {}, {}
    for key, val in nearby_reps['L'].items():
        ingens = [k for k, v in val.items() if len(v) == 0]
        if key not in lefty:
            lefty[key] = []
        # If a genome is added to lefty['repin'] it means
        # that genome DOES NOT HAVE a corresponding REPIN
        lefty[key].extend(ingens)
    for key, val in nearby_reps['R'].items():
        ingens = [k for k, v in val.items() if len(v) == 0]
        if key not in righty:
            righty[key] = []
        # If a genome is added to righty['repin'] it means
        # that genome DOES NOT HAVE a corresponding REPIN
        righty[key].extend(ingens)

    allgens = list(set([rep.split("_")[0]
                        for val in clusters.values() for rep in val]))
    no_repins = {}
    emptynest = {}
    for key, val in clusters.items():
        gens = [rep.split("_")[0] for rep in val]
        no_repins[key] = list(set(allgens) - set(gens))
        egs_exists = []
        for rep in val:
            try:
                egs_exists.extend(lefty[rep])
            except Exception:
                pass
            try:
                egs_exists.extend(righty[rep])
            except Exception:
                pass
        egs_exists = list(set(egs_exists))
        egs_norepin = list(set(egs_exists) - set(gens))

        emptynest[key] = len([x for x in egs_norepin if x != "chlPCL1606"])

    emptynest = {k: v for k, v in emptynest.items() if v > 0}

    plt.hist(emptynest.values(), bins=42)
    plt.xticks(range(1, 43))
    plt.xlabel("Genomes having EGS but not REPIN")
    plt.ylabel("Number of Clusters")
    plt.title("Distribution of extragenic spaces not having a REPIN")
    plt.show()


def graph4():
    # swap = {"right": "left", "left": "right"}
    lone, rone = {}, {}
    for key, val in nearby_reps['L'].items():
        items = [v for v in val.values() if len(v) > 0]
        lone[key] = {}
        for its in items:
            item = its[0]
            if key in item:
                continue
            gen = item[0].split("_")[0]
            nums = item[0].split("_")[1:]
            lone[key][gen] = [int(x) for x in nums]
    for key, val in nearby_reps['R'].items():
        items = [v for v in val.values() if len(v) > 0]
        rone[key] = {}
        for its in items:
            item = its[0]
            if key in item:
                continue
            gen = item[0].split("_")[0]
            nums = item[0].split("_")[1:]
            rone[key][gen] = [int(x) for x in nums]

    unreqL = {}
    unreqR = {}
    trashcan = []
    EGSWIDTH = 600
    for key, val in clusters.items():
        key = 95
        val = clusters[key]
        print(key, len(val))
        print(val)
        unreqL[key] = []
        unreqR[key] = []
        for x in range(len(val)):
            if val[x] not in rone or val[x] not in lone:
                trashcan.append(val[x])
                continue
            for y in range(x + 1, len(val)):
                if val[y] not in rone or val[y] not in lone:
                    trashcan.append(val[y])
                    continue
                genx = val[x].split("_")[0]
                geny = val[y].split("_")[0]
                xnums = [int(p) for p in val[x].split("_")[1:]]
                ynums = [int(p) for p in val[y].split("_")[1:]]
                if genx == geny:
                    diff3 = min(abs(xnums[0] - ynums[0]), abs(xnums[0] - ynums[1]),
                                abs(xnums[1] - ynums[0]), abs(xnums[1] - ynums[1]))
                    if diff3 > EGSWIDTH:
                        print(
                            f"Error in {key} with {genx}:{diff3}")

                checker, diff1, diff2, diff3 = "", 0, 0, 0
                if genx in lone[val[y]]:
                    dnum1 = lone[val[y]][genx]
                    diff1 = min(abs(xnums[0] - dnum1[0]), abs(xnums[0] - dnum1[1]),
                                abs(xnums[1] - dnum1[0]), abs(xnums[1] - dnum1[1]))
                    if diff1 <= EGSWIDTH:
                        checker += 'y'
                if geny in lone[val[x]]:
                    dnum2 = lone[val[x]][geny]
                    diff2 = min(abs(ynums[0] - dnum2[0]), abs(ynums[0] - dnum2[1]),
                                abs(ynums[1] - dnum2[0]), abs(ynums[1] - dnum2[1]))
                    if diff2 <= EGSWIDTH:
                        checker += 'x'
                if checker == 'x':
                    unreqL[key].append(genx)
                if checker == 'y':
                    unreqL[key].append(geny)

                checker, diff1, diff2, diff3 = "", 0, 0, 0
                if genx in rone[val[y]]:
                    dnum1 = rone[val[y]][genx]
                    diff1 = min(abs(xnums[0] - dnum1[0]), abs(xnums[0] - dnum1[1]),
                                abs(xnums[1] - dnum1[0]), abs(xnums[1] - dnum1[1]))
                    if diff1 <= EGSWIDTH:
                        checker += 'y'
                if geny in rone[val[x]]:
                    dnum2 = rone[val[x]][geny]
                    diff2 = min(abs(ynums[0] - dnum2[0]), abs(ynums[0] - dnum2[1]),
                                abs(ynums[1] - dnum2[0]), abs(ynums[1] - dnum2[1]))
                    if diff2 <= EGSWIDTH:
                        checker += 'x'
                if checker == 'x':
                    unreqR[key].append(genx)
                if checker == 'y':
                    unreqR[key].append(geny)
        exit()

    bygen = {}
    for key, val in clusters.items():
        unreq_both = set(unreqL[key]).union(set(unreqR[key]))
        unreq_both -= {"chlPCL1606"}
        if key == 95:
            print(unreq_both)
            print(unreqL[95], unreqR[95])
            exit()
        if len(unreq_both) > 0:
            print(key, len(val), len(unreq_both))
        for gen in unreq_both:
            if gen not in bygen:
                bygen[gen] = 0
            bygen[gen] += 1

    exit()
    bygen = {k[3:]: v for k, v in bygen.items()}
    del bygen["PCL1606"]

    plt.bar(range(len(bygen)), bygen.values(), tick_label=list(bygen.keys()))
    plt.xticks(rotation=90)
    plt.xlabel("Genomes")
    plt.ylabel("Number of unreciprocated connection")
    plt.title("Number of best hits of a flanking sequence that are not reciprocated")
    plt.tight_layout()
    plt.show()


def main():
    read_input()
    # graph3()
    graph4()


if __name__ == "__main__":
    main()
