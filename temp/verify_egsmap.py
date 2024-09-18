#!/usr/local/bin/python3.10

MAPFILE = "clusters_extragenic_map.csv"

with open(MAPFILE, "r") as f:
    mpgen = {}
    mpall = {}
    for y in f.read().split("\n"):
        if len(y) < 1:
            continue
        line = y.split(",")
        num = int(line[0])
        gen = line[1]
        if num not in mpgen:
            mpgen[num] = []
            mpall[num] = []
        mpgen[num].append(gen)
        mpall[num].append(line)

# How many entries per cluster?
_clussize = {k: len(v) for k,v in mpgen.items()}
clussize = set(_clussize.values())
if len(clussize) == 1:
    print(f"---1- All clusters have the same number of entries {clussize}")
else:
    clussize_counter = {x:list(_clussize.values()).count(x) for x in clussize}
    _smaller = min(clussize_counter, key=clussize_counter.get)
    smaller = [k for k,v in _clussize.items() if v == _smaller]
    print(f"---1- The following clusters {smaller} have reduced entries ({_smaller})")
 

