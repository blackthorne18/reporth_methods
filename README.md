## Visualising Clusters on the Phylogenetic Tree
Run `visualise_clusters.r` to visualise REPINs found in a particular cluster on the phylogenetic tree of the 42 genomes. The variable `whichclus` is responsible for identifying the cluster number that you want to visualise.

## What is runtest?
Runs the clustering algorithm on the input file. Note this is a local instance of the clustering algorithm

## How to use graph_main.py
The figures are generated as per the function name.
`KEEPTYPE` refers to whether you want the analysis to be done on one type of REPIN or on all. Choose the appopriate repin type

The `temp/*.p` files are `pickle` files stored from the clustering program running. These are typically deleted on completion of the program but storing it makes it faster to work with the data.
If you want the clustering algorithm to not delete the temporary files, comment out line 730 of `rc_entry.p` and it will not delete these temporary files. They are needed for `graph_main.py` to run

## How to use Output_Graphs?

Run as 'python3 output_graphs.py --help' for details on the commands and options
Note: The <METAFILE_NAME> is the file in `output/cluster_output_Sep20_875/meta_cluster_Sep20.txt`

### 1. You want to look at the Network of connections between REPIN flanking sequences of a particular cluster?

`python3 output_graphs.py --file <METAFILE_NAME> --cluster <CLUSTER_NUMBER>`

### 2. If you want to look at histogram summary of the Quantification of Connectivity of Flanking Sequences of each and every cluster?
#### a. Based on Number of Cliques

`python3 output_graphs.py --file <METAFILE_NAME> --summary 0 --summary_types 1`

#### b. Based on Number of Connectivity

`python3 output_graphs.py --file <METAFILE_NAME> --summary 0 --summary_types 1`

#### Summary for each cluster, rather than a histogram distribution

`python3 output_graphs.py --file <METAFILE_NAME> --summary 1 --summary_types <0_or_1>`

### Caching
Once the program has been run, a directory is created where the temporary files are stored, making it faster to rerun the code for the same dataset. 

If you want it to rerun the code in its entirity without caching before plotting, use the '--recache 1' option at the end of the code.
