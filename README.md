## Visualising Clusters on the Phylogenetic Tree
Run `visualise_clusters.r` to visualise REPINs found in a particular cluster on the phylogenetic tree of the 42 genomes. The variable `whichclus` is responsible for identifying the cluster number that you want to visualise.

## What is runtest?
Runs the clustering algorithm on the input file. Note this is a local instance of the clustering algorithm
> From the rarefan output folder in this directory.
`python3.10 rc_entry.py --repin rarefan_output/ --genomes ./input/genomes --reptypes 0,1,2`

> Once the above code is run, a file called 'sortedrepins.txt' is created
> in ./rarefan_output/ 
> This sortedrepins.txt has been moved to the input folder so that the following
> code can be run. This file contains all the REPINs found by RAREFAN from all genomes
> and is converting the rarefan output into a more readable form for this program
`python3 rc_entry.py --repin ./input/og_replist.txt --genomes ./input/genomes --reptypes 0,1,2`

> IMPORTANT
> The output from RAREFAN here contains the genome 'chlPCL1606' this has not been included in the analysis. From the
> 'sortedrepins.txt' file as mentioned above, all REPINs from 'chlPCL1606' has been removed. I have not edited the 
> RAREFAN output file to show the original dataset unmodified

A local version of REPORTH is used here for ease of testing, REPORTH installed from PyPI can also directly be used to generate clusters.

## How to use graph_main.py
Directly run the script as `python3 graph_main.py`. Each function generates a graph corresponding to the figure number in the REPORTH paper. The generated plots are stored in `temp/images/`.
All global variables such as those defining paths to files are all listed at the top. of the script.

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
