library(tidyverse)
library(ggtree)
library(glue)
library(treeio)
library(ggstance)
library(hash)
library(ggtreeExtra)

# Where are the input files stored?
TREEPATH <- "input/tree2.nwk"
CLUSPATH <- "temp/cluster_output_Sep20_875/clusters_Sep20.txt"
fortytwo <- c(12, 66, 74, 79, 96, 111, 141, 150, 159, 873, 880, 906, 912, 924, 929, 934, 936, 945, 949, 952)
clusters <- read.csv(CLUSPATH, sep=' ', header=FALSE)
colnames(clusters) <- c('cluster', 'genome', 'start', 'end', 'color', 'seq')
clusters <- clusters[c('cluster', 'genome', 'color')]
clusters <- clusters[clusters$cluster %in% fortytwo, ]
clusters$genome <- gsub('chl', '', clusters$genome)
# IMPORTANT NOTE
# THE CURRENT NEWICK FILE DOES NOT HAVE GENOMES STFRB508 AND 06. THIS WILL
# CAUSE ERRORS SINCE THE REPINS DATASET HAS IT BUT THE TREE DOES NOT
clusters <- clusters[clusters$genome != 'StFRB508', ]
clusters <- clusters[clusters$genome != 'O6', ]
clusters$cluster <- as.factor(clusters$cluster)
# Input data
tree <- read.tree(TREEPATH)
rayttree <- read.csv('./input/rayts.csv', sep=',')

# ---------------------------------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------------------------
xv <- 50
plot1 <- ggtree(tree, branch.length = 'none')  +
  geom_tiplab(align=TRUE, linesize=.5) + xlim(0, xv) + 
  geom_fruit(data=rayttree,geom=geom_tile,
             mapping = aes(y=genome, x=Color, fill=rayt),
             offset=0.2,pwidth=0.05) +
  geom_treescale()

for (cn in fortytwo){
  cname <- as.character(cn)
  dclus <- clusters[clusters$cluster == cname, ]
  plot1 <- plot1 + geom_fruit(data=dclus,geom=geom_tile,
               mapping = aes(y=genome, x=color, fill=color),
               offset=0.05,
               pwidth=0.02,
               axis.params=list(
                 axis="x",
                 title = cname,
                 title.size=6,
                 text.size = 0,
                 line.size = 0,
                 vjust = 0
               ))
}

plot1
# ---------------------------------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------------------------
