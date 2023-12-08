library(tidyverse)
library(ggtree)
library(glue)
library(treeio)
library(ggstance)
library(hash)
library(ggtreeExtra)

# Where are the input files stored?
TREEPATH <- "input/tree2.nwk"
CLUSPATH <- "input/clusters.txt"
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
clusters <- clusters[clusters$genome != 'PCL1606', ]
clusters$cluster <- as.factor(clusters$cluster)

# Testing cases where there are multiple REPINs of the same type in a cluster
# clusters[nrow(clusters) + 1,] = c('952', 'TAMOak81', 'type2')
# clusters[nrow(clusters) + 1,] = c('952', '189', 'type2')


# Input data
tree <- read.tree(TREEPATH)
rayttree <- read.csv('./input/rayts.csv', sep=',')

# ---------------------------------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------------------------
xv <- 55
plot1 <- ggtree(tree, branch.length = 'none')  +
  geom_tiplab(align=TRUE, linesize=.5) + xlim(0, xv) + 
  geom_fruit(data=rayttree,geom=geom_tile,
             mapping = aes(y=genome, x=Color, fill=rayt),
             pwidth=0.1, offset=0.25) +
  geom_treescale()

#fortytwo=c(924, 936, 952)
numClus=length(fortytwo)
for (cn in fortytwo){
  cname <- as.character(cn) #141 DOES NOT WORK
  dclus <- clusters[clusters$cluster == cname, ]
  
  genomelist <- unique(dclus$genome)
  for (gen in genomelist){
    dclus[dclus$genome == gen, ]$color <- ave(as.character(dclus[dclus$genome == gen, ]$color),
                                              dclus[dclus$genome == gen, ]$color, 
                                              FUN=function(x) if (length(x)>1) paste0(x[1], '_', seq_along(x)) else x[1])
  }
  
  maxR=max(max(table(dclus[,2])),length(unique(dclus[,3])))
  
  if(maxR==1){
    pw=0.4
  }else if(maxR==2){
    pw=0.05*maxR/3
  }else if(maxR==3){
    pw=0.1*2/3
  }else{
    pw=0.09
  }
  
  plot1 <- plot1 + geom_fruit(data=dclus,geom=geom_tile,
                              mapping = aes(y=genome, x=color, fill=color),pwidth=pw, offset=0.049,
                              axis.params=list(
                                axis="x",
                                title = cname,
                                title.size=6,
                                title.angle=-90,
                                title.height=0.05,
                                text.size = 0,
                                line.size = 0,
                              ))
}

color_values <- c("type0"="green4", "type0_1"="green4", "type0_2"="green",
                  "type1"="red4", "type1_1"="red4", "type1_2"="red",
                  "type2"="steelblue", "type2_1"="steelblue", "type2_2"="blue4")

plot1+scale_fill_manual(values=color_values)  + 
  vexpand(0.04)

# ---------------------------------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------------------------