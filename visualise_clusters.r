library(tidyverse)
library(ggtree)
library(glue)
library(treeio)
library(ggstance)
library(hash)
library(ggtreeExtra)
library(showtext)
library(patchwork)
font_add(family = "dejavusans", regular = "input/dejavusans.ttf")
showtext_auto()

repgreen <- '#008b00'
repred   <- '#8b0000'
repblue  <- '#4682b4'

# Where are the input files stored?
TREEPATH <- "input/tree.nwk"
CLUSPATH <- "input/clusters.txt"
fortytwo <- c(13, 17, 79, 87, 92, 109, 125, 157, 159, 170, 179, 859, 878, 884, 896, 901, 904, 906, 914, 917)
clusters <- read.csv(CLUSPATH, sep=' ', header=FALSE)
colnames(clusters) <- c('cluster', 'genome', 'start', 'end', 'color', 'seq')
clusters <- clusters[c('cluster', 'genome', 'color')]
clusters <- clusters[clusters$cluster %in% fortytwo, ]
clusters$genome <- gsub('chl', '', clusters$genome)
# IMPORTANT NOTE
# THE CURRENT NEWICK FILE DOES NOT HAVE GENOMES STFRB508 AND 06. THIS WILL
# CAUSE ERRORS SINCE THE REPINS DATASET HAS IT BUT THE TREE DOES NOT
# clusters <- clusters[clusters$genome != 'StFRB508', ]
# clusters <- clusters[clusters$genome != 'O6', ]
clusters$cluster <- as.factor(clusters$cluster)

# Testing cases where there are multiple REPINs of the same type in a cluster
# clusters[nrow(clusters) + 1,] = c('952', 'TAMOak81', 'type2')
# clusters[nrow(clusters) + 1,] = c('952', '189', 'type2')

# Input data
tree <- read.tree(TREEPATH)
rayttree <- read.csv('./input/rayts.csv', sep=',')

# Adding 'empty spaces with PCL1606
# tree <- drop.tip(tree, 'PCL1606')
clusters <- clusters[clusters$genome != 'PCL1606', ]
for (clus in fortytwo){
  clus <- as.character(clus)
  clusters[nrow(clusters) + 1, ] = c(clus, 'PCL1606', 'type0')
  clusters[nrow(clusters) + 1, ] = c(clus, 'PCL1606', 'type1')
  clusters[nrow(clusters) + 1, ] = c(clus, 'PCL1606', 'type2')
}

# ---------------------------------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------------------------
plot1 <- ggtree(tree, branch.length = 0.8)  +
  geom_tiplab(align=TRUE, linesize=.5, size=6) + 
  geom_fruit(data=rayttree,geom=geom_tile,
             mapping = aes(y=genome, x=Color, fill=rayt),
             pwidth=0.1, offset=0.45) +
  geom_treescale()
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
    pw=0.004
  }else if(maxR==2){
    pw=0.05*maxR/3
  }else if (maxR==3) {
    pw=0.1*2/3
  }else{
    pw=0.09
  }

  plot1 <- plot1 + geom_fruit(data=dclus,geom=geom_tile,
                              mapping = aes(y=genome, x=color, fill=color),pwidth=pw,
                              offset=0.05,
                              axis.params=list(
                                axis="x",
                                title = cname,
                                title.size=6,
                                title.angle=-90,
                                title.height=0.05,
                                text.size = 0,
                                line.size = 0
                              ))
}

# Originally it was GBR xnow I am changing it to BRG
color_values <- c("type2"="green4", "type0_1"="green4", "type0_2"="green",
                  "type1"="red4", "type1_1"="red4", "type1_2"="red",
                  "type0"="steelblue", "type2_1"="steelblue", "type2_2"="blue4")

plot1 <- plot1+scale_fill_manual(values=color_values)  +
  vexpand(0.04) +
  theme(text = element_text(family = "dejavusans", size=22))


# ---------------------------------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------------------------

repcount <- data.frame(genome=unique(clusters$genome))
repcount$type0 <- rep(0, length(repcount$genome))
repcount$type1 <- rep(0, length(repcount$genome))
repcount$type2 <- rep(0, length(repcount$genome))
rownames(repcount) <- repcount$genome
for (clus in fortytwo){
  dclus <- clusters[clusters$cluster == clus,]
  for (i in 1:nrow(dclus)){
    row <- dclus[i, ]
    repcount[row$genome, row$color] = repcount[row$genome, row$color] + 1
  }
}

plot2 <- ggtree(tree)  +
  geom_tiplab(align=TRUE, linesize=.5, size=6) + 
  geom_fruit(data=rayttree,geom=geom_tile,
             mapping = aes(y=genome, x=Color, fill=rayt),
             pwidth=0.1, offset=0.9) +
  geom_treescale() + 
  scale_fill_manual(values=color_values) +
  geom_facet(mapping=aes(x=type0), fill=color_values['type0'], data=repcount, panel='Type 0',
             geom=geom_barh, stat='identity') +
  geom_facet(mapping=aes(x=type1), fill=color_values['type1'], data=repcount, panel='Type 1',
             geom=geom_barh, stat='identity') +
  geom_facet(mapping=aes(x=type2), fill=color_values['type2'], data=repcount, panel='Type 2',
             geom=geom_barh, stat='identity') +
  xlim_expand(xlim=c(0,0.2), panel='Tree') +
  xlim_expand(xlim=c(0,42), panel='Type 0') +
  xlim_expand(xlim=c(0,42), panel='Type 1') +
  xlim_expand(xlim=c(0,42), panel='Type 2') +
  theme_tree2(text = element_text(family = "dejavusans", size=22)) 


layout <- "
BBBBBBBB
#AAAAAAA
"
plot2/plot1 + plot_layout(design=layout)


