if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
library(tidyverse)
library(ggtree)
library(glue)
library(treeio)
library(ggstance)
library(hash)

# Where are the input files stored?
TREEPATH <- "input/tree.nwk"
CLUSPATH <- "temp/cluster_output_Sep20_875/clusters_Sep20.txt"
GENORDER <- c('ChPhzTR38', 'ChPhzS24', 'ChPhzTR18', 'ChPhzTR39', 'PA23', 'ChPhzS23', '66',
              'Lzh-T5', 'O6', 'ChPhzTR36', '6698', 'C50', 'P2', '2210', '19603', 'CW2',
              'Q16', '464', '449', 'M12', 'StFRB508', 'JD37', 'M71', 'K27', 'Pb-St2',
              'B25', 'SLPH10', 'ChPhzS135', 'DTR133', 'ToZa7', 'ChPhzTR44', 'PCL1607',
              'PCL1391', 'ChPhzS140', '17809', '17411', 'ZJU60', '21509', '189', '17415',
              '50083', 'TAMOak81')
colorkey <- hash()
colorkey[['type0']] <- 'red'
colorkey[['type1']] <- 'blue'
colorkey[['type2']] <- 'green'

# Input data
tree <- read.tree(TREEPATH)
clus <- read.table(file = CLUSPATH, sep = ' ', 
                   header = FALSE)
colnames(clus) <- c("cluster_number", "genome", "repin_start",
                    "repin_end", "type", "seq")
clus$genome <- gsub("chl", "", as.character(clus$genome))
nodekey <- data.frame(name=GENORDER)
nodekey <- summarize(group_by(nodekey, name),
                     number=which(tree$tip.label== name))

# Which cluster do you want to see?
whichclus <- 982
thisclus <- clus[clus$cluster_number == whichclus, ]
checklist <- data.frame(genome=thisclus$genome, type=thisclus$type)
checklist <- summarize(group_by(checklist, genome, type),
          gnumber=which(tree$tip.label== genome))
checklist$gnumber <- as.numeric(checklist$gnumber)

# Basic Tree
xv <- 50
plot1<- ggtree(tree, branch.length = 'none')  +
  geom_tiplab(align=TRUE, linesize=.5) + xlim(0, xv) +
  geom_cladelabel(label="", node=43, color="black", 
                  offset=xv * 0.065, align=TRUE, barsize = 1.5) +
  geom_cladelabel(label="", node=47, color="blue", 
                  offset=xv * 0.065, align=TRUE, barsize = 1.5) +
  geom_cladelabel(label="", node=58, color="red", 
                  offset=xv * 0.065, align=TRUE, barsize = 1.5) +
  geom_cladelabel(label="", node=70, color="green", 
                  offset=xv * 0.065, align=TRUE, barsize = 1.5)


for (i in 1:nrow(checklist)){
  row <- checklist[i, ]
  plot1 <- plot1 + geom_cladelabel(label=row$type, node=row$gnumber, 
                                   color=colorkey[[row$type]], 
                                   offset= xv * 0.07, align=TRUE, barsize = 1.5)
}
# Note the offsets work when you zoom on the graph/open it via the zoom option in RStudio
plot1 + xlab(glue("Cluster {whichclus}"))

# --------------------------------------------------------------------------------


# # --------------------------------------------------------------------------------
# # Which cluster do you want to see?
# whichclus <- 982
# thisclus <- clus[clus$cluster_number == whichclus, ]
# checklist <- data.frame('genome' = GENORDER,
#                         'present' = rep(0, length(GENORDER)))
# checklist$present <- as.numeric(checklist$genome %in% thisclus$genome)
# checklist[checklist == 0] <- NA
# # checklist$genome <- as.character(
# #   factor(
# #     checklist$genome, 
# #     levels = nodekey$name, 
# #     labels = nodekey$number
# #   )
# # )
# 
# # Basic Tree
# xv <- 20
# plot1 <- ggtree(tree, branch.length = 'none')  +
#   geom_tiplab(align=TRUE, linesize=.5) + xlim(0, xv) +
#   geom_cladelabel(label="", node=43, color="black", offset=xv * 0.12, align=TRUE, barsize = 1.5) +
#   geom_cladelabel(label="", node=47, color="blue", offset=xv * 0.12, align=TRUE, barsize = 1.5) +
#   geom_cladelabel(label="", node=58, color="red", offset=xv * 0.12, align=TRUE, barsize = 1.5) +
#   geom_cladelabel(label="", node=70, color="green", offset=xv * 0.12, align=TRUE, barsize = 1.5)
# plot1
# 
# # Plotting the REPINs on a separate graph
# facet_plot(plot1, panel = "Clusters", data=checklist, geom=geom_point,
#              aes(x=present))
# 
# # --------------------------------------------------------------------------------