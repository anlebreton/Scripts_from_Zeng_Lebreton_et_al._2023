## Import packages
library(phyloseq)
library(ggplot2)
library(ggpubr)
library(plotly)

##  source all extra function from a github repo
 source("https://raw.githubusercontent.com/mahendra-mariadassou/phyloseq-extended/master/load-extra-functions.R")

## Setting variables
## The Phyloseq object (format rdata)
phyloseq <- load("../1_data/Yulong_ITS2.rdata")



## The experiment variable that you want to analyse
varExp <-  "Tree_species"

## The methods of beta diversity you want to compute
## to see all available distance methods, type distanceMethodList
## The most common one are : cc for Jaccard, bray for Bray-Curtis, unifrac and wunifrac for Unifrac and weighted Unifrac
## N.B. if the tree is not available in your RData, you cannot choose Unifrac or Weighted Unifrac
## You may precise multiple distance by separating them by a comma
 methods <- "cc,bray,unifrac,wunifrac"

## Create input and parameters dataframe
params <- data.frame( "phyloseq" = phyloseq, "varExp" = varExp, "methods" = methods)

## store methods in list
methods <- as.list(strsplit(params$methods, ",")[[1]])

## Order samples according to grouping variable
sampleOrder <- levels(reorder(sample_names(data), as.numeric(get_variable(data, params$varExp)))) 

# Heatmap plot of the beta distance 
for (method in methods){
  dist.a <- distance(data, method = method)
  a <- as.matrix(dist.a)
  write.table(a, paste(sep="", method, ".tsv"), sep="\t", quote=FALSE, col.names=NA)
  pa <- plot_dist_as_heatmap(dist.a, order = sampleOrder, title = paste("Heatmap plot of the beta distance :",method)) +
    theme(plot.title = element_text(hjust = 0.5))
  plot(pa)
}


#-------------------------------------------------------------
#____________ TREE association ____________________________
#-------------------------------------------------------------

## The beta diversity distance matrix file
 distance <- "cc.tsv"

## The ordination method you want to use
## You can choose between "MDS" (for MDS/PCoA), "NMDS" or "DPCoA"
 method <- "NMDS"

## The experiment variable that you want to analyse
 varExp <- "Tree_species"

## Create input and parameters dataframe
 params <- data.frame( "phyloseq" = phyloseq, "distance" = distance, "method" = method, "varExp" = varExp)

## the distance matrix file
A    <- read.table(file=params$distance, row.names=1)
dist <- as.dist(A)

## Add SampleID variable to physeq metadata
sample_data(data)$SampleID <- sample_names(data)

ord  <- ordinate(data, method = params$method, distance = dist)
var  <- get_variable(data, params$varExp)

data@sam_data$Tree_species <- factor(data@sam_data$Tree_species, levels=c("Abies","Picea","Quercus"))

myColor=c("Abies"="#d64933","Picea"="#3fb52f","Quercus"="#5a92ed")


## NMDS 
p1t <- plot_ordination(data, ord, color = params$varExp) + 
  geom_text(aes(label = SampleID), alpha = 0) + ## add invisible labels
  theme_bw() + ggtitle(paste(params$method,"Jaccard")) +
  scale_color_manual(values=myColor)+
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.x=element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none")  
ggplotly(p1t, tooltip = c("colour", "label"))



distance <- "bray.tsv"
params <- data.frame( "phyloseq" = phyloseq, "distance" = distance, "method" = method, "varExp" = varExp)
A    <- read.table(file=params$distance, row.names=1)
dist <- as.dist(A)
ord  <- ordinate(data, method = params$method, distance = dist)
var  <- get_variable(data, params$varExp)

p2t <- plot_ordination(data, ord, color = params$varExp) + 
  geom_text(aes(label = SampleID), alpha = 0) + ## add invisible labels
  theme_bw() +# ggtitle(paste(params$method,"Bray-Curtis")) +
  geom_point(size = 1)+
  ggtitle(" ") +
  scale_color_manual(values=myColor)+
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.x=element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none"
        )  +
  stat_ellipse( type = "t", linetype = 1, level = 0.95) 

ggplotly(p2t, tooltip = c("colour", "label"))
  


distance <- "unifrac.tsv"
params <- data.frame( "phyloseq" = phyloseq, "distance" = distance, "method" = method, "varExp" = varExp)
A    <- read.table(file=params$distance, row.names=1)
dist <- as.dist(A)
ord  <- ordinate(data, method = params$method, distance = dist)
var  <- get_variable(data, params$varExp)

p3t <- plot_ordination(data, ord, color = params$varExp) + 
  geom_text(aes(label = SampleID), alpha = 0) + ## add invisible labels
  theme_bw() + ggtitle(paste(params$method,"Unifrac")) +
  scale_color_manual(values=myColor)+
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.x=element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none"
  )  
ggplotly(p3t, tooltip = c("colour", "label"))

distance <- "wunifrac.tsv"
params <- data.frame( "phyloseq" = phyloseq, "distance" = distance, "method" = method, "varExp" = varExp)
A    <- read.table(file=params$distance, row.names=1)
dist <- as.dist(A)
ord  <- ordinate(data, method = params$method, distance = dist)
var  <- get_variable(data, params$varExp)

p4t <- plot_ordination(data, ord, color = params$varExp) + 
  geom_text(aes(label = SampleID), alpha = 0) + ## add invisible labels
  theme_bw() + ggtitle(paste(params$method,"WUnifrac")) +
  scale_color_manual(values=myColor)+
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.x=element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none"
  )  
ggplotly(p4t, tooltip = c("colour", "label"))




ggarrange(p1t, p2t, p3t, p4t, #+ rremove("x.text"), 
          labels = c("A", "B", "C","D"),
          ncol = 2, nrow = 2,
          common.legend = TRUE, legend="bottom") 
# [5,17]
#size= [6,8] PDF
#size = [800 - 575]


#-------------------------------------------------------------
#__________________ soil Layer    ____________________________
#-------------------------------------------------------------

## The beta diversity distance matrix file
distance <- "cc.tsv"

## The ordination method you want to use
## You can choose between "MDS" (for MDS/PCoA), "NMDS" or "DPCoA"
method <- "NMDS"

## The experiment variable that you want to analyse
varExp <- "Layer"

## Create input and parameters dataframe
params <- data.frame( "phyloseq" = phyloseq, "distance" = distance, "method" = method, "varExp" = varExp)

## Load data
## the phyloseq object
#load(params$phyloseq)

## the distance matrix file
A    <- read.table(file=params$distance, row.names=1)
dist <- as.dist(A)

## Add SampleID variable to physeq metadata
sample_data(data)$SampleID <- sample_names(data)

ord  <- ordinate(data, method = params$method, distance = dist)
var  <- get_variable(data, params$varExp)

data@sam_data$Layer <- factor(data@sam_data$Layer, levels=c("OM","OS"))

myColor=c("OS"="#6B3A2C","OM"="#dbcb81")



## plot(p1)
p1l <- plot_ordination(data, ord, color = params$varExp) + 
  geom_text(aes(label = SampleID), alpha = 0) + ## add invisible labels
  theme_bw() + ggtitle(paste(params$method,"Jaccard")) +
  scale_color_manual(values=myColor)+
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.x=element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none")  
ggplotly(p1l, tooltip = c("colour", "label"))



distance <- "bray.tsv"
params <- data.frame( "phyloseq" = phyloseq, "distance" = distance, "method" = method, "varExp" = varExp)
A    <- read.table(file=params$distance, row.names=1)
dist <- as.dist(A)
ord  <- ordinate(data, method = params$method, distance = dist)
var  <- get_variable(data, params$varExp)

p2l <- plot_ordination(data, ord, color = params$varExp) + 
  geom_text(aes(label = SampleID), alpha = 0) +## add invisible labels
  stat_ellipse( type = "t", linetype = 1, level = 0.95) + 
  geom_point(size = 1)+
  theme_bw() + # ggtitle(paste(params$method,"Bray-Curtis")) +
  ggtitle(" ") +
  scale_color_manual(values=myColor)+
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.x=element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none"
  )  
p2l
ggplotly(p2l, tooltip = c("colour", "label"))


distance <- "unifrac.tsv"
params <- data.frame( "phyloseq" = phyloseq, "distance" = distance, "method" = method, "varExp" = varExp)
A    <- read.table(file=params$distance, row.names=1)
dist <- as.dist(A)
ord  <- ordinate(data, method = params$method, distance = dist)
var  <- get_variable(data, params$varExp)

p3l <- plot_ordination(data, ord, color = params$varExp) + 
  geom_text(aes(label = SampleID), alpha = 0) + ## add invisible labels
  theme_bw() + ggtitle(paste(params$method,"Unifrac")) +
  scale_color_manual(values=myColor)+
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.x=element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none"
  )  
ggplotly(p3l, tooltip = c("colour", "label"))

distance <- "wunifrac.tsv"
params <- data.frame( "phyloseq" = phyloseq, "distance" = distance, "method" = method, "varExp" = varExp)
A    <- read.table(file=params$distance, row.names=1)
dist <- as.dist(A)
ord  <- ordinate(data, method = params$method, distance = dist)
var  <- get_variable(data, params$varExp)

p4l <- plot_ordination(data, ord, color = params$varExp) + 
  geom_text(aes(label = SampleID), alpha = 0) + ## add invisible labels
  theme_bw() + ggtitle(paste(params$method,"WUnifrac")) +
  scale_color_manual(values=myColor)+
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.x=element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none"
  )  
ggplotly(p4l, tooltip = c("colour", "label"))




ggarrange(p1l, p2l, p3l, p4l, #+ rremove("x.text"), 
          labels = c("A", "B", "C","D"),
          ncol = 2, nrow = 2,
          common.legend = TRUE, legend="bottom") 

#size= [6,8] PDF
#size = [800 - 575]



#-------------------------------------------------------------
#__________________ Season   ____________________________
#-------------------------------------------------------------


## The beta diversity distance matrix file
distance <- "cc.tsv"

## The ordination method you want to use
## You can choose between "MDS" (for MDS/PCoA), "NMDS" or "DPCoA"
method <- "NMDS"

## The experiment variable that you want to analyse
varExp <- "Season"

## Create input and parameters dataframe
params <- data.frame( "phyloseq" = phyloseq, "distance" = distance, "method" = method, "varExp" = varExp)

## Load data
## the phyloseq object
#load(params$phyloseq)

## the distance matrix file
A    <- read.table(file=params$distance, row.names=1)
dist <- as.dist(A)

## Add SampleID variable to physeq metadata
sample_data(data)$SampleID <- sample_names(data)

ord  <- ordinate(data, method = params$method, distance = dist)
var  <- get_variable(data, params$varExp)

data@sam_data$Season <- factor(data@sam_data$Season, levels=c("2019_wet", "2020_dry", "2020_wet"))
myColor=c("2019_wet"="#35978f", "2020_dry"="#dfc27d","2020_wet"="navyblue")


## plot(p1)
p1s <- plot_ordination(data, ord, color = params$varExp) + 
  geom_text(aes(label = SampleID), alpha = 0) + ## add invisible labels
  theme_bw() + ggtitle(paste(params$method,"Jaccard")) +
  scale_color_manual(values=myColor)+
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.x=element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none")  
ggplotly(p1s, tooltip = c("colour", "label"))



distance <- "bray.tsv"
params <- data.frame( "phyloseq" = phyloseq, "distance" = distance, "method" = method, "varExp" = varExp)
A    <- read.table(file=params$distance, row.names=1)
dist <- as.dist(A)
ord  <- ordinate(data, method = params$method, distance = dist)
var  <- get_variable(data, params$varExp)

p2s <- plot_ordination(data, ord, color = params$varExp) + 
  geom_text(aes(label = SampleID), alpha = 0) + ## add invisible labels
  stat_ellipse( type = "t", linetype = 1, level = 0.95) +
  theme_bw() + #ggtitle(paste(params$method,"Bray-Curtis")) +
  ggtitle(" ") +
  geom_point(size = 1)+
  scale_color_manual(values=myColor)+
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.x=element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none"
        
  )  
ggplotly(p2s, tooltip = c("colour", "label"))


distance <- "unifrac.tsv"
params <- data.frame( "phyloseq" = phyloseq, "distance" = distance, "method" = method, "varExp" = varExp)
A    <- read.table(file=params$distance, row.names=1)
dist <- as.dist(A)
ord  <- ordinate(data, method = params$method, distance = dist)
var  <- get_variable(data, params$varExp)

p3s <- plot_ordination(data, ord, color = params$varExp) + 
  geom_text(aes(label = SampleID), alpha = 0) + ## add invisible labels
  theme_bw() + ggtitle(paste(params$method,"Unifrac")) +
  scale_color_manual(values=myColor)+
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.x=element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none"
  )  
ggplotly(p3s, tooltip = c("colour", "label"))

distance <- "wunifrac.tsv"
params <- data.frame( "phyloseq" = phyloseq, "distance" = distance, "method" = method, "varExp" = varExp)
A    <- read.table(file=params$distance, row.names=1)
dist <- as.dist(A)
ord  <- ordinate(data, method = params$method, distance = dist)
var  <- get_variable(data, params$varExp)

p4s <- plot_ordination(data, ord, color = params$varExp) + 
  geom_text(aes(label = SampleID), alpha = 0) + ## add invisible labels
  theme_bw() + ggtitle(paste(params$method,"WUnifrac")) +
  scale_color_manual(values=myColor)+
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.x=element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none"
  )  
ggplotly(p4s, tooltip = c("colour", "label"))




ggarrange(p1s, p2s, p3s, p4s, #+ rremove("x.text"), 
          labels = c("A", "B", "C","D"),
          #ncol = 2, nrow = 2,
          ncol = 2, nrow = 2,
          common.legend = TRUE, legend="bottom") 

#size= [6,8] PDF
#size = [800 - 575]




#-------------------------------------------------------------
#__________________ combined plot   _________________________
#-------------------------------------------------------------

ggarrange(p2t, p2s, p2l, #+ rremove("x.text"), 
          labels = c("Tree species", "Season", "Layer"),
          font.label = list(size = 12, color = "grey20"  , face = "plain",family = NULL) ,
          ncol = 3, nrow = 1,
          #common.legend = TRUE,
          legend="bottom") 

#12*4

