library(phyloseq)
library(ggplot2)
library(tidyr)
#install.packages("ggpubr")
#library(ggpubr)
library(ggsignif)
#library(vegan)
#library(microbiome)

#install.packages("remotes")
#remotes::install_github("microbiome/microbiome")

## Import packages
#library(phyloseq)
#library(ggplot2)
#source(file.path(params$libdir, "graphical_methods.R"))
## Alternative to source all extra function from a github repo
source("https://raw.githubusercontent.com/mahendra-mariadassou/phyloseq-extended/master/load-extra-functions.R")

setwd("/Users/annie/Documents/Projets/Yulong/metabarcoding/Yulong_R_scripts")


meta <- read.csv(file="../3_make-own-Alpha/ITS1_meta.txt","sep"="\t")
names(meta)[names(meta) == 'X'] <- 'sample'

## Setting variables
## The Phyloseq object (format rdata)
phylose <- load("../3_make-own-Alpha/ITS1.rdata")

## The experiment variable that you want to analyse
varExp <- "plot_name"

## "The alpha diversity indices to compute. Multiple indice may be indicated by separating them by a comma.
## Available indices are : Observed, Chao1, Shannon, InvSimpson, Simpson, ACE, Fisher
measures <- "InvSimpson,Observed,Shannon,"

## Create input and parameters dataframe
params <- data.frame( "phyloseq" = phylose, "measures" = measures, "varExp" = varExp)

## Load data
#load(params$phyloseq)

## Convert measures to list
measures <- as.list(strsplit(params$measures, ",")[[1]])

## Compute numeric values of alpha diversity indices
alpha.diversity <- estimate_richness(data, measures = measures)
alpha.diversity$sample <- row.names(alpha.diversity)


df <- gather(alpha.diversity, "indice", "value", -"sample")
df <- merge(df, meta)

# alpha.diversity$X <- row.names(alpha.diversity)
# alpha.diversity <- merge(alpha.diversity, meta, by="X")
# rownames(alpha.diversity) <- alpha.diversity$X

#alpha.diversity$plot_name <- factor(alpha.diversity$plot_name, levels=c("Abies","Picea","Quercus"))

#myColor=list("layer"=c("OS"="grey","OM"="beige"),
#"plot_name"=c("Abies"="#ffb3b3","Picea"="#c6ecd9","Quercus"="#66ccff"),
#"Year_Season"=c("2019_wet"="#35978f", "2020_dry"="#dfc27d","2020_wet"="navyblue"))

#data@sam_data$plot_name <- factor(data@sam_data$plot_name, levels=c("Abies","Picea","Quercus"))


df$plot_name <- factor(df$plot_name, levels=c("Abies","Picea","Quercus"))
df$indice <- factor(df$indice, levels=c("Observed","Shannon","InvSimpson"))


myColor=c("Abies"="#ffb3b3","Picea"="#c6ecd9","Quercus"="#66ccff")
myColor=c("Abies"="red","Picea"="forestGreen","Quercus"="blue")

ggplot(data=df, aes(x=plot_name, y=value, colour=plot_name) ) +
  #geom_col() + 
  geom_boxplot(outlier.shape = NA)+
  geom_point(alpha=0.7, size=0.01)+
  theme_bw()+
  #geom_jitter(cex=0.2)+
  facet_wrap(~ indice, scales = "free") + 
  #scale_fill_manual(values=c("#999999",rev(c("#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7","springgreen4","mediumorchid4","cadetblue","lightpink","khaki","turquoise1","black","steelblue")))) +
  theme(axis.text.x=element_text(angle = 55,hjust=1))+
  geom_signif(
    comparisons = list(c("Abies", "Picea"), c( "Quercus", "Picea"),
                       c("Abies", "Quercus")),
    step_increase = 0.1,
    vjust = 1.5,
    map_signif_level = TRUE,
    size = 0.2,
    test = "wilcox.test"
  )+
  scale_color_manual(values=myColor)+
  scale_shape_manual(values=c(1,1,1))+
  theme(plot.title = element_text(hjust = 0.5), axis.title.x=element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none"
  ) 

#  size = [389 - 276]   SVG
# size = [3 - 3.5] PDF
# 
# 
# 
# p <- plot_richness(data, x = params$varExp, color = params$varExp, measures = measures) +
#   #geom_point(size=0.05, alpha=0.7)+
#   #ggtitle(paste("Alpha diversity distribution in function of tree type"))+ 
#   theme_bw()+
#   # using `ggsignif` to display comparison of interest
#   geom_signif(
#     comparisons = list(c("Abies", "Picea"), c( "Quercus", "Picea"),
#                        c("Abies", "Quercus")),
#     step_increase = 0.1,
#     vjust = 1.5,
#     map_signif_level = TRUE,
#     size = 0.2,
#     test = "wilcox.test"
#   )+
#   scale_color_manual(values=myColor)+
#   scale_shape_manual(values=c(1,1,1))+
#   theme(plot.title = element_text(hjust = 0.5), axis.title.x=element_blank()
#         )
# plot(p)
# 
# 
# p <- p + geom_boxplot(alpha = 0.2) +
#   geom_point()+ #theme_grey() +
#   #theme(axis.text.x = element_text(angle=90, hjust=1)) +
#   theme(plot.title = element_text(hjust = 0.5)) 
# plot(p)
# 



anova_data <- cbind(sample_data(data), alpha.diversity)
anova_data$Depth <- sample_sums(data)

variables <- paste(sep=" + ", "Depth", params$varExp )

## Perform ANOVA on observed richness, which effects are significant
for (m in measures){
  f <- paste(m," ~ ", variables)
  cat(sep = "", "###############################################################\n#Perform ANOVA on ",m,", which effects are significant\nanova.",m," <-aov( ",f,", anova_data)\nsummary(anova.",m,")\n")
  anova_res <- aov( as.formula(f), anova_data)
  res <- summary(anova_res)
  print(res)
  cat("\n\n")
}




