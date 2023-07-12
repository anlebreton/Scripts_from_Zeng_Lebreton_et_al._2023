library(phyloseq)
library(ggplot2)
library(tidyr)
library(dplyr)
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


meta <- read.csv(file="../0_Metadata/Yulong_sampling_bacteria_metadata.txt","sep"="\t")
meta <- read.csv(file="../0_Metadata/Yulong_sampling_ITS2_metadata.txt","sep"="\t")

names(meta)[names(meta) == 'X'] <- 'sample'
meta$sample <- gsub("-", ".", meta$sample)

## Setting variables
## The Phyloseq object (format rdata)
phylose <- load("../1_data/Yulong_bacteria.rdata")
phylose <- load("../1_data/Yulong_ITS2.rdata")

## The experiment variable that you want to analyse
varExp <- "Tree_species"

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
df <- df %>% mutate(Layer=case_when(Layer != "OS" ~ Layer, TRUE ~ "OS: 0 - 5cm")) 
df <- df %>% mutate(Layer=case_when(Layer != "OM" ~ Layer, TRUE ~ "OM: 5 - 25cm"))
df$Layer <- factor(df$Layer, levels=c("OS: 0 - 5cm","OM: 5 - 25cm" ))

# alpha.diversity$X <- row.names(alpha.diversity)
# alpha.diversity <- merge(alpha.diversity, meta, by="X")
# rownames(alpha.diversity) <- alpha.diversity$X

#alpha.diversity$Tree_species <- factor(alpha.diversity$Tree_species, levels=c("Abies","Picea","Quercus"))

#myColor=list("layer"=c("OS"="grey","OM"="beige"),
#"Tree_species"=c("Abies"="#ffb3b3","Picea"="#c6ecd9","Quercus"="#66ccff"),
#"Year_Season"=c("2019_wet"="#35978f", "2020_dry"="#dfc27d","2020_wet"="navyblue"))

#data@sam_data$Tree_species <- factor(data@sam_data$Tree_species, levels=c("Abies","Picea","Quercus"))


df$Tree_species <- factor(df$Tree_species, levels=c("Abies","Picea","Quercus"))
df$indice <- factor(df$indice, levels=c("Observed","Shannon","InvSimpson"))
#df$Layer <- factor(df$Layer, levels=c("OS","OM"))
df$Season <- factor(df$Season, levels=c("2019_wet", "2020_dry", "2020_wet"))



#####
df <- df %>% filter(indice == "Observed")
# plot_max = 5000
# red = 400
# red2 = 150
# plot_min = 2000
 plot_max = 800
 plot_min = 100
 red = 100
 red2=50
#####

myColor=c("Abies"="#ffb3b3","Picea"="#c6ecd9","Quercus"="#66ccff")
myColor=c("Abies"="red","Picea"="forestGreen","Quercus"="blue")
myColor=c("Abies"="#d64933","Picea"="#3fb52f","Quercus"="#5a92ed")

## plot(p1)

p1 <- ggplot(data=df, aes(x=Tree_species, y=value) ) +
  #geom_col() + 
  ggtitle("Tree species")+
  ylim(plot_min, plot_max)+
  geom_boxplot(outlier.shape = NA, aes( colour= Tree_species))+
  geom_point(alpha=0.7, size=0.01, aes( colour= Tree_species))+
  #ggtitle("Tree species")+
  theme_bw()+
  #geom_jitter(cex=0.2)+
  #facet_wrap(~ indice, scales = "free") + 
  #scale_fill_manual(values=c("#999999",rev(c("#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7","springgreen4","mediumorchid4","cadetblue","lightpink","khaki","turquoise1","black","steelblue")))) +
  theme(axis.text.x=element_text(angle = 55,hjust=1))+
  geom_signif(
    comparisons = list(c("Abies", "Picea"), c( "Quercus", "Picea"),
                       c("Abies", "Quercus")),
    step_increase = 0.05,
    y_position = plot_max - red,
    vjust = 1.5,
    tip_length =0,
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


myColor=c("OS: 0 - 5cm"="#6B3A2C","OM: 5 - 25cm"="#dbcb81")

p2 <- ggplot(data=df, aes(x=Layer, y=value) ) +
  #geom_col() + 
  ggtitle("Layer")+
  ylim(plot_min, plot_max)+
  geom_boxplot(outlier.shape = NA, aes( colour=Layer))+
  geom_point(alpha=0.7, size=0.01, aes( colour=Layer))+
  theme_bw()+
  #geom_jitter(cex=0.2)+
  #facet_wrap(~ indice, scales = "free") + 
  #scale_fill_manual(values=c("#999999",rev(c("#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7","springgreen4","mediumorchid4","cadetblue","lightpink","khaki","turquoise1","black","steelblue")))) +
  theme(axis.text.x=element_text(angle = 55,hjust=1))+
  geom_signif(
    comparisons = list(c("OS: 0 - 5cm","OM: 5 - 25cm")),
    step_increase = 0.05,
    y_position = plot_max - red2,
    vjust = 1.5,
    tip_length =0,
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


 



myColor=c("2019_wet"="#35978f", "2020_dry"="#dfc27d","2020_wet"="navyblue")

p3 <- ggplot(data=df, aes(x=Season, y=value) ) +
  ggtitle("Season")+
  ylim(plot_min, plot_max)+
  #geom_col() + 
  geom_boxplot(outlier.shape = NA, aes( colour=Season))+
  geom_point(alpha=0.7, size=0.01, aes( colour=Season))+
  theme_bw()+
  #geom_jitter(cex=0.2)+
  #facet_wrap(~ indice, scales = "free") + 
  #scale_fill_manual(values=c("#999999",rev(c("#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7","springgreen4","mediumorchid4","cadetblue","lightpink","khaki","turquoise1","black","steelblue")))) +
  theme(axis.text.x=element_text(angle = 55,hjust=1))+
  geom_signif(
    comparisons = list(c("2019_wet", "2020_dry"),c("2019_wet", "2020_wet"),c("2020_dry", "2020_wet")),
    step_increase = 0.05,
    tip_length =0,
    vjust = 1.5,
    map_signif_level = TRUE,
    y_position = plot_max - red,
    size = 0.2,
    test = "wilcox.test"
  )+
  scale_color_manual(values=myColor)+
  scale_shape_manual(values=c(1,1,1))+
  theme(plot.title = element_text(hjust = 0.5), axis.title.x=element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none"
  ) 



library(ggplotify)
library(ggpubr)
ggarrange(p1, p3, p2, #+ rremove("x.text"), 
          #labels = c("A", "B", "C","D"),
          align = "h",
          ncol = 3, nrow = 1,
          common.legend = FALSE#, 
          #legend=NA
          ) 
#(4.5*4.5)

 
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

# 
# DF <- df[df$indice == "Observed",]
# 
# a <- DF %>% select("sample","value","Season")
# colnames(a) <- c("sample","value","condition")
# a$ID <- paste(a$sample,"_1", sep="")
# a$variable <- "Season"
# b <- DF %>% select("sample","value","Layer")
# colnames(b) <- c("sample","value","condition")
# b$ID <- paste(a$sample,"_2", sep="")
# b$variable <- "Layer"
# c <- DF %>% select("sample","value","Tree_species")
# colnames(c) <- c("sample","value","condition")
# c$ID <- paste(a$sample,"_3", sep="")
# c$variable <- "Tree_species"
# 
# DF <- rbind(c,b,a)
# 
# ColorCondition <- c("2019_wet"="#35978f", "2020_dry"="#dfc27d","2020_wet"="navyblue", 
#                     "OS"="#6B3A2C","OM"="#dbcb81", "Abies"="#d64933","Picea"="#3fb52f","Quercus"="#5a92ed") 
# 
# DF$variable <-factor(DF$variable, levels= c("Tree_species", "Season", "Layer")) 
# DF$condition <- factor(DF$condition, levels=c("Abies","Picea","Quercus","2019_wet","2020_dry","2020_wet", "OS","OM"))
# 
# ggplot(data=DF, aes(x=condition, y=value, colour=condition) ) +
#   #geom_col() + 
#   geom_boxplot(outlier.shape = NA)+
#   geom_point(alpha=0.7, size=0.01)+
#   theme_bw()+
#   #geom_jitter(cex=0.2)+
#   facet_wrap(~ variable, scales = "free") + 
#   #scale_fill_manual(values=c("#999999",rev(c("#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7","springgreen4","mediumorchid4","cadetblue","lightpink","khaki","turquoise1","black","steelblue")))) +
#   theme(axis.text.x=element_text(angle = 55,hjust=1))+
#   geom_signif(
#     comparisons = list(c("2019_wet", "2020_dry"),c("2019_wet", "2020_wet"),c("2020_dry", "2020_wet"),
#                        c("Abies", "Picea"), c( "Quercus", "Picea"), c("Abies", "Quercus") ,
#                        c("OS","OM")),
#     step_increase = 0.1,
#     vjust = 1.5,
#     map_signif_level = TRUE,
#     size = 0.2,
#     test = "wilcox.test"
#   )+
#   scale_color_manual(values=ColorCondition)+
#   scale_shape_manual(values=c(1,1,1))+
#   theme(plot.title = element_text(hjust = 0.5), axis.title.x=element_blank(),
#         axis.title.y = element_blank(),
#         legend.position = "none"
#   ) 
# 
# 
# 


