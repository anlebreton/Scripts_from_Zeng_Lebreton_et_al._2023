library(phyloseq)
library(ggplot2)
library(tidyr)
library(dplyr)
library(ggsignif)
library(ggplotify)
library(ggpubr)

meta <- read.csv(file="../0_Metadata/Yulong_sampling_ITS2_metadata.txt","sep"="\t")

names(meta)[names(meta) == 'X'] <- 'sample'
meta$sample <- gsub("-", ".", meta$sample)

## Setting variables
## The Phyloseq object (format rdata)
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


df$Tree_species <- factor(df$Tree_species, levels=c("Abies","Picea","Quercus"))
df$indice <- factor(df$indice, levels=c("Observed","Shannon","InvSimpson"))
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

###### ANOVA #########
 
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
 
############### 
 
 
 ### Tree association part of the graph
 
 
myColor=c("Abies"="#d64933","Picea"="#3fb52f","Quercus"="#5a92ed")

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

p1
#  size = [389 - 276]   SVG
# size = [3 - 3.5] PDF

### Layer part of the graph

myColor=c("OS: 0 - 5cm"="#6B3A2C","OM: 5 - 25cm"="#dbcb81")

p2 <- ggplot(data=df, aes(x=Layer, y=value) ) +
  ggtitle("Layer")+
  ylim(plot_min, plot_max)+
  geom_boxplot(outlier.shape = NA, aes( colour=Layer))+
  geom_point(alpha=0.7, size=0.01, aes( colour=Layer))+
  theme_bw()+
  #facet_wrap(~ indice, scales = "free") + 
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


### seasons part of the graph

myColor=c("2019_wet"="#35978f", "2020_dry"="#dfc27d","2020_wet"="navyblue")

p3 <- ggplot(data=df, aes(x=Season, y=value) ) +
  ggtitle("Season")+
  ylim(plot_min, plot_max)+
  geom_boxplot(outlier.shape = NA, aes( colour=Season))+
  geom_point(alpha=0.7, size=0.01, aes( colour=Season))+
  theme_bw()+
  #facet_wrap(~ indice, scales = "free") + 
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


### create a figure with the 3 subfigures. 
ggarrange(p1, p3, p2,
          align = "h",
          ncol = 3, nrow = 1,
          common.legend = FALSE
          ) 
#(4.5*4.5)



