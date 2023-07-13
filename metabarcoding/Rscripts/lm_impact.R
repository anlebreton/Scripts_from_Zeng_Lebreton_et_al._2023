library(phyloseq)
library(ggplot2)
library(tidyr)
library(dplyr)
library(ggsignif)
library(vegan)


GetMyColor <- function(colDF){
  mycolorsCAZymes <- read.csv("../../../mycolorFile.txt",sep="\t")
  mycolors <- mycolorsCAZymes %>% filter(Name %in% colDF)
  mycolors <- mycolors %>% arrange(myLevel)
  mycolor <- mycolors$color
  names(mycolor) <- mycolors$Name
  #colDF <- factor(colDF, levels=mycolors$Name
  return(mycolor)
}

meta <- read.csv(file="../../0_Metadata/Yulong_sampling_ITS2_metadata.txt","sep"="\t")

names(meta)[names(meta) == 'X'] <- 'Sample'
meta$sample <- gsub("-", ".", meta$Sample)
meta$Tree_layer <-  paste(meta$Tree_species,meta$Layer,sep="_")


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

## Convert measures to list
measures <- as.list(strsplit(params$measures, ",")[[1]])

## Compute numeric values of alpha diversity indices
alpha.diversity <- estimate_richness(data, measures = measures)
alpha.diversity$sample <- row.names(alpha.diversity)


df <- gather(alpha.diversity, "indice", "value", -"sample")
df <- merge(df, meta)
#fix names
df <- df %>% mutate(Layer=case_when(Layer != "OS" ~ Layer, TRUE ~ "OS: 0 - 5cm")) 
df <- df %>% mutate(Layer=case_when(Layer != "OM" ~ Layer, TRUE ~ "OM: 5 - 25cm"))
df$Layer <- factor(df$Layer, levels=c("OS: 0 - 5cm","OM: 5 - 25cm" ))

#fix factor level
df$Tree_species <- factor(df$Tree_species, levels=c("Abies","Picea","Quercus"))
df$indice <- factor(df$indice, levels=c("Observed","Shannon","InvSimpson"))
df$Season <- factor(df$Season, levels=c("2019_wet", "2020_dry", "2020_wet"))


DF <- df %>% filter(indice == "Observed")
DF <- DF %>% rename("Corg/N"=ratioCorg_N)

# Impact of the soil caracteristic ?
summary(lm(value ~ Org, DF))
summary(lm(value ~ N, DF))
summary(lm(value ~ log(`Corg/N`), DF))
summary(lm(value ~ pH, DF))
summary(lm(value ~ log(Ca), DF))
summary(lm(value ~ CEC, DF))
summary(lm(value ~ P, DF))
summary(lm(value ~ DBH, DF))


tmp <- !is.na(DF)
# correlation among variables ?
cor(tmp[, c('pH','N','CEC','P','Ca','Corg', "Corg/N" , "DBH")])

myColor <- GetMyColor(DF$Tree_species)

# Plots 
p1 <- ggplot(DF, aes(x = Org , y = value, colour = Tree_species)) +
  geom_point(size = 1) +
  scale_colour_manual(values=myColor)+
  geom_smooth(method = "lm", colour = "black", fill = "grey85", linewidth=0.5) +
  xlab(bquote("SOM" ~ (g/kg^-1))) + 
  labs(y=NULL) + 

  theme_bw()
p1

p2 <- ggplot(DF, aes(x = N, y = value, colour = Tree_species)) +
  geom_point(size = 1) +
  scale_colour_manual(values=myColor)+
  geom_smooth(method = "lm", colour = "black", fill = "grey85", linewidth=0.5) +
  #scale_colour_viridis_d(option = "magma", begin = 0.2, end = 0.8) +
  xlab(bquote("Nitrogen" ~ (g/kg^-1))) + 
  #ylab("fungal alpha diversity") +
  labs(y=NULL) + #ylim(1.1,3) +
  theme_bw()

p2


p3 <- ggplot(DF, aes(x = log(`Corg/N`), y = value, colour = Tree_species)) +
  geom_point(size = 1) +
  scale_colour_manual(values=myColor)+
  #geom_smooth(method = "lm", colour = "black", fill = "grey85", linewidth=0.5) +
  #scale_colour_viridis_d(option = "magma", begin = 0.2, end = 0.8) +
  xlab(bquote("log(Corg/N)" )) + 
  #ylab("fungal alpha diversity") +
  labs(y=NULL) + #ylim(1.1,3) +
  theme_bw()

p3


p4 <- ggplot(DF, aes(x = pH, y = value, colour = Tree_species)) +
  geom_point(size = 1) +
  scale_colour_manual(values=myColor)+
  geom_smooth(method = "lm", colour = "black", fill = "grey85", linewidth=0.5) +
  #scale_colour_viridis_d(option = "magma", begin = 0.2, end = 0.8) +
  xlab(bquote("pH" )) +  labs(y=NULL) + 
  #ylab("fungal alpha diversity") +
  #ylim(1.1,3) +
  theme_bw()
p4


p5 <- ggplot(DF, aes(x = log(Ca), y = value, colour = Tree_species)) +
  geom_point(size = 1) +
  scale_colour_manual(values=myColor)+
  geom_smooth(method = "lm", colour = "black", fill = "grey85", linewidth=0.5) +
  #scale_colour_viridis_d(option = "magma", begin = 0.2, end = 0.8) +
  xlab(bquote("Calcium" ~ log(mg/kg^-1))) + 
  labs(y=NULL) + 
  #ylab("fungal alpha diversity") +
  #ylim(1.1,3) +
  theme_bw()
p5



p6 <- ggplot(DF, aes(x = CEC , y = value, colour = Tree_species)) +
  geom_point(size = 1) +
  scale_colour_manual(values=myColor)+
  geom_smooth(method = "lm", colour = "black", fill = "grey85", linewidth=0.5) +
  #scale_colour_viridis_d(option = "magma", begin = 0.2, end = 0.8) +
  xlab(bquote("CEC" ~ (cmol/kg^-1))) + 
  labs(y=NULL) + 
  #ylab("fungal alpha diversity") +
  #ylim(1.1,3) +
  theme_bw()
p6

p7 <- ggplot(DF, aes(x = P, y = value, colour = Tree_species)) +
  geom_point(size = 1) +
  scale_colour_manual(values=myColor)+
  #geom_smooth(method = "lm", colour = "black", fill = "grey85", linewidth=0.5) +
  #scale_colour_viridis_d(option = "magma", begin = 0.2, end = 0.8) +
  xlab(bquote("Phosphorus" ~ (mg/kg^-1))) + 
  #ylab("fungal alpha diversity") +
  labs(y=NULL) + 
  #ylim(1.1,3) +
  theme_bw()
p7

p8 <- ggplot(DF, aes(x = DBH, y = value, colour = Tree_species)) +
  geom_point(size = 1) +
  scale_colour_manual(values=myColor)+
  #geom_smooth(method = "lm", colour = "black", fill = "grey85", linewidth=0.5) +
  #scale_colour_viridis_d(option = "magma", begin = 0.2, end = 0.8) +
  xlab(bquote("DBH" ~ (mm))) + 
  #ylab("fungal alpha diversity") +
  labs(y=NULL) + 
  #ylim(1.1,3) +
  theme_bw()
p8


library(ggplotify)
library(ggpubr)
ggarrange(p1, p2, p3, p4,p5,p6, p7,p8, #+ rremove("x.text"), 
          #labels = c("A", "B", "C","D"),
          align = "h",
          ncol = 3, nrow = 3,
          common.legend = TRUE, 
          legend="bottom"
) 
#(4.5*4.5)



adonis2(formula = value ~  Tree_species * Season*Layer , data = DF, 
       permutations = 9999)


#############################################
#         Beta div - Bray curtis
#######
## Import packages
library(phyloseq)
library(ggplot2)
library(ggpubr)
## Alternative to source all extra function from a github repo
source("https://raw.githubusercontent.com/mahendra-mariadassou/phyloseq-extended/master/load-extra-functions.R")
meta <- read.csv(file="../0_Metadata/Yulong_sampling_ITS2_metadata.txt","sep"="\t")
names(meta)[names(meta) == 'X'] <- 'Sample'
meta$sample <- gsub("-", ".", meta$Sample)
meta$Tree_layer <-  paste(meta$Tree_species,meta$Layer,sep="_")

phyloseq <- load("../1_data/Yulong_ITS2.rdata")
varExp <-  "Tree_species"
methods <- "bray"

## Create input and parameters dataframe
params <- data.frame( "phyloseq" = phyloseq, "varExp" = varExp, "methods" = methods)

methods <- as.list(strsplit(params$methods, ",")[[1]])

## Order samples according to grouping variable
sampleOrder <- levels(reorder(sample_names(data), as.numeric(get_variable(data, params$varExp)))) 
# 
for (method in methods){
  dist.a <- distance(data, method = method)
  a <- as.matrix(dist.a)
  write.table(a, paste(sep="", method, ".tsv"), sep="\t", quote=FALSE, col.names=NA)
  pa <- plot_dist_as_heatmap(dist.a, order = sampleOrder, title = paste("Heatmap plot of the beta distance :",method)) +
    theme(plot.title = element_text(hjust = 0.5))
  plot(pa)
}
distance <- "bray.tsv"
params <- data.frame( "phyloseq" = phyloseq, "distance" = distance, "method" = method, "varExp" = varExp)
A    <- read.table(file=params$distance, row.names=1)
dist <- as.dist(A)
ord  <- ordinate(data, method = params$method, distance = dist)
var  <- get_variable(data, params$varExp)


t <- meta %>% filter(Sample %in% row.names(A) )

adonis2(formula = dist ~ Layer * Tree_species * Season, data = t, 
       permutations = 9999)

#library(vegan)
#simper(data@otu_table, data@sam_data$Tree_species, permutations = 999)


##
meta <- meta %>% rename("Corg/N"=ratioCorg_N)

meta.z <- meta %>% filter(!is.na(pH)) %>% 
filter(!is.na(N)) %>% filter(!is.na(P)) %>% filter(!is.na(Org)) %>% filter(!is.na(CEC)) %>% filter(!is.na(Ca)) 

meta.z$pH <- (meta.z$pH - mean(meta.z$pH)/sd(meta.z$pH))
meta.z$N <- (meta.z$N - mean(meta.z$N)/sd(meta.z$N))
meta.z$P <- (meta.z$P - mean(meta.z$P)/sd(meta.z$P))
meta.z$Org <- (meta.z$Org - mean(meta.z$Org)/sd(meta.z$Org))
meta.z$CEC <- (meta.z$CEC - mean(meta.z$CEC)/sd(meta.z$CEC))
meta.z$Ca <- (log(meta.z$Ca) - mean(log(meta.z$Ca))/sd(log(meta.z$Ca)))
meta.z$Corg <- (meta.z$Corg - mean(meta.z$Corg)/sd(meta.z$Corg))
meta.z$`Corg/N`  <- (meta.z$`Corg/N` - mean(meta.z$`Corg/N`)/sd(meta.z$`Corg/N`))
meta.z$DBH  <- (meta.z$DBH - mean(meta.z$DBH)/sd(meta.z$DBH))


#meta.z <- meta

A2 <- A %>% filter(row.names(A) %in% meta.z$Sample) %>%
  select(meta.z$sample)
dist2 <- as.dist(A2)
# construct full model and calculate VIF
dbRDA.full <- capscale(dist2 ~ pH+N+P+Org+CEC+Ca+DBH,
                       meta.z)
vif.cca(dbRDA.full)
#removing selected terms from the model until all VIF scores are < 10.

# test overall significance of the analysis
anova(dbRDA.full)
# test significance of each environmental variable
anova(dbRDA.full, by = "terms")

# summary of dbRDA model to extract total variance constrained and axis scores
summary(dbRDA.full)


######
smry <- summary(dbRDA.full)
scrs <- scores(dbRDA.full)
df1  <- data.frame(smry$sites[,1:2]) # site scores for RDA1 and RDA2
df1$Sample <- rownames(df1)  #add site names
df1 <- df1 %>% left_join(meta)
df2  <- data.frame(smry$biplot[,1:2])  # mapping environmental variables
myColor <- GetMyColor(df1$Tree_species)

rda.plot <- ggplot(df1, aes(x=CAP1, y=CAP2, colour = Tree_species)) + 
  geom_point(size = 1) +
  scale_colour_manual(values=myColor)+
  geom_segment(data=df2, aes(x=0, xend=CAP1, y=0, yend=CAP2), 
               color="grey20", arrow=arrow(length=unit(0.01,"npc"))) +
  geom_text(data=df2, aes(x=CAP1,y=CAP2,label=rownames(df2),
                          hjust=0.5*(1-sign(CAP1)),vjust=0.5*(1-sign(CAP2))
                          ), 
            color="grey20", size=3.5) +

  #scale_colour_viridis_d(option = "magma", begin = 0.2, end = 0.8) +
  # geom_text(aes(label=rownames(df1),
  #               hjust=0,vjust=1.5), colour = "black",size=3) +
  geom_hline(yintercept=0, linetype="dotted") +
  geom_vline(xintercept=0, linetype="dotted") +
  #xlim(-1.05, 1.65) +
  #ylim(-1.5, 1.05) +
  xlab("RDA1 (60,9%)") + # this percentage comes from the CAP1 'importance of components:' proportion explained, which can be found in summary(dbRDA.mat) 
  ylab("RDA2 (9.0%)") + # this percentage comes from the CAP2 'importance of components:' proportion explained, which can be found in summary(dbRDA.mat) 
  #coord_fixed() +
  theme_bw() +
  stat_ellipse( type = "t", linetype = 1, level = 0.95) 
  
  
rda.plot
