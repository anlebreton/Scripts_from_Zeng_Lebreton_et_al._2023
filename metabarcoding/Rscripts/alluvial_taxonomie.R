#library("readr")
library("plyr");library("dplyr");library("ggplot2");library("tidyr");library("tidyverse")
library("RColorBrewer");library("ggpubr");library("vegan")
library(ggalluvial)


GetMyColor2 <- function(colDF){
  mycolorsFile <- read.csv("../0_Metadata/mycolorFile_v2.txt",sep="\t")
  mycolors <- mycolorsFile %>% filter(Name %in% colDF)
  mycolors <- mycolors %>% arrange(myLevel)
  mycolor <- mycolors$color
  names(mycolor) <- mycolors$Name
  #colDF <- factor(colDF, levels=mycolors$Name
  return(mycolor)
}
GetMyColor <- function(colDF){
  mycolorsFile <- read.csv("../../mycolorFile.txt",sep="\t")
  mycolors <- mycolorsFile %>% filter(Name %in% colDF)
  mycolors <- mycolors %>% arrange(myLevel)
  mycolor <- mycolors$color
  names(mycolor) <- mycolors$Name
  #colDF <- factor(colDF, levels=mycolors$Name
  return(mycolor)
}

# FROGS BIOM to TSV results with the taxonomy (one column per level) added
abundance <- read.csv("../1_data/Yulong_ITS2_normalised_abundance_addCols.tsv", sep="\t" )
abundance$ID <- paste(abundance$observation_name, abundance$Species, sep="_")

#fix some names
abundance$ID <- str_replace_all(abundance$ID , "-", ".")
abundance$ID <- str_replace_all(abundance$ID , " ", ".")
abundance$ID <- str_replace_all(abundance$ID , "[(]", ".")
abundance$ID <- str_replace_all(abundance$ID , "[)]", ".")

#fix some names
samples_metadata <- read.table("../0_Metadata/Yulong_sampling_ITS2_metadata.txt")
samples_metadata$Sample <- rownames(samples_metadata)
samples_metadata$Sample <- str_replace_all(samples_metadata$Sample , "-", ".")

#obtained from the supplementary of the Fungal trait article
ecology <-read.csv(file="../../0_Metadata/ecology_fungaltrait_genus.txt", sep = "\t")

# selecting purpose, here all the sample names start by the amplicon tag
amplicon="ITS2"

# simplify the full table
df2 <- abundance %>% select(c("ID", starts_with(amplicon))) 

tmp <- df2
df2$ID <- NULL
rownames(df2) <- tmp$ID

df2 <- as.data.frame(t(df2))
df2$Sample <- rownames(df2)
df2 <- df2 %>% left_join(samples_metadata )

# here measures = OTU names
measures <- df2 %>% select(starts_with("Cluster")) 
measures <- unique(colnames(measures))

measures<- str_replace_all(measures , "-", ".")
variables <- c("Layer", "Tree_species","Season")
v1 <- c("Tree_species","Season")
v2 <- c("Layer")


# perform a Wilcoxon test to identify differences among Layer (2 categories) 
# perform a Kruskal test to identify differences among Seasons and Tree associations
# FDR corrections for both test
DF <- data.frame()
res_test <- data.frame(matrix(ncol = length(variables) , nrow = length(measures) ))
colnames(res_test) <- variables
rownames(res_test) <- measures

for (v in v1){
  for (m in measures){
    f <- paste(m," ~ ", v)
    kruskal_res <- kruskal.test(as.formula(f), data=df2 )
    DF[m,v]<- kruskal_res$p.value
  }
  res_test[,v] <- p.adjust(DF[,v], method = "fdr")
}


for (v in v2){
  for (m in measures){
    f <- paste(m," ~ ", v)
    if(f=="Cluster  ~  Layer") next
    #print(f)
    wilcox.test <- wilcox.test(as.formula(f), data=df2 )
    DF[m,v]<- wilcox.test$p.value
  }
  res_test[,v] <- p.adjust(DF[,v], method = "fdr")
}

#  cleanup
rm(DF, m,v,v1,v2,f,kruskal_res,wilcox.test,df,tmp)

# filter to remove OTU with no significant differences
res_test$ID <- rownames(res_test)
res_test_filtered <- res_test %>% filter(ID != "Cluster") %>%
  filter_at(vars(all_of(variables)), any_vars(. < 0.00001 ))

#######################################@

variablesDetails1 <- c("Quercus","Abies","Picea")
variablesDetails2 <- c(unique(samples_metadata$Season))

# compute FC
res_test_FC <- abundance %>% 
  select(c("ID", starts_with(amplicon))) %>%
  pivot_longer(cols=-c("ID"), names_to="Sample", values_to="count") %>% 
  filter(ID %in% res_test_filtered$ID) %>%
  left_join(samples_metadata ) %>% ungroup()%>%
  group_by(ID, Tree_species) %>% 
  summarise(median=median(count),.groups = 'keep') %>% ungroup()%>%
  group_by(ID) %>% 
  summarise(max=max(median),min=min(median), highestFC=max/(min+1)) %>% 
  left_join(res_test_filtered)


# select top 20 
res <- res_test_FC %>% filter(Tree_species < 10e-5) %>% arrange(desc(max)) %>%
  slice_head(n=20)


res2 <- abundance %>% filter(ID %in% res$ID) %>% 
  select(c("ID", starts_with(amplicon),"Phylum", "Order")) %>%
  pivot_longer(cols=-c("ID","Phylum","Order"), names_to="Sample", values_to="count") %>% 
  left_join(samples_metadata ) %>%
  group_by(ID, Tree_species, Phylum,Order ) %>% 
  summarise(median=median(count))

tmp <- abundance %>% filter(ID %in% res$ID) %>% select("ID","Phylum","Order")
#write.table(tmp, file="tmp.tsv", sep="\t", row.names = F)

# check format
is_alluvia_form(as.data.frame(res2), axes = 1:2, silent = TRUE)

#get colors
myColorTree = GetMyColor(res2$Tree_species)
myColorOrder = GetMyColor2(res2$Order)
myColorPhylum =  GetMyColor2(res2$Phylum)

# sort factor
res2$Phylum <- factor(res2$Phylum, levels =names(myColorPhylum ))
res2$Order <- factor(res2$Order, levels =names(myColorOrder) )

# alluvial plot
ggplot(as.data.frame(res2),
       aes(y = median, axis1=ID, axis2 = Tree_species)) +
  geom_alluvium(aes(fill = Tree_species), width = 1/12) +
  geom_stratum(width = 1/12, color = "grey") +
  scale_fill_manual(values = myColorTree) +
  geom_text(
    aes(
      label = after_stat(stratum),
      hjust=1,
      color = "grey"
    ),
    stat = "stratum", fontface = "bold", 
    size = 3
  )+
  scale_x_discrete(limits = c("OTU", "Tree species"), expand = c(0.5, 0.5)) +
  scale_y_continuous(breaks = NULL) +
  theme_minimal() +
  theme(
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_blank()
  ) +
  theme(legend.position = "none") +
  ylab(NULL)

