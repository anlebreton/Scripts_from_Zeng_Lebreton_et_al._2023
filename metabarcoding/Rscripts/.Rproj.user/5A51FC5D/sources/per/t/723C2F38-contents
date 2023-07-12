#library("readr")
library("plyr");library("dplyr");library("ggplot2");library("tidyr");library("tidyverse")
library("RColorBrewer");library("ggpubr");library("vegan")


abundance <- read.csv("../1_data/Yulong_ITS1_normalised_abundance_addCols.tsv", sep="\t" )
abundance$ID <- paste(abundance$observation_name, abundance$Species, sep="_")
#abundance$ID <- paste(abundance$observation_name, abundance$Order, sep="_")

abundance$ID <- str_replace_all(abundance$ID , "-", ".")
abundance$ID <- str_replace_all(abundance$ID , " ", ".")
abundance$ID <- str_replace_all(abundance$ID , "[(]", ".")
abundance$ID <- str_replace_all(abundance$ID , "[)]", ".")


samples_metadata <- read.table("../0_Metadata/Yulong_sampling_ITS1_metadata.txt")
samples_metadata$Sample <- rownames(samples_metadata)
samples_metadata$Sample <- str_replace_all(samples_metadata$Sample , "-", ".")

ecology <-read.csv(file="../../ecology_fungaltrait_genus.txt", sep = "\t")
amplicon="ITS1"

df <-  abundance %>% left_join(ecology) %>% filter(grepl(pattern = "sapro", primary_lifestyle)) 

#%>% filter(observation_sum >50)
#df <-  abundance %>% filter(observation_name %in% DE$observation_name)
#df$ID <- paste(df$observation_name, df$Species, sep="_")
#df$ID <- str_replace_all(df$ID , "-", ".")


df2 <- df %>% select(c("ID", starts_with(amplicon))) #%>%
#pivot_longer(cols=-c("ID"), names_to="samples", values_to="count") %>%
#group_by("ID","samples") %>% summarise("counts"=sum(count)) 

tmp <- df2
df2$ID <- NULL
rownames(df2) <- tmp$ID

df2 <- as.data.frame(t(df2))
df2$Sample <- rownames(df2)
df2 <- df2 %>% left_join(samples_metadata )
#df2$layer_tree <- paste( df2$layer,df2$tree_species, sep="_")

#colnames(df)
measures <- df2 %>% select(starts_with("Cluster")) 
measures <- unique(colnames(measures))

measures<- str_replace_all(measures , "-", ".")
variables <- c("Layer", "Tree_species","Season")
v1 <- c("Tree_species","Season")
v2 <- c("Layer")

DF <- data.frame()
res_test <- data.frame(matrix(ncol = length(variables) , nrow = length(measures) ))
colnames(res_test) <- variables
rownames(res_test) <- measures

for (v in v1){
  for (m in measures){
    f <- paste(m," ~ ", v)
    kruskal_res <- kruskal.test(as.formula(f), data=df2 )
    #kruskal_res<- kruskal.test(AA3 ~ tree_species, data=df2 )
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
    #kruskal_res<- kruskal.test(AA3 ~ tree_species, data=df2 )
    DF[m,v]<- wilcox.test$p.value
  }
  res_test[,v] <- p.adjust(DF[,v], method = "fdr")
}

rm(DF, m,v,v1,v2,f,kruskal_res,wilcox.test,df,tmp)
res_test$ID <- rownames(res_test)
res_test_filtered <- res_test %>% filter(ID != "Cluster") %>%
  filter_at(vars(all_of(variables)), any_vars(. < 0.00001 ))

#######################################@

variablesDetails1 <- c("Quercus","Abies","Picea")
variablesDetails2 <- c(unique(samples_metadata$Season))

res_test_FC <- abundance %>% #filter(observation_sum >50) %>%
  select(c("ID", starts_with(amplicon))) %>%
  pivot_longer(cols=-c("ID"), names_to="Sample", values_to="count") %>% 
  filter(ID %in% res_test_filtered$ID) %>%
  left_join(samples_metadata ) %>%
  group_by(ID, Tree_species) %>% 
  summarise(median=median(count),.groups = 'keep') %>% ungroup()%>%
  group_by(ID) %>% 
  summarise(max=max(median),min=min(median), highestFC=max/(min+1)) %>% 
  left_join(res_test_filtered)



# res_test_FC <- abundance %>% #filter(observation_sum >50) %>%
#   select(c("ID", starts_with(amplicon))) %>%
#   pivot_longer(cols=-c("ID"), names_to="Sample", values_to="count") %>%
#   filter(ID %in% res_test_filtered$ID) %>%
#   left_join(samples_metadata ) %>%
#   group_by(ID, Layer) %>%
#   summarise(median=median(count),.groups = 'keep') %>% ungroup()%>%
#   group_by(ID) %>%
#   summarise(max=max(median),min=min(median), highestFC=max/(min+1)) %>%
#   left_join(res_test_filtered)



df <- res_test_FC %>% filter(highestFC > 1.5)

#group_by(PD_Order) %>% summarise(max=max(tpm_counts), mean=mean(tpm_counts), median=median(tpm_counts)) 

#group_by(vars(all_of(variablesDetails1),  )) %>% 


#####################################


mycolorsFullList <- read.csv("../..//mycolorFile.txt",sep="\t")


m <- abundance %>% filter(ID %in% df$ID) %>%
  select(c("ID", starts_with(amplicon)))

rownames(m) <- m$ID
m$ID <- NULL

library(pheatmap)
#pheatmap(m)

annSample <- samples_metadata %>% select(Sample, Layer, Season, #Cluster,
                                         Tree_species) %>% unique()
# rownames(annSample) <-  paste(annSample$Cluster, annSample$Season, sep="_" )
rownames(annSample) <- annSample$Sample
annSample$Sample <- NULL
#samples_metadata$Sample <- str_replace_all(samples_metadata$Sample , "-", ".")


mycolors <-  mycolorsFullList %>% filter(Name %in% annSample$Season)
mycolors <- mycolors %>% arrange(myLevel)
annotation_colors_season <- mycolors$color
names(annotation_colors_season) <- mycolors$Name

mycolors <-  mycolorsFullList %>% filter(Name %in% annSample$Tree_species)
mycolors <- mycolors %>% arrange(myLevel)
annotation_colors_tree <- mycolors$color
names(annotation_colors_tree) <- mycolors$Name

mycolors <-  mycolorsFullList %>% filter(Name %in% annSample$Cluster)
mycolors <- mycolors %>% arrange(myLevel)
annotation_colors_Cluster <- mycolors$color
names(annotation_colors_Cluster) <- mycolors$Name

mycolors <-  mycolorsFullList %>% filter(Name %in% annSample$Layer)
mycolors <- mycolors %>% arrange(myLevel)
annotation_colors_Layer <- mycolors$color
names(annotation_colors_Layer) <- mycolors$Name


topOTUorder <- abundance %>% select( ID,Order ) %>%
  mutate(Order=case_when((Order !="unknown order") ~ Order,
                         TRUE ~ "Unidentified" )) %>%
  mutate(Order=case_when((Order !="Multi-affiliation" ) ~ Order,
                         TRUE ~ "Unidentified" )) %>%
  filter(ID %in% df$ID) %>% group_by(Order) %>%
  summarise(nb=n()) %>% 
  slice_max(nb, n=3)

topOTUphylum <- abundance %>% select( ID,Phylum ) %>%
  filter(ID %in% df$ID) %>% group_by(Phylum) %>%
  summarise(nb=n()) %>% 
  slice_max(nb, n=8)

annOTU <- abundance %>% #select( ID,Genus,Order ) %>%
  select(ID,Phylum, Genus) %>%
  filter(ID %in% df$ID) %>%
  # mutate(Phylum=case_when((Phylum %in% topOTUphylum$Phylum) ~ Phylum,
  #                         TRUE ~ "Others" )) %>%
  left_join(ecology) %>% 
  #mutate(primary_lifestyle=case_when(!(is.na(primary_lifestyle)) ~ primary_lifestyle,
  #                                   TRUE ~ "Unidentified" )) %>%
   mutate(Order=case_when((Order !="unknown order") ~ Order,
    TRUE ~ "Unidentified" )) %>%
   mutate(Order=case_when((Order !="Multi-affiliation" ) ~ Order,
    TRUE ~ "Unidentified" )) %>%
  mutate(Order=case_when(Order %in% topOTUorder$Order ~ Order, TRUE ~ "Others" ))%>% 
    select(ID, Order,primary_lifestyle) 
  #select(ID,Phylum,primary_lifestyle)



rownames(annOTU) <- annOTU$ID  
annOTU$ID <- NULL

mycolors <-  mycolorsFullList %>% filter(Name %in% annOTU$primary_lifestyle)
mycolors <- mycolors %>% arrange(myLevel)
annotation_colors_ecology <- mycolors$color
names(annotation_colors_ecology) <- mycolors$Name

mycolors <-  mycolorsFullList %>% filter(Name %in% annOTU$Order)
mycolors <- mycolors %>% arrange(myLevel)
annotation_colors_order <- mycolors$color
names(annotation_colors_order) <- mycolors$Name

mycolors <-  mycolorsFullList %>% filter(Name %in% annOTU$Phylum)
mycolors <- mycolors %>% arrange(myLevel)
annotation_colors_Phylum <- mycolors$color
names(annotation_colors_Phylum) <- mycolors$Name



annotation_colors <- list("Tree_species"=annotation_colors_tree, 
                          "Cluster"=annotation_colors_Cluster,
                          "Season"=annotation_colors_season,
                          "Layer"=annotation_colors_Layer,
                          "Phylum"=annotation_colors_Phylum,
                          "Order"=annotation_colors_order,
                          "primary_lifestyle"=annotation_colors_ecology
)

#df$Phylum <- factor(df$Phylum, levels=mycolors$Name)





pheatmap(log(m+1),colorRampPalette(c("white",brewer.pal(9, "Blues"),"black"))(100),
         #colorRampPalette(c("white", "#e6f5ff", "navyblue"), space = "Lab")(100)#,
         annotation_colors = annotation_colors,
         annotation_col = annSample,
         annotation_row = annOTU,
         #cutree_rows = 2,
         cutree_cols = 2,
         clustering_method = "ward.D2",
         show_colnames = FALSE,
         #fontsize_col = 1,
         fontsize = 7 ,
         #fontsize_row = 3,
         angle_col = 45
)

