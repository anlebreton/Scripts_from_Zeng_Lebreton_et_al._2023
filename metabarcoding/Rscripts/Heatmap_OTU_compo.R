#library("readr")
library("plyr");library("dplyr");library("ggplot2");library("tidyr");library("tidyverse")
library("RColorBrewer");library("ggpubr");library("vegan")
library(pheatmap)

abundance <- read.csv("../1_data/Yulong_ITS1_normalised_abundance_addCols.tsv", sep="\t" )
#abundance$ID <- paste(abundance$observation_name, abundance$Species, sep="_")
abundance$ID <- paste(abundance$observation_name, abundance$Order, sep="_")

abundance$ID <- str_replace_all(abundance$ID , "-", ".")
abundance$ID <- str_replace_all(abundance$ID , " ", ".")
abundance$ID <- str_replace_all(abundance$ID , "[(]", ".")
abundance$ID <- str_replace_all(abundance$ID , "[)]", ".")

samples_metadata <- read.table("../0_Metadata/Yulong_sampling_ITS1_metadata.txt")
samples_metadata$Sample <- rownames(samples_metadata)
samples_metadata$Sample <- str_replace_all(samples_metadata$Sample , "-", ".")

ecology <-read.csv(file="../../ecology_fungaltrait_genus.txt", sep = "\t")
amplicon="ITS1"
eco <- ecology %>% select(Genus, primary_lifestyle)
abundance <- left_join(abundance,eco)

# 
# df <-  abundance #%>% filter(observation_sum >50)
# #df <-  abundance %>% filter(observation_name %in% DE$observation_name)
# #df$ID <- paste(df$observation_name, df$Species, sep="_")
# #df$ID <- str_replace_all(df$ID , "-", ".")
# 
# 
# df2 <- df %>% select(c("ID", starts_with(amplicon))) #%>%
# #pivot_longer(cols=-c("ID"), names_to="samples", values_to="count") %>%
# #group_by("ID","samples") %>% summarise("counts"=sum(count)) 
# 
# tmp <- df2
# df2$ID <- NULL
# rownames(df2) <- tmp$ID
# 
# df2 <- as.data.frame(t(df2))
# df2$Sample <- rownames(df2)
# df2 <- df2 %>% left_join(samples_metadata )
# #df2$layer_tree <- paste( df2$layer,df2$tree_species, sep="_")
# 
##################
abundance$SpeciesID <- paste(abundance$observation_name, abundance$Species, sep="_")
abundance$observation_perc <- (abundance$observation_sum/sum(abundance$observation_sum))*100
number_samples <- length(samples_metadata$Sample)


# res <-  abundance %>%  #filter(percSample > 50) #%>%
#   arrange(desc(observation_perc)) %>%
#   mutate(cumul_observation_perc=cumsum(observation_perc)) %>%
#   filter(cumul_observation_perc <= 90)
# lenValue <- length(res$SpeciesID)
# 
# res <-  abundance %>%  #filter(percSample > 50) #%>%
#   arrange(desc(observation_perc)) %>%
#   mutate(cumul_observation_perc=cumsum(observation_perc)) %>%
#   slice_head(n=lenValue+1)
# 
# 
# abundance <- res

################

# res_test_FC <- abundance %>% #filter(observation_sum >50) %>%
#   select(c("ID", starts_with(amplicon))) %>%
#   pivot_longer(cols=-c("ID"), names_to="Sample", values_to="count") %>% 
#   #filter(ID %in% res_test_filtered$ID) %>%
#   left_join(samples_metadata ) %>%
#   group_by(ID, Tree_species) %>% 
#   summarise(median=median(count),.groups = 'keep') %>% ungroup()%>%
#   group_by(ID) %>% 
#   summarise(max=max(median),min=min(median), highestFC=max/(min+1)) #%>% 
#   #left_join(res_test_filtered)
# 
# df <- res_test_FC
df <- abundance
#######################################@
mycolorsFullList <- read.csv("../..//mycolorFile.txt",sep="\t")

m <- abundance %>% #filter(ID %in% df$ID) %>%
  select(c("ID", starts_with(amplicon)))



rownames(m) <- m$ID
m$ID <- NULL
#######

############


annSample <- samples_metadata %>% select(Sample, Layer, Season, Cluster, Tree_species) %>% unique()
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


topOTUorder <- abundance %>% 
  filter(primary_lifestyle %in% c("dung_saprotroph","litter_saprotroph","soil_saprotroph",
                                                               "unspecified_saprotroph","wood_saprotroph")) %>%
  
  select( ID,Order ) %>%
  mutate(Order=case_when((Order !="unknown order") ~ Order,
                         TRUE ~ "Unidentified" )) %>%
  mutate(Order=case_when((Order !="Multi-affiliation" ) ~ Order,
                         TRUE ~ "Unidentified" )) %>%
  filter(ID %in% df$ID) %>% group_by(Order) %>%
  summarise(nb=n()) %>% 
  slice_max(nb, n=10)

topOTUphylum <- abundance %>%
  filter(primary_lifestyle %in% c("dung_saprotroph","litter_saprotroph","soil_saprotroph",
                                  "unspecified_saprotroph","wood_saprotroph")) %>%
  select( ID,Phylum ) %>%
  filter(ID %in% df$ID) %>% group_by(Phylum) %>%
  summarise(nb=n()) %>% 
  slice_max(nb, n=3)

topOTUorder <- abundance %>% 
  filter(primary_lifestyle %in% c("dung_saprotroph","litter_saprotroph","soil_saprotroph",
                                  "unspecified_saprotroph","wood_saprotroph")) %>%
  select( ID,Order ) %>%
  filter(ID %in% df$ID) %>% group_by(Order) %>%
  summarise(nb=n()) %>% 
  slice_max(nb, n=5)

topOTUclass <- abundance %>% select( ID,Class ) %>%
  filter(ID %in% df$ID) %>% group_by(Class) %>%
  summarise(nb=n()) %>% 
  slice_max(nb, n=3)

annOTU <- abundance %>% #select( ID,Genus,Order ) %>%
  select(ID,Phylum, Class, Order, primary_lifestyle) %>%
  filter(ID %in% df$ID) %>%
  mutate(Phylum=case_when((Phylum %in% topOTUphylum$Phylum) ~ Phylum,
                          TRUE ~ "Others" )) %>%
  mutate(Class=case_when((Class %in% topOTUclass$Class) ~ Class, TRUE ~ "Others" ))%>% 
  #left_join(ecology) %>%
  mutate(primary_lifestyle=case_when(!(is.na(primary_lifestyle)) ~ primary_lifestyle,
                                     TRUE ~ "Unidentified" )) %>%
   mutate(Order=case_when((Order !="unknown order") ~ Order,
    TRUE ~ "Unidentified" )) %>%
   mutate(Order=case_when((Order !="Multi-affiliation" ) ~ Order,
    TRUE ~ "Unidentified" )) %>%
  mutate(Order=case_when(Order %in% topOTUorder$Order ~ Order, TRUE ~ "Others" ))%>% 
    #select(ID, primary_lifestyle,Order) %>%
  select(ID, Order,Class,primary_lifestyle) 
  #select(ID,Phylum)


rownames(annOTU) <- annOTU$ID  
annOTU$ID <- NULL

#mycolors <-  mycolorsFullList %>% filter(Name %in% annOTU$primary_lifestyle)
mycolors <-  mycolorsFullList %>% filter(Name %in% c("dung_saprotroph","litter_saprotroph","soil_saprotroph",
  "unspecified_saprotroph","wood_saprotroph"))
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

mycolors <-  mycolorsFullList %>% filter(Name %in% annOTU$Class)
mycolors <- mycolors %>% arrange(myLevel)
annotation_colors_Class <- mycolors$color
names(annotation_colors_Class) <- mycolors$Name


#######

m2 <- abundance %>% #filter(observation_sum >50) %>%
  #filter(primary_lifestyle == "ectomycorrhizal") %>%
  filter(primary_lifestyle %in% c("dung_saprotroph","litter_saprotroph","soil_saprotroph",
                                  "unspecified_saprotroph","wood_saprotroph")) %>%
  select(c("ID", starts_with(amplicon))) %>%
  pivot_longer(cols=-c("ID"), names_to="Sample", values_to="count") %>% 
  filter(ID %in% df$ID) %>%
  left_join(samples_metadata ) %>%
  group_by(ID, Cluster) %>% 
  summarise("counts"=sum(count)) %>%
  pivot_wider(names_from = Cluster, values_from =counts )

tmp <- m2
m2$ID <- NULL
rownames(m2) <- tmp$ID
rm(tmp)

#annSample <- samples_metadata %>% select(Cluster, Tree_species) %>% unique()
annSample <- samples_metadata %>% select(Cluster) %>% unique()
rownames(annSample) <- annSample$Cluster
#annSample$Cluster <- NULL

annotation_colors <- list(#"Tree_species"=annotation_colors_tree,
                          "Cluster"=annotation_colors_Cluster,
                          "primary_lifestyle"=annotation_colors_ecology,
                          "Order"=annotation_colors_order,
                          "Class"=annotation_colors_Class)


#,"navyblue","black","black"
pheatmap(log(m2+1),colorRampPalette(c(brewer.pal(9, "Blues"),"navyblue","navyblue"))(100),
         #colorRampPalette(c("white", "#e6f5ff", "navyblue"), space = "Lab")(100)#,
         annotation_colors = annotation_colors,
         annotation_col = annSample,
         annotation_row = annOTU,
         #cutree_rows = 3,
         #clustering_method = "ward.D2",
         clustering_method = "complete",
         #show_rownames = FALSE,
         fontsize_col = 7,
         fontsize_row = 1 ,
         fontsize=7,
         angle_col = 45
         )



