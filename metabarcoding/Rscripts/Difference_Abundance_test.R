#library("readr")
library("plyr");library("dplyr");library("ggplot2");library("tidyr");library("tidyverse")
library("RColorBrewer");library("ggpubr");library("vegan")


abundance <- read.csv("../1_data/Yulong_ITS1_normalised_abundance_addCols.tsv", sep="\t" )
abundance$ID <- paste(abundance$observation_name, abundance$Species, sep="_")
abundance$ID <- str_replace_all(abundance$ID , "-", ".")

samples_metadata <- read.table("../0_Metadata/Yulong_sampling_ITS1_metadata.txt")
samples_metadata$Sample <- rownames(samples_metadata)
samples_metadata$Sample <- str_replace_all(samples_metadata$Sample , "-", ".")

ecology <-read.csv(file="../../ecology_fungaltrait_genus.txt", sep = "\t")
amplicon="ITS1"

df <-  abundance #%>% filter(observation_sum >50)
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



res_test_FC <- abundance %>% #filter(observation_sum >50) %>%
  select(c("ID", starts_with(amplicon))) %>%
  pivot_longer(cols=-c("ID"), names_to="Sample", values_to="count") %>% 
  filter(ID %in% res_test_filtered$ID) %>%
  left_join(samples_metadata ) %>%
  group_by(ID, Season) %>% 
  summarise(median=median(count),.groups = 'keep') %>% ungroup()%>%
  group_by(ID) %>% 
  summarise(max=max(median),min=min(median), highestFC=max/(min+1)) %>% 
  left_join(res_test_filtered)



df <- res_test_FC %>% filter(highestFC > 1.5 )

#group_by(PD_Order) %>% summarise(max=max(tpm_counts), mean=mean(tpm_counts), median=median(tpm_counts)) 

  #group_by(vars(all_of(variablesDetails1),  )) %>% 


#####################################


mycolorsFullList <- read.csv("../..//mycolorFile.txt",sep="\t")


m <- abundance %>% filter(ID %in% df$ID) %>%
  select(c("ID", starts_with(amplicon)))

rownames(m) <- m$ID
m$ID <- NULL

library(pheatmap)
pheatmap(m)

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



annOTU <- abundance %>% select( ID,Genus,Order ) %>%
                        filter(ID %in% df$ID) %>%
                        left_join(ecology) %>%
                        mutate(primary_lifestyle=case_when(!(is.na(primary_lifestyle)) ~ primary_lifestyle,
                                                           TRUE ~ "Unidentified" )) %>%
                        select(ID, primary_lifestyle, Order)

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





annotation_colors <- list("Tree_species"=annotation_colors_tree, 
                          "Cluster"=annotation_colors_Cluster,
                          "Season"=annotation_colors_season,
                          "Layer"=annotation_colors_Layer,
                          "primary_lifestyle"=annotation_colors_ecology,
                          "Order"=annotation_colors_order)

#df$Phylum <- factor(df$Phylum, levels=mycolors$Name)


pheatmap(log(m+1),colorRampPalette(c("white",brewer.pal(9, "Blues"),"black"))(100),
         #colorRampPalette(c("white", "#e6f5ff", "navyblue"), space = "Lab")(100)#,
          annotation_colors = annotation_colors,
            annotation_col = annSample,
         annotation_row = annOTU,
         #cutree_rows = 3,
         clustering_method = "ward.D2",
         show_colnames = FALSE,
         #fontsize_col = 1,
         fontsize = 7 ,
         angle_col = 45
)





##############################################################


m2 <- abundance %>% #filter(observation_sum >50) %>%
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

annSample <- samples_metadata %>% select(Cluster, Tree_species) %>% unique()
rownames(annSample) <- annSample$Cluster
annSample$Cluster <- NULL

annotation_colors <- list("Tree_species"=annotation_colors_tree,
                          "primary_lifestyle"=annotation_colors_ecology,
                          "Order"=annotation_colors_order)



pheatmap(log(m2+1),colorRampPalette(c("white",brewer.pal(9, "Blues"),"black"))(100),
         #colorRampPalette(c("white", "#e6f5ff", "navyblue"), space = "Lab")(100)#,
         annotation_colors = annotation_colors,
         annotation_col = annSample,
         annotation_row = annOTU,
         #cutree_rows = 3,
         clustering_method = "ward.D2",
         #show_colnames = FALSE,
         #fontsize_col = 1,
         fontsize = 7 ,
         angle_col = 45
)


#####################



m2 <- abundance %>% #filter(observation_sum >50) %>%
  select(c("ID", starts_with(amplicon))) %>%
  pivot_longer(cols=-c("ID"), names_to="Sample", values_to="count") %>% 
  filter(ID %in% df$ID) %>%
  left_join(samples_metadata ) %>%
  #mutate(Season_cluster=paste(Cluster,Layer,sep="_")) %>%
  group_by(ID, Layer) %>% 
  summarise("counts"=sum(count)) %>%
  pivot_wider(names_from = Layer, values_from =counts )

tmp <- m2
m2$ID <- NULL
rownames(m2) <- tmp$ID
rm(tmp)

annSample <- samples_metadata %>% 
              #mutate(Season_cluster=paste(Cluster,Layer,sep="_")) %>%
              select(Season_cluster, Layer,Cluster, Tree_species) %>% unique()

rownames(annSample) <- annSample$Season_cluster
annSample$Season_cluster <- NULL

annotation_colors <- list("Tree_species"=annotation_colors_tree,
                          "Layer"=annotation_colors_Layer,
                          "Cluster"=annotation_colors_Cluster,
                          "primary_lifestyle"=annotation_colors_ecology,
                          "Order"=annotation_colors_order)



pheatmap((m2+1),colorRampPalette(c("white",brewer.pal(9, "Blues"),"black"))(100),
         #colorRampPalette(c("white", "#e6f5ff", "navyblue"), space = "Lab")(100)#,
         annotation_colors = annotation_colors,
         #annotation_col = annSample,
         annotation_row = annOTU,
         #cutree_rows = 3,
         clustering_method = "ward.D2",
         #show_colnames = FALSE,
         #fontsize_col = 1,
         fontsize = 7 ,
         angle_col = 45
)



# 
# 
# 
# 
# 
# 
# 
# m <- df2
# m$ID <- NULL
# rownames(m) <- df2$ID
# 
# 
# df3 <- merge(df,filteredFungalDB2, all.x = T)
# df3$guild_fg[is.na(df3$guild_fg)] <- "Unidentified"
# 
# ann <- df3 %>% select(guild_fg, Order,Phylum )
# rownames(ann) <- df3$ID
# 
# topNelement <- df3 %>%  group_by(Order) %>% mutate(nbObs =n() ) %>%
#   select(Order, nbObs) %>% unique() %>% ungroup() %>%
#   filter(nbObs > 4)
#   #slice_max(nbObs, n = 7)
# 
# ann$Order[!(ann$Order %in% topNelement$Order) & !(ann$Order == "unidentified")] <- "Other"
# 
# #annSample <- samples_metadata %>% select(Year_Season, plot_name, layer)
# annSample <- samples_metadata %>% select(Year_Season, plot_name, cluster_number) %>% unique()
# rownames(annSample) <-  paste(annSample$cluster_number, annSample$Year_Season, sep="_" )
# annSample$cluster_number <- NULL
# 
# pheatmap(log(m+1),
#          colorRampPalette(c("#e6f5ff", "navyblue"), space = "Lab")(100),
#          clustering_method = "ward.D2",
#          annotation_row = ann,
#          annotation_col = annSample,
#          annotation_colors = annotation_colors,
#          #show_colnames = FALSE,
#          fontsize = 8,
# angle_col = 45,
# cutree_rows = 3,
# #cutree_cols =7
# #fontsize_row = 3
# )
# 
# 
# ?pheatmap
# 
# 
# 
# pheatmap(log(m+1))
# 
# 
# 
# 
# 
# ###  Kruskal ---------------------------------------------------
# 
# Fungi_ratio <- Fungi_TPM %>%
#   dplyr::select(-starts_with("NR_"),-starts_with("Myc"),-ID,-ctg_length, -sum_contig_count) %>% 
#   filter(!is.na(CAZy_family)) %>%
#   #filter(grepl('saprotroph',LifeStyle) ) %>%
#   #mutate(CAZy_family = case_when(CAZy_family %in% top10CAZymes$CAZy_family ~ CAZy_family, TRUE ~ "Others")) %>% 
#   group_by(CAZy_family) %>% summarise_if(is.numeric, sum) 
# 
# df <- t(Fungi_ratio)
# colnames(df) <- df[1,]
# df <- df[-1,]
# df2 <- data.frame(df)
# df2$samples <- rownames(df)
# df2 <- df2 %>% left_join(treatment_metadata)
# df2$layer_tree <- paste( df2$layer,df2$tree_species, sep="_")
# # %>% pivot_longer(cols=-c(CAZy_family), names_to="samples", values_to="tpm_counts") %>%
# #  left_join(treatment_metadata)
# 
# #colnames(df)
# measures <- unique(colnames(df))
# measures<- str_replace_all(measures , "-", ".")
# variables <- c("layer", "tree_species","layer_tree")
# 
# DF <- data.frame()
# DF2 <- data.frame(matrix(ncol = length(variables) , nrow = length(measures) ))
# colnames(DF2) <- variables
# rownames(DF2) <- measures
# 
# for (v in variables){
#   for (m in measures){
#     f <- paste(m," ~ ", v)
#     kruskal_res <- kruskal.test(as.formula(f), data=df2 )
#     #kruskal_res<- kruskal.test(AA3 ~ tree_species, data=df2 )
#     DF[m,v]<- kruskal_res$p.value
#   }
#   DF2[,v] <- p.adjust(DF[,v], method = "fdr")
# }
# 
# rm(DF, m,v,f,kruskal_res)
# DF2$CAZyme_family <- rownames(DF2)
# DF <- DF2 %>% filter_at(vars(variables), any_vars(. < 0.01 ))
# #DF <- DF %>% filter(tree_species < 0.01 )
