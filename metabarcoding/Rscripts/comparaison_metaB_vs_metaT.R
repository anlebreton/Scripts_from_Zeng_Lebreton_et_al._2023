library("plyr");library("dplyr");library("ggplot2");library("tidyr");library("tidyverse")
library("RColorBrewer");library("ggpubr");library("vegan")
#library(pheatmap)
#library(RColorBrewer)

mycolorsFullList <- read.csv("../..//mycolorFile.txt",sep="\t")

samples_metadata <- read.table("../../metabarcoding/0_metadata/Yulong_sampling_ITS2_metadata.txt")
samples_metadata$samples <- rownames(samples_metadata)
samples_metadata$samples <- str_replace_all(samples_metadata$samples , "-", ".") 
samples_metadata <- samples_metadata %>% mutate(sample=paste(Layer,Tree_number, sep=".") )
  

data <- read.csv("../../metabarcoding/1_data/Yulong_ITS2_normalised_abundance_addCols.tsv", sep='\t')
data$ID <- paste(data$observation_name, data$Species, sep="_")
ecology <-read.csv(file="../../ecology_fungaltrait_genus.txt", sep = "\t")

amplicon="ITS2"

OTUinfo <- data %>% select(ID, Genus) %>% unique() %>%
                left_join(ecology) %>%
                select(ID,Genus,primary_lifestyle)

############################### get only wet season 2019 ########
samples_metadata_subset <- samples_metadata %>% select(Tree_species, Cluster, Layer, sample) %>% unique()
rownames(samples_metadata_subset) <- samples_metadata_subset$sample

abundance <- data %>%
  #dplyr::select(-starts_with("NR_"),-starts_with("Myc"),-ID,-ctg_length, -sum_contig_count) %>% 
  dplyr::select(starts_with(amplicon), ID, Genus) %>%
  #group_by(Phylum) %>% summarise_if(is.numeric, sum) %>%
  pivot_longer(cols=-c(ID,Genus), names_to="samples", values_to="counts")%>%
  left_join(., samples_metadata, by="samples") %>% 
  filter(Season == "2019_wet") %>% 
  mutate(sample=paste(Layer,Tree_number, sep=".") )%>% 
  group_by(ID, sample) %>% 
  summarise("reads_count"=sum(counts), .groups="keep") %>% 
  left_join(OTUinfo) %>%
  left_join(samples_metadata_subset) %>%
  mutate(primary_lifestyle=case_when(!(is.na(primary_lifestyle)) ~ primary_lifestyle, TRUE ~ "Unidentified"))
  

toplifestyle <- abundance %>% select( ID,primary_lifestyle, reads_count ) %>%
  #filter(ID %in% df$ID) %>% 
  group_by(primary_lifestyle) %>%
  summarise(val=sum(reads_count)) %>% 
  slice_max(val, n=8)

df <- abundance %>% group_by(Layer,Tree_species, sample,primary_lifestyle )%>% 
  mutate(primary_lifestyle=case_when(primary_lifestyle %in% toplifestyle$primary_lifestyle ~ primary_lifestyle, TRUE ~ "Others"))%>%
  summarise("count"=sum(reads_count),.groups="keep")

mycolors <-  mycolorsFullList %>% filter(Name %in% df$primary_lifestyle)
mycolors <- mycolors %>% arrange(myLevel)
myColorEco <- mycolors$color
names(myColorEco) <- mycolors$Name
df$primary_lifestyle<- factor(as.character(df$primary_lifestyle), levels=mycolors$Name)
# 

ggplot(data=df, aes(x=sample, y=count, fill=primary_lifestyle)) +
  geom_col() + facet_wrap(~ Layer*Tree_species, scales = "free_x") +
  scale_fill_manual(values=myColorEco,limits =names(myColorEco) )+
  #c("#999999",rev(c("#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7","springgreen4","mediumorchid4","cadetblue","lightpink","khaki","turquoise1")))) +
  theme(axis.text.x=element_text(angle = 55,hjust=1))
#c(grey="#999999", orange="#E69F00", lightblue="#56B4E9", green="#009E73", yellow="#F0E442", blue="#0072B2", red="#D55E00", purple="#CC79A7", darkgreen="springgreen4", darkpurple="mediumorchid4")


df %>% group_by(Layer,Tree_species, sample)%>% 
  mutate(sum_per_group=sum(count))%>% ungroup()%>%
  #group_by(Phylum, Cluster,Tree_species) %>%
  #summarise("reads_count"=sum(counts)) %>%
  group_by(Layer,Tree_species, sample,primary_lifestyle) %>%
  mutate(perc=(count/sum_per_group)*100) %>% ungroup()%>%
  group_by(primary_lifestyle) %>%
  summarise( max=max(perc), mean=mean(perc), median=median(perc)) #%>%
  #slice_max(.,order_by = max,n=12)
  
  
#6x8  
#ITS2abundance <- abundance
#ITS1abundance<- abundance
#AMFabundance <- abundance
  
####################################### metaTranscriptomic ######################  
  
  Fungi_TPM <- read_delim("../../metatranscriptomique/ALL_SITES_Fungi_TPM.txt", "\t", escape_double = FALSE, col_names = TRUE, trim_ws = TRUE)
  treatment_metadata <- as.data.frame(read_delim("../../metatranscriptomique/metadata_all.txt", "\t", escape_double = FALSE, trim_ws = TRUE))
  colnames(treatment_metadata) <- c("sample", "samplenames", "Tree_number", "Layer", "Cluster","Tree_species" )
  
  abundance <- Fungi_TPM %>%
    dplyr::select(-starts_with("NR_"),-starts_with("Myc"),-ID,-ctg_length, -sum_contig_count) %>% 
    #filter(grepl('saprotroph',LifeStyle) ) %>%
    group_by(PD_Phylum,PD_Class, PD_Order,PD_Family, PD_Genus, LifeStyle) %>% summarise_if(is.numeric, sum) %>% pivot_longer(cols=-c(PD_Phylum,PD_Class, PD_Order,PD_Family, PD_Genus, LifeStyle), names_to="sample", values_to="tpm_counts") %>%
    left_join(treatment_metadata) %>% #left_join(ecology, by=c("PD_Genus"="Genus") ) %>%
    mutate(LifeStyle=case_when(!(is.na(LifeStyle)) ~ LifeStyle, TRUE ~ "Unidentified"))
  
#metaTabundance <- abundance
  
  toplifestyle <- abundance %>% #select( ID,primary_lifestyle, tmp_counts ) %>%
    #filter(ID %in% df$ID) %>% 
    group_by(LifeStyle) %>%
    summarise(val=sum(tpm_counts)) %>% 
    slice_max(val, n=12)
  
  df <- abundance %>% group_by(Layer,Tree_species, sample,LifeStyle )%>% 
    mutate(LifeStyle=case_when(LifeStyle %in% toplifestyle$LifeStyle ~ LifeStyle, TRUE ~ "Others"))%>%
    summarise("count"=sum(tpm_counts),.groups="keep")
  
  mycolors <-  mycolorsFullList %>% filter(Name %in% df$LifeStyle)
  mycolors <- mycolors %>% arrange(myLevel)
  myColorEco <- mycolors$color
  names(myColorEco) <- mycolors$Name
  df$LifeStyle<- factor(as.character(df$LifeStyle), levels=mycolors$Name)
  # 
  
  ggplot(data=df, aes(x=sample, y=count, fill=LifeStyle)) +
    geom_col() + facet_wrap(~ Layer*Tree_species, scales = "free_x") +
    scale_fill_manual(values=myColorEco,limits =names(myColorEco) )+
    #c("#999999",rev(c("#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7","springgreen4","mediumorchid4","cadetblue","lightpink","khaki","turquoise1")))) +
    theme(axis.text.x=element_text(angle = 55,hjust=1))
  #c(grey="#999999", orange="#E69F00", lightblue="#56B4E9", green="#009E73", yellow="#F0E442", blue="#0072B2", red="#D55E00", purple="#CC79A7", darkgreen="springgreen4", darkpurple="mediumorchid4")
  
  
 df %>% group_by(Layer,Tree_species, sample)%>% 
    mutate(sum_per_group=sum(count))%>% ungroup()%>%
    #group_by(Phylum, Cluster,Tree_species) %>%
    #summarise("reads_count"=sum(counts)) %>%
    group_by(Layer,Tree_species, sample,LifeStyle) %>%
    mutate(perc=(count/sum_per_group)*100) %>% ungroup()%>%
    group_by(LifeStyle) %>%
    summarise( max=max(perc), mean=mean(perc), median=median(perc)) #%>%
  #slice_max(.,order_by = max,n=12)
  
  
#######################################  
  
library("ggvenn")
 
ITS1 <- ITS1abundance #%>% filter(grepl('saprotroph',primary_lifestyle))
 
ITS2 <- ITS2abundance #%>% filter(grepl('saprotroph',primary_lifestyle))

AMF <- AMFabundance #%>% filter(grepl('saprotroph',primary_lifestyle)) 

metaB <- unique(c(unique(ITS1$Genus), unique(ITS2$Genus), unique(AMF$Genus)))

metaT <- abundance #%>% filter(grepl('saprotroph',LifeStyle) )
   
x <- list("ITS1"=unique(ITS1$Genus),
          
          "ITS2"=unique(ITS2$Genus),
          "18S"=unique(AMF$Genus),
          #"metabarcoding"=metaB,
          #"metatranscriptomic"= unique(metaT$Genus)
           "metatranscriptomic"= unique(metaT$PD_Genus)
          )

# Graphique par défaut
ggvenn(x,
       fill_color = c("beige", "lightblue","lightgreen","red"
                     ),
       stroke_size = 0.5,
       set_name_size = 4)

##########################

overall <- sum(ITS1abundance$reads_count)
ITS1 <- ITS1abundance %>% group_by(Genus) %>%
  filter(!(grepl('saprotroph',primary_lifestyle))) %>%
  filter(!(grepl('mycorrhizal',primary_lifestyle))) %>%
             mutate(perctmp=(reads_count/overall)*100) %>%
             select(Genus,perctmp) %>% group_by(Genus)%>% summarise(perc=sum(perctmp))
ITS1$cond <- "ITS1"


overall <- sum(ITS2abundance$reads_count)
ITS2 <- ITS2abundance %>% group_by(Genus) %>%
  filter(!(grepl('mycorrhizal',primary_lifestyle))) %>%
  filter(!(grepl('saprotroph',primary_lifestyle))) %>%
  mutate(perctmp=(reads_count/overall)*100) %>%
  select(Genus,perctmp) %>% group_by(Genus)%>% summarise(perc=sum(perctmp))

ITS2$cond <- "ITS2"

overall <- sum(AMFabundance$reads_count)
AMF <- AMFabundance %>% group_by(Genus) %>%
  filter(!(grepl('saprotroph',primary_lifestyle))) %>%
  filter(!(grepl('mycorrhizal',primary_lifestyle))) %>%
  mutate(perctmp=(reads_count/overall)*100) %>%
  select(Genus,perctmp) %>% group_by(Genus)%>% summarise(perc=sum(perctmp))

AMF$cond <- "18S"

overall <- sum(metaTabundance$tpm_counts)
metaT <- metaTabundance %>% group_by(PD_Genus) %>% 
  filter(!(grepl('saprotroph',LifeStyle))) %>%
  filter(!(grepl('mycorrhizal',LifeStyle))) %>%
            mutate(perctmp=(tpm_counts /overall)*100) %>%
           select(PD_Genus,perctmp) %>% group_by(PD_Genus)%>% summarise(perc=sum(perctmp))

colnames(metaT) <- c( "Genus","perc") 
metaT$cond <- "metatranscriptomic"

df <- rbind(ITS1,ITS2,AMF,metaT)

library("pheatmap")

topGenus <- df %>% filter(!is.na(Genus)) %>%
  filter(perc > 0.1)

m <- df %>% filter(!is.na(Genus)) %>%
  #filter(perc < 0.1)%>%
  filter((Genus %in% topGenus$Genus )) %>%
  #mutate(perc= case_when(perc ==0 ~ perc, TRUE ~ 1 ))%>%
  mutate(perc= case_when(perc < 0.5 ~ perc, TRUE ~ 0.5 ))%>%
  pivot_wider(names_from = cond, values_from = perc, values_fill=0)
  
a <- m
m$Genus <- NULL
rownames(m) <- a$Genus

toplifestyle <- abundance %>% #select( ID,primary_lifestyle, tmp_counts ) %>%
  #filter(ID %in% df$ID) %>% 
  group_by(LifeStyle) %>%
  summarise(val=sum(tpm_counts)) %>% 
  slice_max(val, n=15)

annOTU <- as.data.frame(rownames(m))
colnames(annOTU) <- "Genus"
annOTU <- annOTU %>% left_join(ecology) %>%
  select(primary_lifestyle) %>% 
  mutate(primary_lifestyle= case_when(primary_lifestyle != "" ~ primary_lifestyle, TRUE ~ "Unidentified" )) #%>%
  #mutate(primary_lifestyle=case_when(primary_lifestyle %in% toplifestyle$LifeStyle ~ primary_lifestyle, TRUE ~ "Others"))
  

rownames(annOTU) <- rownames(m)





mycolors <-  mycolorsFullList %>% filter(Name %in% annOTU$primary_lifestyle)
mycolors <- mycolors %>% arrange(myLevel)
annotation_colors_ecology <- mycolors$color
names(annotation_colors_ecology) <- mycolors$Name

annotation_colors <- list("primary_lifestyle"=annotation_colors_ecology)

pheatmap(m, #color=c("white","black"),
         colorRampPalette(c("white",brewer.pal(9, "Blues"),"black"))(100),
         annotation_row = annOTU,

         fontsize_row=8,
         annotation_colors = annotation_colors
         )

#?pheatmap
