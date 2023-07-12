library("readr")
library("plyr");library("dplyr");library("ggplot2");library("tidyr");library("tidyverse")
library("RColorBrewer");library("ggpubr");library("vegan")


treatment_metadata <- as.data.frame(read_delim("../0_metadata/metadata_all.txt", "\t", escape_double = FALSE, trim_ws = TRUE))
treatment_metadata$layer_tree <- paste( treatment_metadata$layer,treatment_metadata$tree_species, sep="_")
treatment_metadata$Layer <- treatment_metadata$layer
treatment_metadata <- treatment_metadata %>% mutate(layer=case_when(layer != "OS" ~ layer, TRUE ~ "OS: 0 - 5cm")) 
treatment_metadata <- treatment_metadata %>% mutate(layer=case_when(layer != "OM" ~ layer, TRUE ~ "OM: 5 - 25cm"))
treatment_metadata$layer <- factor(treatment_metadata$layer, levels=c("OS: 0 - 5cm","OM: 5 - 25cm" ))


Fungi_TPM <- read_delim("../1_data/ALL_SITES_Fungi_TPM.txt", "\t", escape_double = FALSE, col_names = TRUE, trim_ws = TRUE) %>%
  filter(tax_kingdom %in% "Fungi", PD_Name %notin% c("Batde5","Crypa2","Pospl1")) %>% 
  select(cols=-c("MEROPS_coverage", "MEROPS_evalue","MEROPS_score" )) 
  
  
Fungi_TPM <-  Fungi_TPM %>% mutate(species= paste(PD_Genus,"aff.",PD_Species,sep=" "))
Fungi_TPM<-  Fungi_TPM %>% mutate(species= case_when(species!="Amanita aff. cecilae" ~ species, TRUE ~  "Amanita aff. griseofolia")) %>%
                            mutate(species= case_when(species!="Tuber aff. borchii" ~ species, TRUE ~  "ATuber zhongdianense")) %>%
                            mutate(species= case_when(species!="Russula aff. senecis" ~ species, TRUE ~  "Russula punctipes")) 
  



Fungi_TPM <- Fungi_TPM %>%  mutate(LifeStyle = case_when(LifeStyle %notin% c("dung_saprotroph") ~ LifeStyle, TRUE ~ "soil_saprotroph")) %>%
                            mutate(LifeStyle = case_when(LifeStyle %notin% c("nectar/tap_saprotroph","unspecified_saprotroph","pollen_saprotroph") ~ 
                                                           LifeStyle, TRUE ~ "other_saprotroph")) %>%
                            mutate(LifeStyle = case_when(!is.na(LifeStyle) ~ LifeStyle, TRUE ~ "unassigned")) %>%
                          mutate(LifeStyle = case_when(LifeStyle %in% c("plant_pathogen",
                                                                        "wood_saprotroph", "litter_saprotroph","soil_saprotroph","other_saprotroph",
                                                                        "ectomycorrhizal", "arbuscular_mycorrhizal",
                                                                        "unassigned") ~ LifeStyle, TRUE ~ "others"))



Fungi_ratio <- Fungi_TPM %>%
  dplyr::select(-starts_with("NR_"),-starts_with("Myc"),-ID,-ctg_length, -sum_contig_count) %>% 
  group_by(PD_Species, PD_Genus, PD_Family, PD_Order, PD_Class, PD_Phylum,LifeStyle) %>% summarise_if(is.numeric, sum) %>% 
  pivot_longer(cols=-c(PD_Species, PD_Genus, PD_Family, PD_Order, PD_Class, PD_Phylum,LifeStyle), names_to="samples", values_to="tpm_counts") %>%
  left_join(treatment_metadata) %>%
  group_by(PD_Species, PD_Genus, PD_Family, PD_Order, PD_Class, PD_Phylum,LifeStyle, layer_tree) %>% summarise(tpm_sum=sum(tpm_counts)) %>%
  pivot_wider(names_from = layer_tree, values_from = tpm_sum, values_fill = 0)

Fungi_ratio <- Fungi_TPM %>%
  dplyr::select(-starts_with("NR_"),-starts_with("Myc"),-ID,-ctg_length, -sum_contig_count) %>%
  filter(PD_Species %in% c("borchii","undulatus","glaucopus","triplex",
                           "muscaria","haematopus","badius","pura","granulatus","saniosus", "sp",
                           "camargensis"), PD_Genus %in% c("Tuber","Basidioascus","Phlegmacium",
                                                           "Geastrum","Amanita","Mycena","Xerocomus",
                                                           "Umbelopsis","Elaphomyces","Cortinarius",
                                                           "Linnemannia")) %>%
  group_by(contig,PD_Species, PD_Genus) %>% summarise_if(is.numeric, sum) %>%
  pivot_longer(cols=-c(contig,PD_Species, PD_Genus), names_to="samples", values_to="tpm_counts") %>%
  left_join(treatment_metadata) %>%
  group_by(PD_Species, PD_Genus,contig, layer_tree) %>% summarise(tpm_sum=sum(tpm_counts))

Fungi_ratio <- Fungi_TPM %>%
  dplyr::select(-starts_with("NR_"),-starts_with("Myc"),-ID,-ctg_length, -sum_contig_count) %>%
  filter(contig %in% c("k141_503747","k141_3748920","k141_2562615","k141_3073207","k141_3703555",
                       "k141_1319161","k141_326632","k141_7646129", "k141_56981", "k141_7162408",
                       "k141_5541758","k141_8777435"))

#write.table(Fungi_ratio,file="../taxonomy_lisadd.tsv", sep="\t", row.names = F)

Fungi_ratio <- Fungi_TPM %>%
  dplyr::select(-starts_with("NR_"),-starts_with("Myc"),-ID,-ctg_length, -sum_contig_count) %>% 
  group_by(LifeStyle) %>% summarise_if(is.numeric, sum) %>% pivot_longer(cols=-c(LifeStyle), names_to="samples", values_to="tpm_counts") %>%
  left_join(treatment_metadata)  %>% 
  filter(!grepl("score", samples))%>% 
  filter(!grepl("evalue", samples))%>% 
  filter(!grepl("coverage", samples))

tmp <- Fungi_ratio %>% ungroup %>% group_by(samples) %>% mutate(tmp_perc=(tpm_counts/sum(tpm_counts))*100) %>% ungroup %>%
  group_by(LifeStyle, samples) %>% summarise(tpm_group=sum(tmp_perc)) %>% group_by(LifeStyle) %>%
  summarise(max=max(tpm_group), mean=mean(tpm_group), min=min(tpm_group)) 

top10Eco <- Fungi_TPM %>%
  dplyr::select(-starts_with("NR_"),-starts_with("Myc"),-ID,-ctg_length, -sum_contig_count) %>% 
  group_by(LifeStyle) %>% summarise_if(is.numeric, sum) %>% pivot_longer(cols=-c(LifeStyle), names_to="samples", values_to="tpm_counts")%>%
  group_by(LifeStyle) %>% summarise(max=max(tpm_counts), mean=mean(tpm_counts), median=median(tpm_counts)) %>% 
  slice_max(.,order_by = mean,n=10)



###### Lifestyle #######
mycolors <- GetMyColor(Fungi_TPM$LifeStyle)
Fungi_TPM$LifeStyle <- factor(Fungi_TPM$LifeStyle, levels=names(mycolors)) 

Fungi_ratio <- get.Fungi_ratio(Fungi_TPM, WhereToTest="LifeStyle",SelectMetadataValue="all", columnMetadataValue="all") 

### barplot per sample ###
ggplot(data=Fungi_ratio, aes(x=samples, y=tpm_counts, fill=LifeStyle)) +
  geom_col() + facet_wrap(. ~ layer*tree_species, scales = "free_x") + 
  scale_fill_manual(values=mycolors)+
  theme(axis.text.x=element_text(angle = 55,hjust=1))
#8*10

### boxplot with Tukey test on Lifestyle ###
# prepare data
Fungi_ratio <- Fungi_ratio %>% 
              mutate(LifeStyle=case_when(!grepl('saprotroph',LifeStyle) ~LifeStyle, TRUE ~ "SAP"  )) %>% # Group saprotrophs together
              mutate(LifeStyle=case_when(LifeStyle !="ectomycorrhizal" ~LifeStyle, TRUE ~ "EM"  )) %>% # rename ectomycorrhizae into EM
              mutate(LifeStyle = case_when(LifeStyle %in% c("EM","SAP","plant_pathogen" ) ~ LifeStyle, TRUE ~ "Others"))  %>%   
              group_by(LifeStyle, samples, tree_species, Layer, layer, layer_tree) %>% summarise(tpm_counts=sum(tpm_counts)) %>% 
              mutate(tree_association_ecology=paste(Layer, tree_species,  LifeStyle)) # create an unique tag for the combinaison tree -layer-ecology
  
library(multcompView)
#anova
model <- aov(tpm_counts~tree_association_ecology, data=Fungi_ratio)
summary(model)
#tukey test
tukey.res <-TukeyHSD(model, conf.level=.95)
#get factors possibility
tmp <-Fungi_ratio %>% ungroup()%>% select(tree_association_ecology,layer_tree,layer, tree_species,LifeStyle ) %>% unique()
# get label (use function in annex file)
LABELS <- generate_label_df(tukey.res , "tree_association_ecology")
LABELS<-  LABELS %>% rename("tree_association_ecology"=treatment) %>% left_join(tmp)

mycolors <- GetMyColor(Fungi_ratio$LifeStyle)
tmp <- GetMyColor(Fungi_ratio$tree_species)
mycolorLabel=rep(tmp,each=4)


ggplot(data=Fungi_ratio, aes(x=tree_association_ecology, y=tpm_counts, colour=LifeStyle)) +
  geom_boxplot(outlier.size = 0.5) + 
  #ggtitle(catToTest)+
  scale_colour_manual(values=mycolors)+
  scale_y_continuous(labels = function(x) format(x, scientific = FALSE))+
  geom_jitter(size=0.5)+
  geom_text(
    data = LABELS,
    aes(
      x = tree_association_ecology, 
      y = max(Fungi_ratio$tpm_counts)+25000, 
      label = Letters
    ))  +  facet_wrap(. ~ layer, scales = "free_x") +
  theme_bw()+
  #theme(legend.position="bottom")+
  theme(axis.text.x=element_text(angle = 55,hjust=1, colour = mycolorLabel), axis.title.x=element_blank())
#8*12






### plot main graph - Fig.4a
Fungi_ratio <- Fungi_TPM %>%
  dplyr::select(-starts_with("NR_"),-starts_with("Myc"),-ID,-ctg_length, -sum_contig_count) %>% 
  group_by(LifeStyle) %>% summarise_if(is.numeric, sum) %>% pivot_longer(cols=-c(LifeStyle), names_to="samples", values_to="tpm_counts") %>%
  left_join(treatment_metadata) %>%
  group_by(LifeStyle, layer, layer_tree) %>% summarise(tpm_mean=mean(tpm_counts))

mycolors <- GetMyColor(Fungi_ratio$LifeStyle)
ggplot(data=Fungi_ratio, aes(x=layer_tree, y=tpm_mean, fill=LifeStyle)) +
  geom_col() + #facet_wrap(. ~ layer*tree_species, scales = "free_x") + 
  facet_wrap(. ~ layer, scales = "free_x") + 
  scale_fill_manual(values=mycolors)+
  #c("#999999",rev(c("#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7","springgreen4","mediumorchid4","cadetblue","lightpink","khaki","turquoise1")))) +
  theme(axis.text.x=element_text(angle = 55,hjust=1), axis.title.x=element_blank())
#5*5


#########################

##-------------------------
## -------- KOG -----------
##-------------------------

KOG_meta <- read_delim("../0_metadata/KOG_list.txt", "\t", col_names = TRUE)
KOG_meta$'KOG class' <- paste(KOG_meta$Letter_code, ": ",KOG_meta$KOG_class )

####### KOG group ############
WhereToTest <- "KOG_group"

Fungi_TPM_filtered_Annot <-Split_Annot_In_Multiple_Rows(Fungi_TPM, WhereToTest)

WhereToTest2 <- "LifeStyle"

Fungi_ratio <- get.Fungi_ratio.2groupments(Fungi_TPM_filtered_Annot, WhereToTest="KOG_group",WhereToTest2="LifeStyle",SelectMetadataValue="all", columnMetadataValue="all")

Fungi_ratio <-Fungi_ratio %>%
  mutate(KOG_group=case_when(KOG_group %in% c("CELLULAR PROCESSES AND SIGNALING", "INFORMATION STORAGE AND PROCESSING", "METABOLISM") ~ KOG_group, TRUE ~ "Unknown"  )) %>%
  mutate(KOG_group=case_when(KOG_group != "CELLULAR PROCESSES AND SIGNALING" ~ KOG_group,TRUE ~ "Cellular P.&Si.")) %>%
  mutate(KOG_group=case_when(KOG_group != "INFORMATION STORAGE AND PROCESSING" ~ KOG_group,TRUE ~ "Info. St.&P.")) %>%
  mutate(KOG_group=case_when(KOG_group != "METABOLISM" ~ KOG_group,TRUE ~ "Metabolism"))%>%
  group_by(KOG_group, LifeStyle, samples) %>% summarise(tpm=sum(tpm_counts))
colnames(Fungi_ratio)[colnames(Fungi_ratio) == "tpm"] <- "tpm_counts"

dfstatTPM <- get_dfstatTPM_LifeStyle(Fungi_ratio, WhereToTest)
dfstatPERC <- get_dfstatPERC_LifeStyle(Fungi_ratio, WhereToTest)

GroupToTest <- unique(Fungi_ratio$KOG_group) 
WhereToTest <- "KOG_group"
# GroupToTest <- GroupToTest[!GroupToTest %in% c("0", NA)]
# byCond= "Layer"
# condValues= c("OS","OM")
#t.res <- t.test_paired_Layer(Fungi_ratio, GroupToTest, WhereToTest)
# test <- Fungi_ratio %>%
#        group_by(KOG_group, LifeStyle, samples) %>% summarise(tpm=sum(tpm_counts)) %>% pivot_wider(names_from = samples, values_from = tpm)
###### plot 
dfTPM <- Fungi_ratio %>%   left_join(treatment_metadata) %>% 
  group_by(KOG_group, LifeStyle,layer,tree_species, layer_tree) %>% summarise(mean_tpm=mean(tpm_counts))
#dfTPM <- dfTPM %>% left_join(sumlayertree) %>% mutate(perc=(tpm/sommeLayerTree)*100)

mycolors <- GetMyColor(dfTPM$LifeStyle)
dfTPM$LifeStyle <- factor(dfTPM$LifeStyle, levels=names(mycolors)) 

ggplot(data=dfTPM, aes(x=layer_tree, y=mean_tpm, fill=LifeStyle  )) +
  geom_col() + scale_fill_manual(values=mycolors)+
  theme_bw()+
  theme(axis.text.x=element_text(angle = 55,hjust=1), axis.title.x=element_blank())+
  #ylim(0, 2000000)+
  facet_wrap(. ~ KOG_group*layer, scales = "free_x", nrow=2) 
#size 7*8

###########

####### KOG class ############

WhereToTest <- "KOG_class"

Fungi_TPM_filtered_Annot <-Split_Annot_In_Multiple_Rows(Fungi_TPM, WhereToTest)
Fungi_TPM_filtered_Annot <- Fungi_TPM_filtered_Annot %>% filter(KOG_group %in% c("CELLULAR PROCESSES AND SIGNALING", "INFORMATION STORAGE AND PROCESSING", "METABOLISM"))
Fungi_TPM_filtered_Annot$KOG_class[Fungi_TPM_filtered_Annot$KOG_class == "Posttranslational modification, protein turnover, chaperones"] <- "Post-translational modification, protein turnover and chaperones"

###### Total KOG_class ######

Fungi_ratio <- get.Fungi_ratio(Fungi_TPM_filtered_Annot, WhereToTest,SelectMetadataValue="all", columnMetadataValue="all")
# Fungi_ratio <- get.Fungi_ratio(Fungi_TPM, WhereToTest,SelectMetadataValue="Quercus", columnMetadataValue="tree_species")

dfstatTPM <- get_dfstatTPM(Fungi_ratio, WhereToTest)
dfstatPERC <- get_dfstatPERC(Fungi_ratio, WhereToTest)
dfstatPERC <- merge(dfstatPERC,KOG_meta)
dfstatPERC <- dfstatPERC %>% select(`KOG class`, max, min, mean, median)
#############
###### Total by Lifestyle KOG_class ######

WhereToTest2 <- "LifeStyle"

Fungi_ratio <- get.Fungi_ratio.2groupments(Fungi_TPM_filtered_Annot, WhereToTest, WhereToTest2="LifeStyle",SelectMetadataValue="all", columnMetadataValue="all")

Fungi_ratio <-Fungi_ratio %>%
   mutate(LifeStyle=case_when(!grepl('saprotroph',LifeStyle) ~LifeStyle, TRUE ~ "saprotroph"  )) %>%
   mutate(LifeStyle = case_when(LifeStyle %in% c("ectomycorrhizal","saprotroph" ) ~ LifeStyle, TRUE ~ "Others")) %>%
  group_by(KOG_class, LifeStyle, samples) %>% summarise(tpm=sum(tpm_counts))
colnames(Fungi_ratio)[colnames(Fungi_ratio) == "tpm"] <- "tpm_counts"

dfstatTPM <- get_dfstatTPM_LifeStyle(Fungi_ratio, WhereToTest)
dfstatPERC <- get_dfstatPERC_LifeStyle(Fungi_ratio, WhereToTest)

###### plot ###
dfTPM <- Fungi_ratio %>%   left_join(treatment_metadata) %>%
  group_by(KOG_class, LifeStyle,layer,tree_species, layer_tree) %>% summarise(mean_tpm=mean(tpm_counts)) %>% join(KOG_meta)
#dfTPM <- dfTPM %>% left_join(sumlayertree) %>% mutate(perc=(tpm/sommeLayerTree)*100)

mycolors <- GetMyColor(dfTPM$LifeStyle)
dfTPM$LifeStyle <- factor(dfTPM$LifeStyle, levels=names(mycolors)) 

ggplot(data=dfTPM, aes(x=Letter_code, y=mean_tpm, fill=LifeStyle , alpha=`KOG class` )) +
  geom_col() + scale_fill_manual(values=mycolors)+
  scale_alpha_manual(values=rep(1,times=23)) +
  theme_bw()+
  #theme(legend.position="bottom")+
  #ylim(0, 2000000)+
  facet_wrap(. ~ layer*tree_species, scales = "free_x") #+ theme(axis.text.x=element_text(angle = 55,hjust=1)) 
#size 7*18
#### transform to heatmap ?

###################
mycolors <- GetMyColor(Fungi_TPM_filtered_Annot$LifeStyle)
Fungi_TPM_filtered_Annot$LifeStyle <- factor(Fungi_TPM_filtered_Annot$LifeStyle, levels=names(mycolors)) 

GroupToTest <- unique(Fungi_TPM_filtered_Annot$KOG_class) 
catToTest <- "Translation, ribosomal structure and biogenesis"  
for(catToTest in GroupToTest){
Fungi_ratio <- Fungi_TPM_filtered_Annot %>%
  dplyr::select(-starts_with("NR_"),-starts_with("Myc"),-ID,-ctg_length, -sum_contig_count) %>%
  filter(KOG_class == catToTest) %>%
  group_by(LifeStyle) %>% summarise_if(is.numeric, sum) %>% pivot_longer(cols=-c(LifeStyle), names_to="samples", values_to="tpm_counts") %>%
  left_join(treatment_metadata) %>%
  #group_by(LifeStyle, layer, layer_tree) %>% summarise(tpm_median=median(tpm_counts))
  group_by(LifeStyle, layer, layer_tree) %>% summarise(tpm_mean=mean(tpm_counts))
  
  
nom <- paste(catToTest,".pdf", sep="")
nom <- gsub("/", "-", nom)
print(nom)

p <- #ggplot(data=Fungi_ratio, aes(x=layer_tree, y=tpm_median, fill=LifeStyle)) +
  ggplot(data=Fungi_ratio, aes(x=layer_tree, y=tpm_mean, fill=LifeStyle)) +
  geom_col() + #facet_wrap(. ~ layer*tree_species, scales = "free_x") +
  facet_wrap(. ~ layer, scales = "free_x") +
  scale_fill_manual(values=mycolors)+
  #c("#999999",rev(c("#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7","springgreen4","mediumorchid4","cadetblue","lightpink","khaki","turquoise1")))) +
  theme(axis.text.x=element_text(angle = 55,hjust=1), axis.title.x=element_blank())+
  ggtitle(catToTest)

pdf(file=nom,width=5, height=5)
print(p)
dev.off()
#5*5

}

#########################

##-------------------------
## --------  Ions/lipids/Carbohydrate -----------
##-------------------------


####### get annotation #####
Fungi_TPM$ID2 <- paste(Fungi_TPM$contig, Fungi_TPM$NR_score , sep = "__")
rownames(Fungi_TPM) <- Fungi_TPM$ID2


tmp_lipases <- Fungi_TPM %>% filter(KOG_group %in% "METABOLISM" &
                                      KOG_class %in% "Lipid transport and metabolism" &
                                      grepl('Lipase|lipase|lase|dehydrogenase|dase|esterase', KOG_defline)) %>%
  mutate(UA="Lipases") %>%
  select(ID2,UA)

tmp_transporter <- Fungi_TPM %>% filter(KOG_group %in% "METABOLISM" &
                                          KOG_class %in% c("Amino acid transport and metabolism",
                                                           "Carbohydrate transport and metabolism", 
                                                           "Inorganic ion transport and metabolism",
                                                           "Lipid transport and metabolism", 
                                                           "Nucleotide transport and metabolism") &
                                          grepl('Transport|transport|Channel|channel|Permease|permease', KOG_defline)) %>%
  mutate(UA="transporter") 

tmp_nitrogen <- tmp_transporter %>% filter(grepl('Ammon|Nitr|nitr|ammon', KOG_defline)) %>%
  mutate(UA_class="Nitrogen-related transporters") %>%
  select(ID2,UA_class)

tmp_nitrogen2 <-  tmp_transporter %>% filter( KOG_class %in% c("Amino acid transport and metabolism",
                                                               "Nucleotide transport and metabolism")) %>%
  mutate(UA_class="Nitrogen-related transporters") %>% select(ID2, UA_class)


tmp_inor <- tmp_transporter %>% filter( KOG_class %in% "Inorganic ion transport and metabolism") %>%
  mutate(UA_class="Inorganic ion transporters") %>%
  select(ID2,UA_class)

tmp_carb <-  tmp_transporter %>% filter(KOG_class %in% c("Carbohydrate transport and metabolism",
                                                         "Lipid transport and metabolism")) %>%
  mutate(UA_class="Carbon-related transporters") %>%
  select(ID2,UA_class)
tmp_transporter <- tmp_transporter %>% select(ID2,UA) 


UA_class_d <- rbind(tmp_nitrogen,tmp_nitrogen2,tmp_inor, tmp_carb)
UA_d <- rbind(tmp_lipases,tmp_transporter)

Fungi_TPM$UA <- NA
Fungi_TPM$UA_class <- NA
Fungi_TPM <- Fungi_TPM %>% mutate(UA= ifelse(ID2 %in% UA_d$ID2, UA, UA_d$UA))
Fungi_TPM <- Fungi_TPM %>% mutate(UA_class= ifelse(ID2 %in% UA_class_d$ID2, UA_class, UA_class_d$UA_class))


##------- boxplots ----------
##----------------------------

kog="Nitrogen-related transporters"
kog="Carbon-related transporters"
kog="Inorganic ion transporters"


Fungi_ratio <- Fungi_TPM %>%
  dplyr::select(-starts_with("NR_"),-starts_with("Myc"),-ID,-ctg_length, -sum_contig_count) %>% 
  #filter(grepl('mycorrh',LifeStyle) ) %>%
  mutate(LifeStyle=case_when(!grepl('saprotroph',LifeStyle) ~LifeStyle, TRUE ~ "SAP"  )) %>%
  mutate(LifeStyle=case_when(LifeStyle != "ectomycorrhizal" ~ LifeStyle, TRUE ~ "EM"  )) %>%
  filter(LifeStyle %in% c("SAP","EM") ) %>%
  filter(UA_class == kog ) %>% 
  group_by(UA_class,LifeStyle) %>% summarise_if(is.numeric, sum) %>% pivot_longer(cols=-c(UA_class,LifeStyle), names_to="samples", values_to="tpm_counts") %>%
  left_join(treatment_metadata)  # %>% left_join(sumlayertree) %>% mutate(perc=(tpm_counts/sommeLayerTree)*100)

mycolors <- GetMyColor(Fungi_ratio$LifeStyle)
Fungi_ratio$LifeStyle <- factor(Fungi_ratio$LifeStyle, levels=names(mycolors)) 

mycolors <- GetMyColor(Fungi_ratio$tree_species)
ggplot(data = Fungi_ratio, aes(x=LifeStyle, y=tpm_counts, fill=tree_species))+
  geom_boxplot()+
  theme_bw()+
  scale_fill_manual(values=mycolors) +
  facet_wrap(. ~ layer, nrow=1)+
  labs(title=kog) +
  ylim(0,4e5)+
  theme(axis.title.x=element_blank()
        #axis.title.y =element_blank()
  )

#5*4

# 
# 
# 
# #######################
# annotSample <- treatment_metadata %>% select(layer, cluster_number, tree_species)
# rownames(annotSample) <- treatment_metadata$samples
# myColorAnnot <- list("tree_species"= GetMyColor(annotSample$tree_species),
#                      "cluster_number"=GetMyColor(annotSample$cluster_number),
#                      "layer"=GetMyColor(annotSample$layer))
# 
# 
# ecology="all"
# kog="Nitrogen-related transporters"
# 
# Fungi_ratio <- Fungi_TPM %>% dplyr::select(-starts_with("NR_"),-starts_with("Myc"),-ID,-ctg_length, -sum_contig_count) %>% 
#   filter(UA_class == kog ) %>% 
#   #filter(KOG_class %in% "Coenzyme transport and metabolism" ) %>%
#   #filter(grepl(ecology ,LifeStyle) ) %>%
#   #mutate(species= paste(PD_Genus,PD_Species,sep="_")) %>%
#   mutate(species= PD_Order) %>%
#   #filter(!(tax_sciname %in% c("56484", "5306", "2461416", "1898205", "*"))) %>%
#   #mutate(tax_sciname = str_trunc(tax_sciname, 30, "right")) %>%
#   group_by(species) %>% summarise_if(is.numeric, sum) %>% pivot_longer(cols=-c(species), names_to="samples", values_to="tpm") %>%
#   left_join(treatment_metadata) %>% group_by(species, cluster_number) %>% summarise(tmp_counts=sum(tpm), .groups = "keep") %>%
#   pivot_wider(names_from  = cluster_number, values_from = "tmp_counts")
# 
# 
# sp_tmp <- Fungi_ratio %>% pivot_longer(cols=-c(species), names_to="cluster_number", values_to="tpm_counts") %>% 
#   group_by(species) %>% summarise_if(is.numeric, sum) %>%
#   filter(tpm_counts < 200)
# sp_preval <- Fungi_ratio %>% pivot_longer(cols=-c(species), names_to="cluster_number", values_to="tpm_counts") %>% 
#   mutate(nb=0)%>%
#   mutate(nb=case_when(tpm_counts < 1 ~  nb, TRUE ~ 1)) %>%
#   group_by(species) %>% summarise(preval=sum(nb)) %>%
#   filter(preval < 4)
# 
# 
# Fungi_ratio <- Fungi_ratio %>% filter(species %notin% sp_tmp$species ) %>% 
#   filter(species %notin% sp_preval$species ) %>%
#   filter(!is.na(species))
# rownames(Fungi_ratio) <- Fungi_ratio$species
# 
# 
# m <- Fungi_ratio %>% ungroup() %>% select(-species)
# rownames(m) <- Fungi_ratio$species
# 
# nom <- paste(ecology, "__", kog,".pdf", sep="")
# nom <- gsub("/", "-", nom)
# 
# nom2 <- paste(ecology,":", kog)
# library(pheatmap)
# 
# p <- pheatmap(log(m+1),
#               #color= c(rev(colorRampPalette(brewer.pal(9, "Blues"))(100)), colorRampPalette(brewer.pal(9, "YlOrRd"))(100), "black"),
#               main=nom2,
#               cellwidth=2,
#               cellheight = 4,
#               fontsize_col=2,
#               fontsize_row=4,
#               annotation_col = annotSample,
#               #annotation_row = annotCAZymes,
#               annotation_colors = myColorAnnot,
#               color= c("white",(colorRampPalette(brewer.pal(9, "OrRd"))(100)),"black"))
# 
# p 
# #h <- 4 + length(Fungi_ratio$species)/40
# h <- 4 + length(Fungi_ratio$species)/20



#########################

##-------------------------------
## -------- KOG class  Heatmaps -----------
##-------------------------------
library(pheatmap)
treatment_metadata$Layer_cluster <- paste(treatment_metadata$Layer,treatment_metadata$cluster_number, sep = " ")

# annotSample <- treatment_metadata %>% select(layer, tree_species,layer_tree) %>% unique() 
# rownames(annotSample) <- annotSample$layer_tree
# annotSample$layer_tree <- NULL

annotSample <- treatment_metadata %>% select(layer, tree_species,Layer_cluster) %>% unique() 
rownames(annotSample) <- annotSample$Layer_cluster
annotSample$Layer_cluster <- NULL

WhereToTest <- "KOG_class"
Fungi_TPM_filtered_Annot <-Split_Annot_In_Multiple_Rows(Fungi_TPM, WhereToTest)
Fungi_TPM_filtered_Annot <- Fungi_TPM_filtered_Annot %>% filter(KOG_group %in% c("CELLULAR PROCESSES AND SIGNALING", "INFORMATION STORAGE AND PROCESSING", "METABOLISM"))
Fungi_TPM_filtered_Annot$KOG_class[Fungi_TPM_filtered_Annot$KOG_class == "Posttranslational modification, protein turnover, chaperones"] <- "Post-translational modification, protein turnover and chaperones"

# 
# myColorAnnot <- list("tree_species"= GetMyColor(annotSample$tree_species),
#                      "layer"=GetMyColor(annotSample$layer))

annotSpecies <- Fungi_TPM_filtered_Annot %>% mutate(LifeStyle=case_when(species!="Unclassified aff. sp" ~ LifeStyle, TRUE ~ "unassigned") ) %>%
  select(species, LifeStyle) %>% unique() 

#mycolors_tmp <- GetMyColor(annotSpecies$LifeStyle)
#library(scales)
#show_col(mycolors_tmp)

ecology=""
kog="Lipid transport and metabolism"
GroupToTest <- unique(Fungi_TPM_filtered_Annot$KOG_class)



for(kog in GroupToTest){
  
  Fungi_ratio <- Fungi_TPM_filtered_Annot %>% dplyr::select(-starts_with("NR_"),-starts_with("Myc"),-ID,-ctg_length, -sum_contig_count) %>% 
    filter( grepl(kog, KOG_class )) %>% 
    group_by(species) %>% summarise_if(is.numeric, sum) %>% pivot_longer(cols=-c(species), names_to="samples", values_to="tpm") %>%
    #left_join(treatment_metadata) %>% group_by(species, layer_tree) %>% summarise(tmp_counts=median(tpm), .groups = "keep") %>%
    left_join(treatment_metadata) %>% group_by(species, Layer_cluster ) %>% summarise(tmp_counts=median(tpm), .groups = "keep") %>%
    #pivot_wider(names_from  = layer_tree, values_from = "tmp_counts")
    pivot_wider(names_from  = Layer_cluster, values_from = "tmp_counts")
  
  sp_tmp <- Fungi_ratio %>% #pivot_longer(cols=-c(species), names_to="layer_tree", values_to="tpm_counts") %>% 
    pivot_longer(cols=-c(species), names_to="Layer_cluster", values_to="tpm_counts") %>% 
    #group_by(layer_tree) %>% 
    group_by(Layer_cluster) %>%
    mutate(sumLayerTree=sum(tpm_counts)) %>% ungroup() %>%
    mutate(observation_perc=(tpm_counts/sumLayerTree)*100  ) %>% 
    #arrange(desc(observation_perc)) %>%
    #group_by(layer_tree) %>% 
    group_by(Layer_cluster) %>%
    #dplyr::arrange(layer_tree, observation_perc) %>% 
    dplyr::arrange(Layer_cluster, observation_perc) %>% 
    dplyr::mutate(cumul_observation_perc=cumsum(observation_perc)) %>%
    #dplyr::arrange(layer_tree, desc(observation_perc)) %>% 
    dplyr::arrange(Layer_cluster, desc(observation_perc)) %>% 
    dplyr::mutate(cumul_observation_perc_desc=cumsum(observation_perc)) %>%
    filter(cumul_observation_perc >= 50)
  
  sp_preval <- Fungi_ratio %>% #pivot_longer(cols=-c(species), names_to="layer_tree", values_to="tpm_counts") %>% 
    pivot_longer(cols=-c(species), names_to="Layer_cluster", values_to="tpm_counts") %>% 
    mutate(nb=0)%>%
    mutate(nb=case_when(tpm_counts < 1 ~  nb, TRUE ~ 1)) %>%
    group_by(species) %>% summarise(preval=sum(nb)) %>%
    filter(preval < 4)
  
  
  Fungi_ratio <- Fungi_ratio %>% filter(species %in% sp_tmp$species ) %>% 
    filter(species %notin% sp_preval$species )
  rownames(Fungi_ratio) <- Fungi_ratio$species
  
  
  m <- Fungi_ratio %>% ungroup() %>% select(-species)
  rownames(m) <- Fungi_ratio$species
  
  #nom <- paste(ecology, "__", kog,".pdf", sep="")
  nom <- paste( kog,".pdf", sep="")
  nom <- gsub("/", "-", nom)
  
  #nom2 <- paste(ecology,":", kog)
  nom2 <- paste(kog)
  
  annotSpeciesSmall <- annotSpecies %>% filter(species %in% Fungi_ratio$species) 
  tmp <- annotSpeciesSmall
  annotSpeciesSmall$species <- NULL
  rownames(annotSpeciesSmall) <- tmp$species
  annotSpeciesSmall <- data.frame(annotSpeciesSmall[rownames(m),])
  rownames(annotSpeciesSmall) <- rownames(m)
  
  myColorAnnot <- list("tree_species"= GetMyColor(annotSample$tree_species),
                       "layer"=GetMyColor(annotSample$layer), 
                       "LifeStyle"=GetMyColor(annotSpeciesSmall$LifeStyle))
  
  quantile_breaks <- function(xs, n = 10) {
    breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
    breaks[!duplicated(breaks)]
  }
  
  mat_breaks <- quantile_breaks(as.matrix(log(m+1)), n = 11)
  
  tryCatch(
    #try to do this
    {
      p <- pheatmap(log(m+1),
                    #color= c(rev(colorRampPalette(brewer.pal(9, "Blues"))(100)), colorRampPalette(brewer.pal(9, "YlOrRd"))(100), "black"),
                    main=nom2,
                    cellwidth=16,
                    cellheight = 4,
                    fontsize_col=8,
                    fontsize_row=4,
                    annotation_col = annotSample,
                    annotation_row = annotSpeciesSmall,
                    annotation_colors = myColorAnnot,
                    color=c("white",brewer.pal(length(mat_breaks) - 3, "OrRd"), "black"),
                    breaks = mat_breaks
                    #color= c("white",(colorRampPalette(brewer.pal(9, "OrRd"))(100)),"black")
                    )
      
      
      #h <- 4 + length(Fungi_ratio$species)/40
      h <- 4 + length(Fungi_ratio$species)/20
      
      pdf(file=nom,width=9, height=h)
      print(p)
      dev.off()
    },
    #if an error occurs, tell me the error
    error=function(e) {
      message('An Error Occurred')
      print(e)
    },
    #if a warning occurs, tell me the warning
    warning=function(w) {
      message('A Warning Occurred')
      print(w)
      return(NA)
    }
  )
  }
# group by ecology /all;  /nb sp; /guild normalised
# heatmap related to fig precedente


#########################

##------------------------------
## --------  CAZymes -----------
##------------------------------
CAZymes_metadata <- read.csv("../0_metadata/metadata_CAZymes_v2.txt", sep="\t")
treatment_metadata$Layer_cluster <- paste(treatment_metadata$Layer,treatment_metadata$cluster_number, sep = " ")

# rm AA7, CE1,CE3 et CE14
Fungi_TPM_CAZymes <- Fungi_TPM %>%
  dplyr::select(-starts_with("NR_"),-starts_with("Myc"),-ID,-ctg_length, -sum_contig_count) %>% 
  filter(!is.na(CAZy_family))

Fungi_TPM_filtered_Annot <- Fungi_TPM_CAZymes %>% separate_rows( "CAZy_family",sep = "-") %>%
                                                  unique() %>%
                                                  filter(CAZy_family %notin%  c("AA7", "CE1","CE3", "CE14")) %>%
                                                  filter(CAZy_family %notin% c("CBM13", "AA6","CBM5","CE4"))


##------------------------------
## -------- New ! 05/07/23 -----------
##------------------------------

Fungi_ratio <- Fungi_TPM_filtered_Annot %>%
  #filter(grepl('ectomycorrhizal',LifeStyle) ) %>%
  #filter(grepl('saprotroph',LifeStyle) ) %>%
  #mutate(CAZy_family = case_when(CAZy_family %in% topCAZymes ~ CAZy_family, TRUE ~ "Others")) %>% 
  #mutate(CAZy_family = case_when(CAZy_family %in% top10CAZymes$CAZy_family ~ CAZy_family, TRUE ~ "Others")) %>% 
  group_by(CAZy_family) %>% summarise_if(is.numeric, sum) %>% 
  pivot_longer(cols=-c(CAZy_family), names_to="samples", values_to="tpm") %>% left_join(CAZymes_metadata) %>%
  group_by(Type) %>% summarise(TPM=sum(tpm)) %>%
  mutate(Type = case_when(!is.na(Type)  ~ Type, TRUE ~ "Others")) %>% ungroup() %>%
  mutate(CAZymeSum=sum(TPM)) %>% group_by(Type) %>% mutate(TPM_perc=(TPM/CAZymeSum))
  #left_join(treatment_metadata) #%>% group_by(CAZy_family, layer_tree) %>% summarise(tpm_median=median(tpm), .groups = "keep") %>%
  #pivot_wider(names_from  = layer_tree, values_from = "tpm_median")
  
Fungi_ratio <- Fungi_ratio %>% ungroup() %>%
arrange(TPM_perc) %>% filter(Type %notin% c("BCW/PCW","MCW"))
Fungi_ratio$Type <- factor(Fungi_ratio$Type, levels=Fungi_ratio$Type) 

toPlot1 <- Fungi_ratio %>% select(Type,TPM_perc ) %>% rename("cat"="Type") %>% mutate(group="CAZyme type")
  

Fungi_ratio <- Fungi_TPM_filtered_Annot %>%
  mutate(LifeStyle=case_when(!grepl('saprotroph',LifeStyle) ~LifeStyle, TRUE ~ "SAP"  )) %>%
  mutate(LifeStyle=case_when(LifeStyle !="ectomycorrhizal" ~LifeStyle, TRUE ~ "EM"  )) %>%
  mutate(LifeStyle = case_when(LifeStyle %in% c("EM","SAP" ) ~ LifeStyle, TRUE ~ "Others")) %>% 
  group_by(LifeStyle) %>% summarise_if(is.numeric, sum) %>% 
  pivot_longer(cols=-c(LifeStyle), names_to="samples", values_to="tpm") %>%# left_join(CAZymes_metadata) %>%
  group_by(LifeStyle) %>% summarise(TPM=sum(tpm)) %>%
  mutate(CAZymeSum=sum(TPM)) %>% group_by(LifeStyle) %>% mutate(TPM_perc=(TPM/CAZymeSum))

Fungi_ratio <- Fungi_ratio %>% ungroup() %>%
  arrange(TPM_perc) 
Fungi_ratio$LifeStyle <- factor(Fungi_ratio$LifeStyle, levels=Fungi_ratio$LifeStyle) 

toPlot2 <- Fungi_ratio %>% select(LifeStyle,TPM_perc ) %>% rename("cat"="LifeStyle") %>% mutate(group="Lifestyle")


Fungi_ratio <- Fungi_TPM_filtered_Annot %>%
  group_by(LifeStyle) %>% summarise_if(is.numeric, sum) %>% 
  pivot_longer(cols=-c(LifeStyle), names_to="samples", values_to="tpm") %>% left_join(treatment_metadata) %>%
  group_by(layer_tree) %>% summarise(TPM=sum(tpm)) %>%
  mutate(CAZymeSum=sum(TPM)) %>% group_by(layer_tree) %>% mutate(TPM_perc=(TPM/CAZymeSum))

Fungi_ratio <- Fungi_ratio %>% ungroup() %>%
  arrange(TPM_perc) 
Fungi_ratio$layer_tree <- factor(Fungi_ratio$layer_tree, levels=Fungi_ratio$layer_tree) 

toPlot3 <- Fungi_ratio %>% select(layer_tree,TPM_perc ) %>% rename("cat"="layer_tree") %>% mutate(group="Location")

toPlot <- rbind(toPlot3,toPlot2,toPlot1)


ggplot(toPlot1, aes(x=cat, y=TPM_perc)) +
  geom_segment( aes(x=cat, xend=cat, y=0, yend=TPM_perc), color="grey") +
  geom_point( color="black", size=3, alpha=0.7) +
  theme_light() +
  coord_flip() +
  labs(y= "CAZymes TPM") +
  theme(
    panel.grid.major.y = element_blank(),
    panel.border = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y=element_blank()
  )+scale_y_continuous(labels = scales::percent) 
# 3*5.83

mycolor <- GetMyColor(toPlot2$cat)
ggplot(toPlot2, aes(x=cat, y=TPM_perc, color=group,)) +
  geom_segment( aes(x=cat, xend=cat, y=0, yend=TPM_perc), color=mycolor) +
  geom_point( size=3, color=mycolor) +
  scale_color_manual(values=mycolor)+
  theme_light() +
  coord_flip() +
  labs(y= "CAZymes TPM") +
  theme(
    panel.grid.major.y = element_blank(),
    panel.border = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y=element_blank()
  )+scale_y_continuous(labels = scales::percent) 

library(patchwork) 
toPlot3 <- Fungi_TPM_filtered_Annot %>%
  group_by(LifeStyle) %>% summarise_if(is.numeric, sum) %>% 
  pivot_longer(cols=-c(LifeStyle), names_to="samples", values_to="tpm") %>% left_join(treatment_metadata) %>%
  group_by(tree_species, Layer,layer_tree) %>% summarise(TPM=sum(tpm)) %>% ungroup()%>%
  mutate(CAZymeSum=sum(TPM)) %>% group_by(layer_tree) %>% mutate(TPM_perc=(TPM/CAZymeSum))

toPlot3 <- toPlot3 %>% ungroup() %>%
  arrange(TPM_perc) 
# Fungi_ratio$layer_tree <- factor(Fungi_ratio$layer_tree, levels=Fungi_ratio$layer_tree) 

mycolor <- GetMyColor(toPlot3$Layer)

toPlot3 %>% 
  ggplot(aes(x=TPM_perc,y=tree_species)) +
  #geom_line(aes(group=tree_species), color="#E7E7E7", linewidth=3.5) + 
  geom_segment( aes(x=0, xend=TPM_perc, y=tree_species, yend=tree_species), color="grey") +
  geom_point(aes(color=Layer), size=3) +
  theme_minimal() +
  theme(legend.position = "bottom",
        #axis.text.y = element_text(color="black"),
        #axis.text.x = element_text(color="#989898"),
        axis.title = element_blank(),
        #panel.grid = element_blank()
  ) +
  scale_color_manual(values=mycolor)+
  scale_x_continuous(labels = scales::percent, limits = c(0,0.23))
  
#######----------

Fungi_ratio <- Fungi_TPM_filtered_Annot %>%
  mutate(LifeStyle=case_when(!grepl('saprotroph',LifeStyle) ~LifeStyle, TRUE ~ "SAP"  )) %>%
  mutate(LifeStyle=case_when(LifeStyle !="ectomycorrhizal" ~LifeStyle, TRUE ~ "EM"  )) %>%
  mutate(LifeStyle = case_when(LifeStyle %in% c("EM","SAP" ) ~ LifeStyle, TRUE ~ "Others")) %>% 
  group_by(CAZy_family, LifeStyle) %>% summarise_if(is.numeric, sum) %>% 
  pivot_longer(cols=-c(CAZy_family,LifeStyle), names_to="samples", values_to="tpm") %>%
  left_join(CAZymes_metadata) %>% ungroup() %>%
  #group_by(Type) %>% 
  mutate(CAZy_Type_sum=sum(tpm)) %>%
  filter(Type=="PCW") %>% ungroup() %>%
  left_join(treatment_metadata) %>% group_by(layer_tree, Layer,tree_species,LifeStyle,CAZy_Type_sum ) %>%
  summarise(TPM=sum(tpm)) %>% mutate(TPM_perc=(TPM/CAZy_Type_sum)*100)

toPlot <- Fungi_ratio %>% ungroup() %>% select(layer_tree, LifeStyle, TPM_perc) %>% pivot_wider(names_from  = LifeStyle, values_from = "TPM_perc") 
colnames(toPlot)
toPlot <- toPlot %>% select(layer_tree,SAP,EM,Others)

m <- toPlot %>% select(-layer_tree)
rownames(m) <- toPlot$layer_tree
breaksList = seq(0,3, by = 0.05)

library(pheatmap)
pheatmap(m,
         #color=myColor,#cutree_rows = 6,
         color= c("white",(colorRampPalette(rev(brewer.pal(9, "RdYlBu")))(length(breaksList)))),
         #color=c("white",brewer.pal(length(mat_breaks)-3, "OrRd"), "black"),
         #breaks = mat_breaks,
         #annotation_col = annotSample,
         #border_color      = "grey80",
         #border_size=.01, 
         #annotation_row = annotCAZymes,
         breaks = breaksList,
         cluster_cols = F,
         cluster_rows = F,
         annotation_colors = myColorAnnot,
         display_numbers = TRUE,
         number_color = "black",
         angle_col=45,
         #cellwidth=16,
         #cellheight = 4,
         fontsize_row = 8,
         fontsize_col = 10#,
         #main="Repartition of FCW CAZymes in the samples (%)"
         #main="PCW"
)
#400 *200

#left_join(treatment_metadata) #%>% group_by(CAZy_family, layer_tree) %>% summarise(tpm_median=median(tpm), .groups = "keep") %>%
#pivot_wider(names_from  = layer_tree, values_from = "tpm_median")


##------------------------------
## -------------------
##------------------------------

#### TPM -- heatmap ------------
Fungi_ratio <- Fungi_TPM_filtered_Annot %>%
  filter(grepl('ectomycorrhizal',LifeStyle) ) %>%
  #filter(grepl('saprotroph',LifeStyle) ) %>%
  #mutate(CAZy_family = case_when(CAZy_family %in% topCAZymes ~ CAZy_family, TRUE ~ "Others")) %>% 
  #mutate(CAZy_family = case_when(CAZy_family %in% top10CAZymes$CAZy_family ~ CAZy_family, TRUE ~ "Others")) %>% 
  group_by(CAZy_family) %>% summarise_if(is.numeric, sum) %>% 
  pivot_longer(cols=-c(CAZy_family), names_to="samples", values_to="tpm") %>%
  left_join(treatment_metadata) %>% group_by(CAZy_family, layer_tree) %>% summarise(tpm_median=median(tpm), .groups = "keep") %>%
  pivot_wider(names_from  = layer_tree, values_from = "tpm_median")

sp_tmp <- Fungi_ratio %>% pivot_longer(cols=-c(CAZy_family), names_to="layer_tree", values_to="tpm_counts") %>% 
  group_by(layer_tree) %>% mutate(sumLayerTree=sum(tpm_counts)) %>% ungroup() %>%
  mutate(observation_perc=(tpm_counts/sumLayerTree)*100  ) %>% 
  #arrange(desc(observation_perc)) %>%
  group_by(layer_tree) %>% 
  dplyr::arrange(layer_tree, observation_perc) %>% 
  dplyr::mutate(cumul_observation_perc=cumsum(observation_perc)) %>%
  dplyr::arrange(layer_tree, desc(observation_perc)) %>% 
  dplyr::mutate(cumul_observation_perc_desc=cumsum(observation_perc)) %>%
  filter(cumul_observation_perc >= 25)

sp_preval <- Fungi_ratio %>% pivot_longer(cols=-c(CAZy_family), names_to="layer_tree", values_to="tpm_counts") %>% 
  mutate(nb=0)%>%
  mutate(nb=case_when(tpm_counts < 1 ~  nb, TRUE ~ 1)) %>%
  group_by(CAZy_family) %>% summarise(preval=sum(nb)) %>%
  filter(preval < 4)


Fungi_ratio <- Fungi_ratio %>% filter(CAZy_family %in% sp_tmp$CAZy_family ) %>% 
  filter(CAZy_family %notin% sp_preval$CAZy_family )
rownames(Fungi_ratio) <- Fungi_ratio$CAZy_family


mycolorsCAZy <- GetMyColorCAZymes(Fungi_ratio$CAZy_family)
#Fungi_ratio$CAZy_family <- factor(Fungi_ratio$CAZy_family, levels=names(mycolorsCAZy)) 

m <- Fungi_ratio %>% ungroup() %>% select(-CAZy_family)
rownames(m) <- Fungi_ratio$CAZy_family

#dat <- data.frame(values = as.numeric(as.matrix(m)))
#ggplot(dat, aes(values)) + geom_density(bw = "SJ")

quantile_breaks <- function(xs, n = 10) {
  breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
  breaks[!duplicated(breaks)]
}

mat_breaks <- quantile_breaks(as.matrix(log(m+1)), n = 11)

library(pheatmap)
#pheatmap(log(m+1))

annotSample <- treatment_metadata %>% select(layer_tree, layer, tree_species) %>% unique()
tmp <- annotSample
annotSample$layer_tree <- NULL
rownames(annotSample) <- tmp$layer_tree

myColorAnnot=list("tree_species"=GetMyColor(annotSample$tree_species),
                  "layer"=GetMyColor(annotSample$layer)#,
                  #"Type"=GetMyColorCAZymes(annotCAZymes$Type),
                  #"Substrate"=GetMyColorCAZymes(annotCAZymes$Substrate)
)

pheatmap(log(m+1),
         #color=myColor,#cutree_rows = 6,
         color= c("white",(colorRampPalette(brewer.pal(9, "OrRd"))(100)),"black"),
         #color=c("white",brewer.pal(length(mat_breaks)-3, "OrRd"), "black"),
         #breaks = mat_breaks,
         annotation_col = annotSample,
         #border_color      = "grey80",
         #border_size=.01, 
         #annotation_row = annotCAZymes,
         annotation_colors = myColorAnnot,
         cellwidth=16,
         cellheight = 4,
         fontsize_row = 4,
         fontsize_col = 8,
         main="CAZymes in Ectomycorrhizal sp."
)


################### boxplots #####################


Fungi_ratio <- Fungi_TPM_filtered_Annot %>%
  left_join(CAZymes_metadata) %>%
  #mutate(Substrate = case_when(!is.na(Substrate) ~ Substrate, TRUE ~ "Unidentified")) %>% 
  mutate(Type = case_when(Type %notin% c("BCW/PCW","MCW", "PCW/FCW") ~ Type, TRUE ~ "Others")) %>% 
  mutate(Type = case_when(!is.na(Type) ~ Type, TRUE ~ "Unidentified")) %>% 
  mutate(LifeStyle=case_when(!grepl('saprotroph',LifeStyle) ~LifeStyle, TRUE ~ "saprotroph"  )) %>%
  mutate(LifeStyle = case_when(LifeStyle %in% c("ectomycorrhizal","saprotroph" ) ~ LifeStyle, TRUE ~ "Others")) %>% 
  #filter(grepl('saprotroph',LifeStyle) ) %>%
  #mutate(CAZy_family = case_when(CAZy_family %in% topCAZymes ~ CAZy_family, TRUE ~ "Others")) %>% 
  #mutate(CAZy_family = case_when(CAZy_family %in% top10CAZymes$CAZy_family ~ CAZy_family, TRUE ~ "Others")) %>% 
  #mutate(TypeLifestyle=paste(Type, LifeStyle, sep="__")) %>%
  group_by(Type,LifeStyle) %>% summarise_if(is.numeric, sum) %>% pivot_longer(cols=-c(Type,LifeStyle), names_to="samples", values_to="tpm_counts") %>%
  left_join(treatment_metadata)  

mycolors <- GetMyColor(Fungi_ratio$LifeStyle)
Fungi_ratio$LifeStyle <- factor(Fungi_ratio$LifeStyle, levels=names(mycolors)) 


ggplot(data=Fungi_ratio, aes(x=Type, y=tpm_counts, colour =LifeStyle)) +
  #geom_col() + 
  geom_boxplot(outlier.size = 1 )+ #outlier.shape = NA)+
  #geom_point(alpha=0.7)+
  #scale_colour_manual(mycolors)+
  theme_bw()+
  #geom_jitter(cex=0.2)+
  facet_wrap(.~ layer*tree_species )+#, scales = "free") + 
  #scale_fill_manual(values=c("#999999",rev(c("#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7","springgreen4","mediumorchid4","cadetblue","lightpink","khaki","turquoise1","black","steelblue")))) +
  theme(axis.text.x=element_text(angle = 55,hjust=1))+ylim(0,4000)


########################################
####--------- Barplot -----------------


  Fungi_ratio <- Fungi_TPM_filtered_Annot %>%
    filter(CAZy_family %notin% c("CBM13", "AA6","CBM5","CE4")) %>% 
    group_by(CAZy_family) %>% summarise_if(is.numeric, sum) %>% pivot_longer(cols=-c(CAZy_family), names_to="samples", values_to="tpm_counts") %>%
    left_join(treatment_metadata) %>%
    left_join(CAZymes_metadata)
  #pivot_wider(names_from = samples, values_from = tpm_counts)

topNb <- 10 

t <- Fungi_ratio %>%
  group_by(CAZy_family) %>% summarise(max=max(tpm_counts), mean=mean(tpm_counts), median=median(tpm_counts)) #%>%
t1 <- t %>% slice_max(.,order_by = max,n=topNb)
t2 <-t  %>% slice_max(.,order_by = median,n=topNb)
t3 <- t %>% slice_max(.,order_by = mean,n=topNb)
#if(t2 > topNb){t2 <- NULL}
topTax <- unique(c(t1$CAZy_family,t2$CAZy_family,t3$CAZy_family))

nbOtherSpecies <- length(t$CAZy_family) - length(topTax)
OtherTag <- paste("Others (n=", nbOtherSpecies,')',sep="")


Fungi_ratio <- Fungi_TPM_filtered_Annot %>%
  mutate(CAZy_family = case_when(CAZy_family %in% topTax ~ CAZy_family, TRUE ~ OtherTag)) %>% 
  group_by(CAZy_family) %>% summarise_if(is.numeric, sum) %>% pivot_longer(cols=-c(CAZy_family), names_to="samples", values_to="tpm_counts") %>%
  left_join(treatment_metadata) %>% 
  left_join(CAZymes_metadata) %>% group_by(layer_tree) %>% mutate(CAZymes_sum=sum(tpm_counts)) %>% ungroup()%>% 
  group_by(layer_tree,CAZy_family,CAZymes_sum) %>% summarise(tpm_counts=sum(tpm_counts))#%>%
  #mutate(tpm_perc=(tpm_counts/CAZymes_sum)*100) 

t <- Fungi_ratio %>% 
  group_by(CAZy_family) %>% summarise(max=max(tpm_counts), mean=mean(tpm_counts), median=median(tpm_counts), min=min(tpm_counts)) #%>%

m <- Fungi_ratio %>% filter(!grepl("Others", CAZy_family)) %>%
  group_by(CAZy_family, layer_tree) %>% summarise(tpm_perc=round(sum(tpm_perc),digits =2)) %>% pivot_wider(names_from=CAZy_family, values_from= tpm_perc)
tmp <- m
m <- m %>% ungroup()%>% select(-layer_tree)
rownames(m) <- tmp$layer_tree
breaksList = seq(0,6, by = 0.1)

library(pheatmap)
pheatmap(m,
         #color=myColor,#cutree_rows = 6,
         color= c("white",(colorRampPalette(rev(brewer.pal(9, "RdYlBu")))(length(breaksList)))),
         #color=c("white",brewer.pal(length(mat_breaks)-3, "OrRd"), "black"),
         #breaks = mat_breaks,
         #annotation_col = annotSample,
         #border_color      = "grey80",
         #border_size=.01, 
         #annotation_row = annotCAZymes,
         breaks = breaksList,
         cluster_cols = F,
       cluster_rows = F,
         annotation_colors = myColorAnnot,
         display_numbers = TRUE,
         number_color = "black",
         angle_col=45,
         #cellwidth=16,
         #cellheight = 4,
         fontsize_row = 8,
         fontsize_col = 10#,
         #main="Repartition of FCW CAZymes in the samples (%)"
         #main="PCW"
)

 #mycolorsCAZy <- GetMyColorCAZymes(Fungi_ratio$CAZy_family)
# Fungi_ratio$CAZy_family <- factor(Fungi_ratio$CAZy_family, levels = names(mycolorsCAZy))


Fungi_ratio2 <- Fungi_ratio %>% filter(!grepl("Others", CAZy_family)) %>% group_by(CAZy_family,layer_tree,layer) %>% summarise(tpm=sum(tpm_counts))

ggplot(data=Fungi_ratio2, aes(x=layer_tree, y=tpm, fill=CAZy_family)) +
  #geom_col() + facet_wrap(.~ layer*tree_species, scales = "free_x") + 
  geom_col() + facet_wrap(.~ layer, scales = "free_x") + 
  #labs(title=nom2) +
  #scale_fill_manual(values=mycolorsCAZy) +
  scale_fill_manual(values=c(brewer.pal(8, "Paired"), brewer.pal(6, "BrBG"), "grey")) +
  #scale_fill_manual(values=c("#999999",rev(c("#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7","springgreen4","mediumorchid4","cadetblue","lightpink","khaki","turquoise1","black","steelblue")))) +
  theme(axis.text.x=element_text(angle = 55,hjust=1))


#####-------------

Fungi_ratio <- Fungi_TPM_filtered_Annot %>%
  filter(CAZy_family %notin% c("CBM13", "AA6","CBM5","CE4")) %>% 
  group_by(CAZy_family) %>% summarise_if(is.numeric, sum) %>% pivot_longer(cols=-c(CAZy_family), names_to="samples", values_to="tpm_counts") %>%
  left_join(treatment_metadata) %>%
  group_by(CAZy_family, layer_tree) %>% summarise(tpm_median=median(tpm_counts), .groups = "keep") #%>%
  #pivot_wider(names_from  = layer_tree, values_from = "tpm_median")
  #left_join(CAZymes_metadata)
#pivot_wider(names_from = samples, values_from = tpm_counts)
  
topNb <- 10 

t <- Fungi_ratio %>%
  group_by(CAZy_family) %>% summarise(max=max(tpm_median), mean=mean(tpm_median), median=median(tpm_median)) #%>%
t1 <- t %>% slice_max(.,order_by = max,n=topNb)
t2 <-t  %>% slice_max(.,order_by = median,n=topNb)
t3 <- t %>% slice_max(.,order_by = mean,n=topNb)
#if(t2 > topNb){t2 <- NULL}
topTax <- unique(c(t1$CAZy_family,t2$CAZy_family,t3$CAZy_family))

nbOtherSpecies <- length(t$CAZy_family) - length(topTax)
OtherTag <- paste("Others (n=", nbOtherSpecies,')',sep="")


Fungi_ratio <- Fungi_TPM_filtered_Annot %>%
  mutate(CAZy_family = case_when(CAZy_family %in% topTax ~ CAZy_family, TRUE ~ OtherTag)) %>% 
  group_by(CAZy_family) %>% summarise_if(is.numeric, sum) %>% pivot_longer(cols=-c(CAZy_family), names_to="samples", values_to="tpm_counts") %>%
  left_join(treatment_metadata) %>%
  group_by(CAZy_family, layer_tree,layer,tree_species) %>% summarise(tpm_median=median(tpm_counts), .groups = "keep")  %>%
  left_join(CAZymes_metadata)


#mycolorsCAZy <- GetMyColorCAZymes(Fungi_ratio$CAZy_family)
# Fungi_ratio$CAZy_family <- factor(Fungi_ratio$CAZy_family, levels = names(mycolorsCAZy))

ggplot(data=Fungi_ratio, aes(x=layer_tree, y=tpm_median, fill=CAZy_family)) +
  geom_col() + facet_wrap(.~ layer, scales = "free_x") + 
  #labs(title=nom2) +
  #scale_fill_manual(values=mycolorsCAZy) +
  scale_fill_manual(values=c(brewer.pal(8, "Paired"), brewer.pal(6, "BrBG"), "grey")) +
  #scale_fill_manual(values=c("#999999",rev(c("#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7","springgreen4","mediumorchid4","cadetblue","lightpink","khaki","turquoise1","black","steelblue")))) +
  theme(axis.text.x=element_text(angle = 55,hjust=1))




##########################
topNb <- 3 
valCheck <- "GH71"
topCAZymes <- c("CBM13", "AA6","CBM5","CE4")
for (valCheck in topCAZymes) {
  
  Fungi_ratio <- Fungi_TPM_filtered_Annot %>% 
    filter(CAZy_family == valCheck)  %>%
    group_by(tax_sciname) %>% summarise_if(is.numeric, sum) %>% pivot_longer(cols=-c(tax_sciname), names_to="samples", values_to="tpm_counts") %>%
    left_join(treatment_metadata)
  #pivot_wider(names_from = samples, values_from = tpm_counts)
  
  t <- Fungi_ratio %>%
    group_by(tax_sciname) %>% summarise(max=max(tpm_counts), mean=mean(tpm_counts), median=median(tpm_counts)) #%>%
  t1 <- t %>% slice_max(.,order_by = max,n=topNb)
  t2 <-t  %>% slice_max(.,order_by = median,n=topNb)
  t3 <- t %>% slice_max(.,order_by = mean,n=topNb)
  if(t2 > topNb){t2 <- NULL}
  topTax <- unique(c(t1$tax_sciname,t2$tax_sciname,t3$tax_sciname))
  
  nbOtherSpecies <- length(t$tax_sciname) - length(topTax)
  OtherTag <- paste("Others (", nbOtherSpecies,' CAZymes)',sep="")
  
  Fungi_ratio <- Fungi_TPM_filtered_Annot %>% 
    filter(CAZy_family == valCheck)  %>%
    mutate(tax_sciname = case_when(tax_sciname %in% topTax ~ tax_sciname, TRUE ~ OtherTag)) %>% 
    group_by(tax_sciname) %>% summarise_if(is.numeric, sum) %>% pivot_longer(cols=-c(tax_sciname), names_to="samples", values_to="tpm_counts") %>%
    left_join(treatment_metadata) 
  
  nom <- paste("CAZyme_", valCheck,".pdf", sep="")
  nom2 <- paste("species composition expressing ", valCheck)
  p <- ggplot(data=Fungi_ratio, aes(x=samples, y=tpm_counts, fill=tax_sciname)) +
    geom_col() + facet_wrap(.~ layer*tree_species, scales = "free_x") + 
    labs(title=nom2) +
    #scale_fill_manual(values=mycolorsCAZy) +
    #scale_fill_manual(values=c("#999999",rev(c("#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7","springgreen4","mediumorchid4","cadetblue","lightpink","khaki","turquoise1","black","steelblue")))) +
    theme(axis.text.x=element_text(angle = 55,hjust=1))
  
  
  pdf(file=nom,width=10, height=7)
  print(p)
  dev.off()
  
}
########## ternary  ----------------
library(ggtern)
Fungi_ratio <- Fungi_TPM_filtered_Annot %>%
  filter(CAZy_family %notin% c("CBM13", "AA6","CBM5","CE4")) %>% 
  group_by(CAZy_family) %>% summarise_if(is.numeric, sum) %>% pivot_longer(cols=-c(CAZy_family), names_to="samples", values_to="tpm_counts") %>%
  left_join(treatment_metadata) 




plotByTree_species<- Fungi_ratio %>% group_by(tree_species,CAZy_family) %>%
  summarise(nb=sum(tpm_counts))%>% ungroup() %>%
  group_by(CAZy_family) %>% mutate(sum_TPM=sum(nb)) %>% ungroup() %>%
  mutate(pcA_per_Name_Tree_species=nb/sum_TPM*100)

pointSize <- plotByTree_species %>% select(CAZy_family, sum_TPM) %>% unique()

plot <- plotByTree_species %>% select(tree_species,pcA_per_Name_Tree_species, CAZy_family ) %>% 
  pivot_wider(names_from = tree_species, values_from = pcA_per_Name_Tree_species, values_fill = 0) %>% 
  left_join(pointSize) %>%
  left_join(CAZymes_metadata) %>%  mutate(Type = case_when(Type %notin% c("BCW/PCW","MCW", "PCW/FCW") ~ Type, TRUE ~ "Others")) 

rownames(plot) <- plot$CAZy_family

mycolorsCAZy <- GetMyColorCAZymes(plot$Type)
plot$Type <- factor(plot$Type, levels=names(mycolorsCAZy)) 

ggtern(data=plot, aes( Quercus,Abies, Picea))+
  geom_point( shape=21, mapping=aes(size=sum_TPM, fill=Type )) + 
  scale_fill_manual(values=mycolorsCAZy)+
  theme_rgbw()

###############
#mycorrhiz
#saprotroph

Fungi_ratio <-  Fungi_TPM_filtered_Annot %>%
  filter(CAZy_family %notin% c("CBM13", "AA6","CBM5","CE4")) %>% 
  filter(grepl('saprotroph',LifeStyle) ) %>%
  group_by(CAZy_family) %>% summarise_if(is.numeric, sum) %>% pivot_longer(cols=-c(CAZy_family), names_to="samples", values_to="tpm_counts") %>%
  left_join(treatment_metadata)

plotByTree_species<- Fungi_ratio %>% group_by(tree_species,CAZy_family) %>%
  summarise(nb=sum(tpm_counts))%>% ungroup() %>%
  group_by(CAZy_family) %>% mutate(sum_TPM=sum(nb)) %>% ungroup() %>%
  mutate(pcA_per_Name_Tree_species=nb/sum_TPM*100)

pointSize <- plotByTree_species %>% select(CAZy_family, sum_TPM) %>% unique()

plot <- plotByTree_species %>% select(tree_species,pcA_per_Name_Tree_species, CAZy_family ) %>% 
  pivot_wider(names_from = tree_species, values_from = pcA_per_Name_Tree_species, values_fill = 0) %>% 
  left_join(pointSize) %>%
  left_join(CAZymes_metadata) %>%  mutate(Type = case_when(Type %notin% c("BCW/PCW","MCW", "PCW/FCW") ~ Type, TRUE ~ "Others")) 

rownames(plot) <- plot$CAZy_family

mycolorsCAZy <- GetMyColorCAZymes(plot$Type)
plot$Type <- factor(plot$Type, levels=names(mycolorsCAZy)) 

SAP <- ggtern(data=plot, aes( Quercus,Abies, Picea))+
  geom_point( shape=21, mapping=aes(size=sum_TPM, fill=Type )) + 
  scale_size_continuous(range = c(1,10))+
  scale_fill_manual(values=mycolorsCAZy)+
  theme_rgbw()+
  guides(fill = guide_legend(override.aes = list(size = 5)))
SAP

##########
Fungi_ratio <-  Fungi_TPM_filtered_Annot %>%
  filter(CAZy_family %notin% c("CBM13", "AA6","CBM5","CE4")) %>% 
  filter(grepl('mycorrhiz',LifeStyle) ) %>%
  group_by(CAZy_family) %>% summarise_if(is.numeric, sum) %>% pivot_longer(cols=-c(CAZy_family), names_to="samples", values_to="tpm_counts") %>%
  left_join(treatment_metadata)

plotByTree_species<- Fungi_ratio %>% group_by(tree_species,CAZy_family) %>%
  summarise(nb=sum(tpm_counts))%>% ungroup() %>%
  group_by(CAZy_family) %>% mutate(sum_TPM=sum(nb)) %>% ungroup() %>%
  mutate(pcA_per_Name_Tree_species=nb/sum_TPM*100)

pointSize <- plotByTree_species %>% select(CAZy_family, sum_TPM) %>% unique()

plot <- plotByTree_species %>% select(tree_species,pcA_per_Name_Tree_species, CAZy_family ) %>% 
  pivot_wider(names_from = tree_species, values_from = pcA_per_Name_Tree_species, values_fill = 0) %>% 
  left_join(pointSize) %>%
  left_join(CAZymes_metadata)%>%  mutate(Type = case_when(Type %notin% c("BCW/PCW","MCW", "PCW/FCW") ~ Type, TRUE ~ "Others")) 

rownames(plot) <- plot$CAZy_family

mycolorsCAZy <- GetMyColorCAZymes(plot$Type)
plot$Type <- factor(plot$Type, levels=names(mycolorsCAZy)) 

Myc <- ggtern(data=plot, aes( Quercus,Abies, Picea))+
  geom_point( shape=21, mapping=aes(size=sum_TPM, fill=Type )) + 
  scale_fill_manual(values=mycolorsCAZy)+
  scale_size_continuous(range = c(0,6))+
  theme_rgbw()+
  guides(fill = guide_legend(override.aes = list(size = 5)))
Myc




########################################
####--------- Heatmap KOG presence/absence  -----------------

WhereToTest <- "KOG_class"
Fungi_TPM_filtered_Annot <-Split_Annot_In_Multiple_Rows(Fungi_TPM, WhereToTest)
Fungi_TPM_filtered_Annot <- Fungi_TPM_filtered_Annot %>% filter(KOG_group %in% c("CELLULAR PROCESSES AND SIGNALING", "INFORMATION STORAGE AND PROCESSING", "METABOLISM"))
Fungi_TPM_filtered_Annot$KOG_class[Fungi_TPM_filtered_Annot$KOG_class == "Posttranslational modification, protein turnover, chaperones"] <- "Post-translational modification, protein turnover and chaperones"

# 
rownames(Fungi_TPM_filtered_Annot) <-paste( Fungi_TPM_filtered_Annot$species, rownames(Fungi_TPM_filtered_Annot), sep="_") ## just to check if it's unique
Fungi_TPM_filtered_Annot$ID <- paste( Fungi_TPM_filtered_Annot$species, rownames(Fungi_TPM_filtered_Annot), sep="_") 


Fungi_ratio <- get.Fungi_ratio.2groupments(Fungi_TPM_filtered_Annot, WhereToTest="KOG_class",WhereToTest2="species",SelectMetadataValue="all", columnMetadataValue="all")
t <- Fungi_ratio %>% select(KOG_class, species,tpm_counts) %>%
  group_by(KOG_class, species) %>% summarise(tpm=sum(tpm_counts))

m <- t %>% select(KOG_class, species,tpm) %>% mutate(tpm2= case_when(tpm == 0 ~ tpm, TRUE ~ 1)) %>% select(-c(tpm))%>%
  pivot_wider(names_from = species, values_from = tpm2, values_fill = 0 )
tmp <- m
m$KOG_class <- NULL
rownames(m) <- tmp$KOG_class

eco <- Fungi_TPM_filtered_Annot %>% select(LifeStyle, species) %>% unique()

# rownames(eco) <- eco$species
# eco <- eco[colnames(m),]
# tmp <- eco
# eco$species <- NULL
# rownames(eco) <- tmp$species


library(pheatmap)
pheatmap(m,         show_colnames=F,
         #annotation_col = eco,
         color = colorRampPalette(brewer.pal(9, "Blues"))(100)
)

sap <- eco %>%   filter(grepl('saprotroph',LifeStyle) )
ecm <- eco %>% filter(LifeStyle == "ectomycorrhizal")
m2 <- m %>% select(sap$species)
rownames(m2) <- rownames(m)
m3 <- m %>% select(ecm$species)
rownames(m3) <- rownames(m)

pheatmap(m2,         show_colnames=F,
         #annotation_col = eco,
         color = colorRampPalette(brewer.pal(9, "Blues"))(100)
)

pheatmap(m3,         #show_colnames=F,
         #annotation_col = eco,
         fontsize_col=3,
         
         color = colorRampPalette(brewer.pal(9, "Blues"))(100)
)

####### CAZyme version #########
#-------------------------------
Fungi_TPM_CAZymes <- Fungi_TPM %>%
  dplyr::select(-starts_with("NR_"),-starts_with("Myc"),-ID,-ctg_length, -sum_contig_count ) %>% 
  filter(!is.na(CAZy_family))

Fungi_TPM_filtered_Annot <- Fungi_TPM_CAZymes %>% separate_rows( "CAZy_family",sep = "-") %>%
  unique() %>%
  filter(CAZy_family %notin%  c("AA7", "CE1","CE3", "CE14","CBM13", "AA6","CBM5","CE4")) 
rownames(Fungi_TPM_filtered_Annot) <-paste( Fungi_TPM_filtered_Annot$species, rownames(Fungi_TPM_filtered_Annot), sep="_") ## just to check if it's unique
Fungi_TPM_filtered_Annot$ID <- paste( Fungi_TPM_filtered_Annot$species, rownames(Fungi_TPM_filtered_Annot), sep="_") 


Fungi_ratio <- get.Fungi_ratio.2groupments(Fungi_TPM_filtered_Annot, WhereToTest="CAZy_family",WhereToTest2="species",SelectMetadataValue="all", columnMetadataValue="all")
t <- Fungi_ratio %>% select(CAZy_family, species,tpm_counts) %>%
  group_by(CAZy_family, species) %>% summarise(tpm=sum(tpm_counts))

m <- t %>% select(CAZy_family, species,tpm) %>% mutate(tpm2= case_when(tpm == 0 ~ tpm, TRUE ~ 1)) %>% select(-c(tpm))%>%
  pivot_wider(names_from = species, values_from = tpm2, values_fill = 0 )
tmp <- m
m$CAZy_family <- NULL
rownames(m) <- tmp$CAZy_family

eco <- Fungi_TPM_filtered_Annot %>% select(LifeStyle, species) %>% unique()

library(pheatmap)
pheatmap(m,         show_colnames=F,
         #annotation_col = eco,
         color = colorRampPalette(brewer.pal(9, "Blues"))(100)
)

sap <- eco %>%   filter(grepl('saprotroph',LifeStyle) )
ecm <- eco %>% filter(LifeStyle == "ectomycorrhizal")
m2 <- m %>% select(sap$species)
rownames(m2) <- rownames(m)
m3 <- m %>% select(ecm$species)
rownames(m3) <- rownames(m)

pheatmap(m2,         show_colnames=F,
         fontsize_row =3,
         
         #annotation_col = eco,
         color = colorRampPalette(brewer.pal(9, "Blues"))(100)
)

pheatmap(m3,         #show_colnames=F,
         #annotation_col = eco,
         fontsize_row =3,
         fontsize_col =4,
         color = colorRampPalette(brewer.pal(9, "Blues"))(100)
)

###########

df <- Fungi_TPM_filtered_Annot %>%  group_by(CAZy_family, species) %>% 
  summarise_if(is.numeric, sum) %>% 
  pivot_longer(cols=-c(CAZy_family,species), names_to="samples", values_to="tpm_counts") %>% 
  filter(tpm_counts > 20) %>% select(CAZy_family,species) %>% unique() %>%
  group_by(CAZy_family) %>% summarise(nb=n())

df2 <-  df %>%
  group_by(group = cut(nb, breaks = c(0,1,5,10,25,50,100))) %>%
  summarise(n = n())

ggplot(data=df2, aes(x=group, y=n)) +
  geom_bar(stat="identity")+
  theme_bw()+
  ylab("Nb of species") + 
  xlab("Nb of CAZymes families")+
  ggtitle("Number of assigned species per CAZymes family detected")+
  geom_text(aes(label=n), vjust=-0.5)

#####

df <- Fungi_TPM_filtered_Annot %>%   filter(grepl('saprotroph',LifeStyle) ) %>%
  group_by(CAZy_family, species) %>% 
  summarise_if(is.numeric, sum) %>% 
  pivot_longer(cols=-c(CAZy_family,species), names_to="samples", values_to="tpm_counts") %>% 
  filter(tpm_counts > 20) %>% select(CAZy_family,species) %>% unique() %>%
  group_by(CAZy_family) %>% summarise(nb=n())

df2 <-  df %>%
  group_by(group = cut(nb, breaks = c(0,1,5,10,25,50,100))) %>%
  summarise(n = n())

ggplot(data=df2, aes(x=group, y=n)) +
  geom_bar(stat="identity")+
  theme_bw()+
  ylab("Nb of species") + 
  xlab("Nb of CAZymes families")+
  ggtitle("Number of assigned SAP species per CAZymes family detected")+
  geom_text(aes(label=n), vjust=-0.5)

#####

df <- Fungi_TPM_filtered_Annot %>%   filter(grepl('ectomycorr',LifeStyle) ) %>%
  group_by(CAZy_family, species) %>% 
  summarise_if(is.numeric, sum) %>% 
  pivot_longer(cols=-c(CAZy_family,species), names_to="samples", values_to="tpm_counts") %>% 
  filter(tpm_counts > 20) %>% select(CAZy_family,species) %>% unique() %>%
  group_by(CAZy_family) %>% summarise(nb=n())

df2 <-  df %>%
  group_by(group = cut(nb, breaks = c(0,1,5,10,25,50,100))) %>%
  summarise(n = n())

ggplot(data=df2, aes(x=group, y=n)) +
  geom_bar(stat="identity")+
  theme_bw()+
  ylab("Nb of species") + 
  xlab("Nb of CAZymes families")+
  ggtitle("Number of assigned EM species per CAZymes family detected")+
  geom_text(aes(label=n), vjust=-0.5)



###### Proteases #########

Fungi_TPM_MEROPS <- Fungi_TPM %>%
  dplyr::select(-starts_with("NR_"),-starts_with("Myc"),-ID,-ctg_length, -sum_contig_count ) %>% 
  filter(MEROPS_Hit != "*")

Fungi_TPM_MEROPS$MEROPS_ID <- sapply(strsplit(Fungi_TPM_MEROPS$MEROPS_description, "#"), "[",2)


Fungi_ratio <- Fungi_TPM_MEROPS %>%
  #filter(grepl('ectomycorrhizal',LifeStyle) ) %>%
  filter(grepl('saprotroph',LifeStyle) ) %>%
  #mutate(CAZy_family = case_when(CAZy_family %in% topCAZymes ~ CAZy_family, TRUE ~ "Others")) %>% 
  #mutate(CAZy_family = case_when(CAZy_family %in% top10CAZymes$CAZy_family ~ CAZy_family, TRUE ~ "Others")) %>% 
  group_by(MEROPS_ID) %>% summarise_if(is.numeric, sum) %>% 
  pivot_longer(cols=-c(MEROPS_ID), names_to="samples", values_to="tpm") %>%
  left_join(treatment_metadata) %>% group_by(MEROPS_ID, layer_tree) %>% summarise(tpm_median=median(tpm), .groups = "keep") %>%
  pivot_wider(names_from  = layer_tree, values_from = "tpm_median")

sp_tmp <- Fungi_ratio %>% pivot_longer(cols=-c(MEROPS_ID), names_to="layer_tree", values_to="tpm_counts") %>% 
  group_by(layer_tree) %>% mutate(sumLayerTree=sum(tpm_counts)) %>% ungroup() %>%
  mutate(observation_perc=(tpm_counts/sumLayerTree)*100  ) %>% 
  #arrange(desc(observation_perc)) %>%
  group_by(layer_tree) %>% 
  dplyr::arrange(layer_tree, observation_perc) %>% 
  dplyr::mutate(cumul_observation_perc=cumsum(observation_perc)) %>%
  dplyr::arrange(layer_tree, desc(observation_perc)) %>% 
  dplyr::mutate(cumul_observation_perc_desc=cumsum(observation_perc)) %>%
  filter(cumul_observation_perc >= 50)

sp_preval <- Fungi_ratio %>% pivot_longer(cols=-c(MEROPS_ID), names_to="layer_tree", values_to="tpm_counts") %>% 
  mutate(nb=0)%>%
  mutate(nb=case_when(tpm_counts < 1 ~  nb, TRUE ~ 1)) %>%
  group_by(MEROPS_ID) %>% summarise(preval=sum(nb)) %>%
  filter(preval < 4)


Fungi_ratio <- Fungi_ratio %>% filter(MEROPS_ID %in% sp_tmp$MEROPS_ID ) %>% 
  filter(MEROPS_ID %notin% sp_preval$MEROPS_ID )
rownames(Fungi_ratio) <- Fungi_ratio$MEROPS_ID


#mycolorsCAZy <- GetMyColorCAZymes(Fungi_ratio$CAZy_family)
#Fungi_ratio$CAZy_family <- factor(Fungi_ratio$CAZy_family, levels=names(mycolorsCAZy)) 

m <- Fungi_ratio %>% ungroup() %>% select(-MEROPS_ID)
rownames(m) <- Fungi_ratio$MEROPS_ID

#dat <- data.frame(values = as.numeric(as.matrix(m)))
#ggplot(dat, aes(values)) + geom_density(bw = "SJ")

quantile_breaks <- function(xs, n = 10) {
  breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
  breaks[!duplicated(breaks)]
}

mat_breaks <- quantile_breaks(as.matrix(log(m+1)), n = 11)

library(pheatmap)
#pheatmap(log(m+1))

annotSample <- treatment_metadata %>% select(layer_tree, layer, tree_species) %>% unique()
tmp <- annotSample
annotSample$layer_tree <- NULL
rownames(annotSample) <- tmp$layer_tree

myColorAnnot=list("tree_species"=GetMyColor(annotSample$tree_species),
                  "layer"=GetMyColor(annotSample$layer)#,
                  #"Type"=GetMyColorCAZymes(annotCAZymes$Type),
                  #"Substrate"=GetMyColorCAZymes(annotCAZymes$Substrate)
)

pheatmap(log(m+1),
         #color=myColor,#cutree_rows = 6,
         color= c("white",(colorRampPalette(brewer.pal(9, "OrRd"))(100)),"black"),
         #color=c("white",brewer.pal(length(mat_breaks)-3, "OrRd"), "black"),
         #breaks = mat_breaks,
         annotation_col = annotSample,
         #border_color      = "grey80",
         #border_size=.01, 
         #annotation_row = annotCAZymes,
         annotation_colors = myColorAnnot,
         cellwidth=16,
         cellheight = 4,
         fontsize_row = 4,
         fontsize_col = 8,
         main="MEROPS in Saprotrophs sp."
)

###########

Fungi_ratio <- Fungi_TPM_MEROPS %>%
  mutate(LifeStyle=case_when(!grepl('saprotroph',LifeStyle) ~LifeStyle, TRUE ~ "SAP"  )) %>%
  mutate(LifeStyle=case_when(LifeStyle !="ectomycorrhizal" ~LifeStyle, TRUE ~ "EM"  )) %>%
  mutate(LifeStyle = case_when(LifeStyle %in% c("EM","SAP" ) ~ LifeStyle, TRUE ~ "Others")) %>% 
  group_by(MEROPS_ID, LifeStyle) %>% summarise_if(is.numeric, sum) %>% 
  pivot_longer(cols=-c(MEROPS_ID,LifeStyle), names_to="samples", values_to="tpm") %>%
  ungroup() %>%
  #group_by(Type) %>% 
  mutate(MEROPS_sum=sum(tpm)) %>%
   ungroup() %>%
  left_join(treatment_metadata) %>% group_by(layer_tree, Layer,tree_species,LifeStyle,MEROPS_sum ) %>%
  summarise(TPM=sum(tpm)) %>% mutate(TPM_perc=(TPM/MEROPS_sum)*100)

toPlot <- Fungi_ratio %>% ungroup() %>% select(layer_tree, LifeStyle, TPM_perc) %>% pivot_wider(names_from  = LifeStyle, values_from = "TPM_perc") 
colnames(toPlot)
toPlot <- toPlot %>% select(layer_tree,SAP,EM,Others)

m <- toPlot %>% select(-layer_tree)
rownames(m) <- toPlot$layer_tree
breaksList = seq(0,10, by = 0.05)

library(pheatmap)
pheatmap(m,
         #color=myColor,#cutree_rows = 6,
         color= c("white",(colorRampPalette(rev(brewer.pal(9, "RdYlBu")))(length(breaksList)))),
         #color=c("white",brewer.pal(length(mat_breaks)-3, "OrRd"), "black"),
         #breaks = mat_breaks,
         #annotation_col = annotSample,
         #border_color      = "grey80",
         #border_size=.01, 
         #annotation_row = annotCAZymes,
         breaks = breaksList,
         cluster_cols = F,
         cluster_rows = F,
         annotation_colors = myColorAnnot,
         display_numbers = TRUE,
         number_color = "black",
         angle_col=45,
         #cellwidth=16,
         #cellheight = 4,
         fontsize_row = 8,
         fontsize_col = 10,
         main="Repartition of proteases 
         among samples (%)"
         #main="PCW"
)
#400 *200


######

Fungi_ratio <- Fungi_TPM_MEROPS %>%
  group_by(MEROPS_ID) %>% summarise_if(is.numeric, sum) %>% pivot_longer(cols=-c(MEROPS_ID), names_to="samples", values_to="tpm_counts") %>%
  left_join(treatment_metadata) %>%
  group_by(MEROPS_ID, layer_tree) %>% summarise(tpm_median=median(tpm_counts), .groups = "keep") #%>%
#pivot_wider(names_from  = layer_tree, values_from = "tpm_median")
#left_join(CAZymes_metadata)
#pivot_wider(names_from = samples, values_from = tpm_counts)

topNb <- 10 

t <- Fungi_ratio %>%
  group_by(MEROPS_ID) %>% summarise(max=max(tpm_median), mean=mean(tpm_median), median=median(tpm_median)) #%>%
t1 <- t %>% slice_max(.,order_by = max,n=topNb)
t2 <-t  %>% slice_max(.,order_by = median,n=topNb)
t3 <- t %>% slice_max(.,order_by = mean,n=topNb)
#if(t2 > topNb){t2 <- NULL}
topTax <- unique(c(t1$MEROPS_ID,t2$MEROPS_ID,t3$MEROPS_ID))

nbOtherSpecies <- length(t$MEROPS_ID) - length(topTax)
OtherTag <- paste("Others (n=", nbOtherSpecies,')',sep="")


Fungi_ratio <- Fungi_TPM_MEROPS %>%
  mutate(MEROPS_ID = case_when(MEROPS_ID %in% topTax ~ MEROPS_ID, TRUE ~ OtherTag)) %>% 
  group_by(MEROPS_ID) %>% summarise_if(is.numeric, sum) %>% pivot_longer(cols=-c(MEROPS_ID), names_to="samples", values_to="tpm_counts") %>%
  left_join(treatment_metadata) %>%
  group_by(MEROPS_ID, layer_tree,layer,tree_species) %>% summarise(tpm=sum(tpm_counts), .groups = "keep")  %>%
  filter(!grepl("Others",MEROPS_ID))


#mycolorsCAZy <- GetMyColorCAZymes(Fungi_ratio$CAZy_family)
Fungi_ratio$MEROPS_ID <- factor(Fungi_ratio$MEROPS_ID, levels = c(topTax,OtherTag) )

ggplot(data=Fungi_ratio, aes(x=layer_tree, y=tpm, fill=MEROPS_ID)) +
  geom_col() + facet_wrap(.~ layer, scales = "free_x") + 
  #labs(title=nom2) +
  #scale_fill_manual(values=mycolorsCAZy) +
  scale_fill_manual(values=c(brewer.pal(8, "Paired"), brewer.pal(6, "BrBG"), "grey")) +
  #scale_fill_manual(values=c("#999999",rev(c("#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7","springgreen4","mediumorchid4","cadetblue","lightpink","khaki","turquoise1","black","steelblue")))) +
  theme(axis.text.x=element_text(angle = 55,hjust=1))




############# TEST with Turkey ###########
WhereToTest <- "KOG_class"
Fungi_TPM_filtered_Annot <-Split_Annot_In_Multiple_Rows(Fungi_TPM, WhereToTest)
Fungi_TPM_filtered_Annot <- Fungi_TPM_filtered_Annot %>% filter(KOG_group %in% c("CELLULAR PROCESSES AND SIGNALING", "INFORMATION STORAGE AND PROCESSING", "METABOLISM"))
Fungi_TPM_filtered_Annot$KOG_class[Fungi_TPM_filtered_Annot$KOG_class == "Posttranslational modification, protein turnover, chaperones"] <- "Post-translational modification, protein turnover and chaperones"

# 
rownames(Fungi_TPM_filtered_Annot) <-paste( Fungi_TPM_filtered_Annot$species, rownames(Fungi_TPM_filtered_Annot), sep="_") ## just to check if it's unique
Fungi_TPM_filtered_Annot$ID <- paste( Fungi_TPM_filtered_Annot$species, rownames(Fungi_TPM_filtered_Annot), sep="_") 


WhereToTest2 <- "LifeStyle"

Fungi_ratio <- get.Fungi_ratio.2groupments(Fungi_TPM_filtered_Annot, WhereToTest, WhereToTest2="LifeStyle",SelectMetadataValue="all", columnMetadataValue="all")

Fungi_ratio <-Fungi_ratio %>%
  mutate(LifeStyle=case_when(!grepl('saprotroph',LifeStyle) ~LifeStyle, TRUE ~ "saprotroph"  )) %>%
  mutate(LifeStyle = case_when(LifeStyle %in% c("ectomycorrhizal","saprotroph" ) ~ LifeStyle, TRUE ~ "Others")) %>%
  group_by(KOG_class, LifeStyle, samples) %>% summarise(tpm=sum(tpm_counts))
colnames(Fungi_ratio)[colnames(Fungi_ratio) == "tpm"] <- "tpm_counts"

dfstatTPM <- get_dfstatTPM_LifeStyle(Fungi_ratio, WhereToTest)
dfstatPERC <- get_dfstatPERC_LifeStyle(Fungi_ratio, WhereToTest)

###### plot ###
dfTPM <- Fungi_ratio %>%   left_join(treatment_metadata) %>%
  group_by(KOG_class, LifeStyle,layer,tree_species, layer_tree) %>% summarise(mean_tpm=mean(tpm_counts)) %>% join(KOG_meta)
#dfTPM <- dfTPM %>% left_join(sumlayertree) %>% mutate(perc=(tpm/sommeLayerTree)*100)

mycolors <- GetMyColor(dfTPM$LifeStyle)
dfTPM$LifeStyle <- factor(dfTPM$LifeStyle, levels=names(mycolors)) 

ggplot(data=dfTPM, aes(x=Letter_code, y=mean_tpm, fill=LifeStyle , alpha=`KOG class` )) +
  geom_col() + scale_fill_manual(values=mycolors)+
  scale_alpha_manual(values=rep(1,times=23)) +
  theme_bw()+
  #theme(legend.position="bottom")+
  #ylim(0, 2000000)+
  facet_wrap(. ~ layer*tree_species, scales = "free_x") #+ theme(axis.text.x=element_text(angle = 55,hjust=1)) 
#size 7*18
#### transform to heatmap ?

###################
mycolors <- GetMyColor(Fungi_TPM_filtered_Annot$LifeStyle)
Fungi_TPM_filtered_Annot$LifeStyle <- factor(Fungi_TPM_filtered_Annot$LifeStyle, levels=names(mycolors)) 

GroupToTest <- unique(Fungi_TPM_filtered_Annot$KOG_class) 
catToTest <- "Translation, ribosomal structure and biogenesis"  
#for(catToTest in GroupToTest){
 

 Fungi_ratio <- Fungi_TPM_filtered_Annot %>%
    dplyr::select(-starts_with("NR_"),-starts_with("Myc"),-ID,-ctg_length, -sum_contig_count) %>%
    filter(KOG_class == catToTest) %>%
    group_by(LifeStyle) %>% summarise_if(is.numeric, sum) %>% pivot_longer(cols=-c(LifeStyle), names_to="samples", values_to="tpm_counts") %>%
    left_join(treatment_metadata) %>%
    group_by(samples, layer, tree_species,layer_tree) %>% summarise(tpm=sum(tpm_counts))
    #group_by(LifeStyle, layer, layer_tree) %>% summarise(tpm_mean=mean(tpm_counts))
# anova
 
  library(multcompView)
  model <- aov(tpm~layer_tree, data=Fungi_ratio)
  summary(model)
  tukey.res <-TukeyHSD(model, conf.level=.95)
  
  generate_label_df <- function(tukey.res, variable){
    # Extract labels and factor levels from Tukey post-hoc 
    Tukey.levels <- tukey.res[[variable]][,4]
    Tukey.labels <- data.frame(multcompLetters(Tukey.levels)['Letters'])
    
    #I need to put the labels in the same order as in the boxplot :
    Tukey.labels$treatment=rownames(Tukey.labels)
    Tukey.labels=Tukey.labels[order(Tukey.labels$treatment) , ]
    return(Tukey.labels)
  }
  
  tmp <-Fungi_ratio %>% ungroup()%>% select(layer_tree,layer, tree_species ) %>% unique()
  # Apply the function on my dataset
  LABELS <- generate_label_df(tukey.res , "layer_tree")
  LABELS<-  LABELS %>% rename("layer_tree"=treatment) %>% left_join(tmp)
  
  
  
  mycolors <- GetMyColor(Fungi_ratio$tree_species)
  
  
  ggplot(data=Fungi_ratio, aes(x=layer_tree, y=tpm, fill=tree_species)) +
    geom_boxplot(outlier.size = 0.5) + #facet_wrap(. ~ layer*tree_species, scales = "free_x") +
    ggtitle(catToTest)+
    scale_fill_manual(values=mycolors)+
    scale_y_continuous(labels = function(x) format(x, scientific = FALSE))+
    geom_text(
      data = LABELS,
      aes(
        x = layer_tree, 
        y = max(Fungi_ratio$tpm)+10000, 
        label = Letters
      ))  +  facet_wrap(. ~ layer, scales = "free_x") +
    #c("#999999",rev(c("#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7","springgreen4","mediumorchid4","cadetblue","lightpink","khaki","turquoise1")))) +
    theme_bw()+
    theme(axis.text.x=element_text(angle = 55,hjust=1), axis.title.x=element_blank())
  
  