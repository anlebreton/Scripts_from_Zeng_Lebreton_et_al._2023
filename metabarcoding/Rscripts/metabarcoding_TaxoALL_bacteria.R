# ## Import functions
library("plyr");library("dplyr");library("ggplot2");library("tidyr");library("tidyverse")
library(pheatmap)
library(RColorBrewer)


samplesMetaData <- read.csv("../0_Metadata/Yulong_sampling_bacteria_metadata.txt", sep='\t')
samplesMetaData <- samplesMetaData %>% mutate(samples= gsub("-", ".", X)) 
data <- read.csv("../1_data/Yulong_bacteria_normalised_abundance_addCols.tsv", sep='\t')
#data[ data == "unidentified" ] <- "Unidentified" 
#data[ data == "Multi-affiliation" ] <- "Unidentified" 

amplicon="bacteria"

mycolorsFullList <- read.csv("../..//mycolorFile.txt",sep="\t")
mycolorVector <- mycolorsFullList$color
names(mycolorVector) <- mycolorsFullList$Name


########## Tree ########

############################### Phylum -grouped 
top10Taxo<- data %>%
  #dplyr::select(-starts_with("NR_"),-starts_with("Myc"),-ID,-ctg_length, -sum_contig_count) %>% 
  dplyr::select(starts_with(amplicon), Phylum) %>%
  group_by(Phylum) %>% summarise_if(is.numeric, sum) %>%
  pivot_longer(cols=-c(Phylum), names_to="samples", values_to="counts")%>%
  left_join(., samplesMetaData, by="samples") %>% 
  group_by(Phylum, Cluster,Tree_species) %>%
  summarise("reads_count"=sum(counts)) %>% 
  group_by(Phylum) %>%
  summarise( max=max(reads_count), mean=mean(reads_count), median=median(reads_count)) %>%
  slice_max(.,order_by = max,n=12)




df <- data %>%
  #dplyr::select(-starts_with("NR_"),-starts_with("Myc"),-ID,-ctg_length, -sum_contig_count) %>% 
  dplyr::select(starts_with(amplicon), Phylum) %>%
  mutate(Phylum = case_when(Phylum %in% top10Taxo$Phylum ~ Phylum, TRUE ~ "Others")) %>% 
  group_by(Phylum) %>% summarise_if(is.numeric, sum) %>%
  pivot_longer(cols=-c(Phylum), names_to="samples", values_to="counts")%>%
  left_join(., samplesMetaData, by="samples") %>% 
  #group_by(Order, Cluster,Season, Layer, Tree_species) %>%
  group_by(Phylum, Cluster,Tree_species) %>%
  summarise("reads_count"=sum(counts)) %>% 
  ungroup() %>%
  group_by(Cluster) %>% mutate("reads_percent"= (reads_count/sum(reads_count))*100 )

# samplesMetaData <- samplesMetaData %>% arrange(Tree_number, Season, Layer, samples)
# t <- samplesMetaData$samples 
# df$samples <- factor(df$samples, levels= t )

mycolors <- mycolorsFullList %>% filter(Name %in% df$Phylum)
mycolors <- mycolors %>% arrange(myLevel)
mycolor <- mycolors$color
names(mycolor) <- mycolors$Name
df$Phylum <- factor(df$Phylum, levels=mycolors$Name)


pp1 <- ggplot(data=df,aes(y=reads_percent,x=Cluster, fill=Phylum))+
  geom_bar(stat="identity")+
  theme_bw()+
  scale_fill_manual(values=mycolor)+
  theme(axis.text.x=element_text(angle = 55,hjust=1#,size=2
                                ), axis.title.x = element_blank()
        )+

  #facet_wrap(~ Season*Tree_species, scales = "free_x") 
  facet_wrap(~ Tree_species, scales = "free_x") + 
  guides(colour = guide_legend(title.position = "right"))#+
  #theme(legend.position="bottom")

pp1

############################### Class - grouped 
top10Taxo<- data %>%
  #dplyr::select(-starts_with("NR_"),-starts_with("Myc"),-ID,-ctg_length, -sum_contig_count) %>% 
  dplyr::select(starts_with(amplicon), Class) %>%
  group_by(Class) %>% summarise_if(is.numeric, sum) %>%
  pivot_longer(cols=-c(Class), names_to="samples", values_to="counts")%>%
  left_join(., samplesMetaData, by="samples") %>% 
  group_by(Class, Cluster,Tree_species) %>%
  summarise("reads_count"=sum(counts)) %>% 
  group_by(Class) %>%
  summarise( max=max(reads_count), mean=mean(reads_count), median=median(reads_count)) %>%
  slice_max(.,order_by = max,n=12)



df <- data %>%
  #dplyr::select(-starts_with("NR_"),-starts_with("Myc"),-ID,-ctg_length, -sum_contig_count) %>% 
  dplyr::select(starts_with(amplicon), Class) %>%
  mutate(Class = case_when(Class %in% top10Taxo$Class ~ Class, TRUE ~ "Others")) %>% 
  group_by(Class) %>% summarise_if(is.numeric, sum) %>%
  pivot_longer(cols=-c(Class), names_to="samples", values_to="counts")%>%
  left_join(., samplesMetaData, by="samples") %>% 
  #group_by(Order, Cluster,Season, Layer, Tree_species) %>%
  group_by(Class, Cluster,Tree_species) %>%
  summarise("reads_count"=sum(counts)) %>% 
  ungroup() %>%
  group_by(Cluster) %>% mutate("reads_percent"= (reads_count/sum(reads_count))*100 )

# samplesMetaData <- samplesMetaData %>% arrange(Tree_number, Season, Layer, samples)
# t <- samplesMetaData$samples 
# df$samples <- factor(df$samples, levels= t )

mycolors <- mycolorsFullList %>% filter(Name %in% df$Class)
#mycolors <- mycolors %>% arrange(Name)
mycolors <- mycolors %>% arrange(myLevel)
mycolor <- mycolors$color
names(mycolor) <- mycolors$Name

df$Class <- factor(df$Class, levels=mycolors$Name)



pc1 <-  ggplot(data=df,aes(y=reads_percent,x=Cluster, fill=Class))+
  geom_bar(stat="identity")+
  theme_bw()+
  scale_fill_manual(values=mycolor)+
  theme(axis.text.x=element_text(angle = 55,hjust=1#,size=2
  ), axis.title.x = element_blank()
  )+
  
  #facet_wrap(~ Season*Tree_species, scales = "free_x") 
  facet_wrap(~ Tree_species, scales = "free_x") + 
  guides(colour = guide_legend(title.position = "right"))

#5*7


pc1


############################  Order - grouped 

top10Taxo<- data %>%
  #dplyr::select(-starts_with("NR_"),-starts_with("Myc"),-ID,-ctg_length, -sum_contig_count) %>% 
  dplyr::select(starts_with(amplicon), Order) %>%
  group_by(Order) %>% summarise_if(is.numeric, sum) %>%
  pivot_longer(cols=-c(Order), names_to="samples", values_to="counts")%>%
  left_join(., samplesMetaData, by="samples") %>% 
  group_by(Order, Cluster,Tree_species) %>%
  summarise("reads_count"=sum(counts)) %>% 
  group_by(Order) %>%
  summarise( max=max(reads_count), mean=mean(reads_count), median=median(reads_count)) %>%
  slice_max(.,order_by = max,n=12)

df <- data %>%
  #dplyr::select(-starts_with("NR_"),-starts_with("Myc"),-ID,-ctg_length, -sum_contig_count) %>% 
  dplyr::select(starts_with(amplicon), Order) %>%
  mutate(Order = case_when(Order %in% top10Taxo$Order ~ Order, TRUE ~ "Others")) %>% 
  mutate(Order = case_when(Order != "unknown order" ~ Order, TRUE ~ "unidentified")) %>% 
  group_by(Order) %>% summarise_if(is.numeric, sum) %>%
  pivot_longer(cols=-c(Order), names_to="samples", values_to="counts")%>%
  left_join(., samplesMetaData, by="samples") %>% 
  #group_by(Order, Cluster,Season, Layer, Tree_species) %>%
  group_by(Order, Cluster,Tree_species) %>%
  summarise("reads_count"=sum(counts)) %>% 
  ungroup() %>%
  group_by(Cluster) %>% mutate("reads_percent"= (reads_count/sum(reads_count))*100 )

# samplesMetaData <- samplesMetaData %>% arrange(Tree_number, Season, Layer, samples)
# t <- samplesMetaData$samples 
# df$samples <- factor(df$samples, levels= t )

mycolors <- mycolorsFullList %>% filter(Name %in% df$Order)
mycolors <- mycolors %>% arrange(myLevel)
mycolor <- mycolors$color
names(mycolor) <- mycolors$Name

df$Order <- factor(df$Order, levels=mycolors$Name)


po1 <- ggplot(data=df,aes(y=reads_percent,x=Cluster, fill=Order))+
  geom_bar(stat="identity")+
  theme_bw()+
  scale_fill_manual(values=mycolor)+
  theme(axis.text.x=element_text(angle = 55,hjust=1#,size=2
  ), axis.title.x = element_blank())+
  #facet_wrap(~ Season*Tree_species, scales = "free_x") 
  facet_wrap(~ Tree_species, scales = "free_x") 

po1 

# 5*6 -samples; Tree_species


############################### Family - grouped 
top10Taxo<- data %>%
  #dplyr::select(-starts_with("NR_"),-starts_with("Myc"),-ID,-ctg_length, -sum_contig_count) %>% 
  dplyr::select(starts_with(amplicon), Family) %>%
  group_by(Family) %>% summarise_if(is.numeric, sum) %>%
  pivot_longer(cols=-c(Family), names_to="samples", values_to="counts")%>%
  left_join(., samplesMetaData, by="samples") %>% 
  group_by(Family, Cluster,Tree_species) %>%
  summarise("reads_count"=sum(counts)) %>% 
  group_by(Family) %>%
  summarise( max=max(reads_count), mean=mean(reads_count), median=median(reads_count)) %>%
  slice_max(.,order_by = max,n=12)


df <- data %>%
  #dplyr::select(-starts_with("NR_"),-starts_with("Myc"),-ID,-ctg_length, -sum_contig_count) %>% 
  dplyr::select(starts_with(amplicon), Family) %>%
  mutate(Family = case_when(Family %in% top10Taxo$Family ~ Family, TRUE ~ "Others")) %>% 
  mutate(Family = case_when(Family != "unknown family" ~ Family, TRUE ~ "unidentified")) %>% 
  group_by(Family) %>% summarise_if(is.numeric, sum) %>%
  pivot_longer(cols=-c(Family), names_to="samples", values_to="counts")%>%
  left_join(., samplesMetaData, by="samples") %>% 
  #group_by(Order, Cluster,Season, Layer, Tree_species) %>%
  group_by(Family, Cluster,Tree_species) %>%
  summarise("reads_count"=sum(counts)) %>% 
  ungroup() %>%
  group_by(Cluster) %>% mutate("reads_percent"= (reads_count/sum(reads_count))*100 )

# samplesMetaData <- samplesMetaData %>% arrange(Tree_number, Season, Layer, samples)
# t <- samplesMetaData$samples 
# df$samples <- factor(df$samples, levels= t )

mycolors <- mycolorsFullList %>% filter(Name %in% df$Family)
#mycolors <- mycolors %>% arrange(Name)
mycolors <- mycolors %>% arrange(myLevel)
mycolor <- mycolors$color
names(mycolor) <- mycolors$Name

df$Family <- factor(df$Family, levels=mycolors$Name)



pf1 <- ggplot(data=df,aes(y=reads_percent,x=Cluster, fill=Family))+
  geom_bar(stat="identity")+
  theme_bw()+
  scale_fill_manual(values=mycolor)+
  theme(axis.text.x=element_text(angle = 55,hjust=1#,size=2
  ), axis.title.x = element_blank()
  )+
  
  #facet_wrap(~ Season*Tree_species, scales = "free_x") 
  facet_wrap(~ Tree_species, scales = "free_x") + 
  guides(colour = guide_legend(title.position = "right"))

#5*7

pf1
############################### Genus - grouped 
top10Taxo<- data %>%
  #dplyr::select(-starts_with("NR_"),-starts_with("Myc"),-ID,-ctg_length, -sum_contig_count) %>% 
  dplyr::select(starts_with(amplicon), Genus) %>%
  group_by(Genus) %>% summarise_if(is.numeric, sum) %>%
  pivot_longer(cols=-c(Genus), names_to="samples", values_to="counts") %>%
  left_join(., samplesMetaData, by="samples") %>% 
  group_by(Genus, Cluster,Tree_species) %>%
  summarise("reads_count"=sum(counts)) %>% 
  group_by(Genus) %>%
  summarise( max=max(reads_count), mean=mean(reads_count), median=median(reads_count)) %>%
  slice_max(.,order_by = max,n=12)


df <- data %>%
  #dplyr::select(-starts_with("NR_"),-starts_with("Myc"),-ID,-ctg_length, -sum_contig_count) %>% 
  dplyr::select(starts_with(amplicon), Genus) %>%
  mutate(Genus = case_when(Genus %in% top10Taxo$Genus ~ Genus, TRUE ~ "Others")) %>% 
  mutate(Genus = case_when(Genus != "unknown genus" ~ Genus, TRUE ~ "unidentified")) %>% 
  group_by(Genus) %>% summarise_if(is.numeric, sum) %>%
  pivot_longer(cols=-c(Genus), names_to="samples", values_to="counts")%>%
  left_join(., samplesMetaData, by="samples") %>% 
  #group_by(Order, Cluster,Season, Layer, Tree_species) %>%
  group_by(Genus, Cluster,Tree_species) %>%
  summarise("reads_count"=sum(counts)) %>% 
  ungroup() %>%
  group_by(Cluster) %>% mutate("reads_percent"= (reads_count/sum(reads_count))*100 )

# samplesMetaData <- samplesMetaData %>% arrange(Tree_number, Season, Layer, samples)
# t <- samplesMetaData$samples 
# df$samples <- factor(df$samples, levels= t )

mycolors <- mycolorsFullList %>% filter(Name %in% df$Genus)
#mycolors <- mycolors %>% arrange(Name)
mycolors <- mycolors %>% arrange(myLevel)
mycolor <- mycolors$color
names(mycolor) <- mycolors$Name

df$Genus <- factor(df$Genus, levels=mycolors$Name)



pg1 <- ggplot(data=df,aes(y=reads_percent,x=Cluster, fill=Genus))+
  geom_bar(stat="identity")+
  theme_bw()+
  scale_fill_manual(values=mycolor)+
  theme(axis.text.x=element_text(angle = 55,hjust=1#,size=2
  ), axis.title.x = element_blank()
  )+
  
  #facet_wrap(~ Season*Tree_species, scales = "free_x") 
  facet_wrap(~ Tree_species, scales = "free_x") + 
  guides(colour = guide_legend(title.position = "right"))

#5*6
pg1


########################## Season ###################
########################## Season ###################
########################## Season ###################
########################## Season ###################
########################## Season ###################


############################### Phylum -grouped 
top10Taxo<- data %>%
  #dplyr::select(-starts_with("NR_"),-starts_with("Myc"),-ID,-ctg_length, -sum_contig_count) %>% 
  dplyr::select(starts_with(amplicon), Phylum) %>%
  group_by(Phylum) %>% summarise_if(is.numeric, sum) %>%
  pivot_longer(cols=-c(Phylum), names_to="samples", values_to="counts")%>%
  left_join(., samplesMetaData, by="samples") %>% 
  group_by(Phylum,Season,Tree_species) %>%
  summarise("reads_count"=sum(counts)) %>% 
  group_by(Phylum) %>%
  summarise( max=max(reads_count), mean=mean(reads_count), median=median(reads_count)) %>%
  slice_max(.,order_by = max,n=12)




df <- data %>%
  #dplyr::select(-starts_with("NR_"),-starts_with("Myc"),-ID,-ctg_length, -sum_contig_count) %>% 
  dplyr::select(starts_with(amplicon), Phylum) %>%
  mutate(Phylum = case_when(Phylum %in% top10Taxo$Phylum ~ Phylum, TRUE ~ "Others")) %>% 
  group_by(Phylum) %>% summarise_if(is.numeric, sum) %>%
  pivot_longer(cols=-c(Phylum), names_to="samples", values_to="counts")%>%
  left_join(., samplesMetaData, by="samples") %>% 
  #group_by(Order, Cluster,Season, Layer, Season) %>%
  group_by(Phylum, Tree_species,Season) %>%
  summarise("reads_count"=sum(counts)) %>% 
  ungroup() %>%
  group_by(Season) %>% mutate("reads_percent"= (reads_count/sum(reads_count))*100 )

# samplesMetaData <- samplesMetaData %>% arrange(Tree_number, Season, Layer, samples)
# t <- samplesMetaData$samples 
# df$samples <- factor(df$samples, levels= t )

mycolors <- mycolorsFullList %>% filter(Name %in% df$Phylum)
mycolors <- mycolors %>% arrange(myLevel)
mycolor <- mycolors$color
names(mycolor) <- mycolors$Name
df$Phylum <- factor(df$Phylum, levels=mycolors$Name)


pp2 <- ggplot(data=df,aes(y=reads_percent,x=Season, fill=Phylum))+
  geom_bar(stat="identity")+
  theme_bw()+
  scale_fill_manual(values=mycolor)+
  theme(axis.text.x=element_text(angle = 55,hjust=1#,size=2
  ), axis.title.x = element_blank()
  )+
  
  #facet_wrap(~ Season*Season, scales = "free_x") 
  facet_wrap(~ Tree_species, scales = "free_x") + 
  guides(colour = guide_legend(title.position = "right"))#+
#theme(legend.position="bottom")

pp2

############################### Class - grouped 
top10Taxo<- data %>%
  #dplyr::select(-starts_with("NR_"),-starts_with("Myc"),-ID,-ctg_length, -sum_contig_count) %>% 
  dplyr::select(starts_with(amplicon), Class) %>%
  group_by(Class) %>% summarise_if(is.numeric, sum) %>%
  pivot_longer(cols=-c(Class), names_to="samples", values_to="counts")%>%
  left_join(., samplesMetaData, by="samples") %>% 
  group_by(Class, Tree_species,Season) %>%
  summarise("reads_count"=sum(counts)) %>% 
  group_by(Class) %>%
  summarise( max=max(reads_count), mean=mean(reads_count), median=median(reads_count)) %>%
  slice_max(.,order_by = max,n=12)



df <- data %>%
  #dplyr::select(-starts_with("NR_"),-starts_with("Myc"),-ID,-ctg_length, -sum_contig_count) %>% 
  dplyr::select(starts_with(amplicon), Class) %>%
  mutate(Class = case_when(Class %in% top10Taxo$Class ~ Class, TRUE ~ "Others")) %>% 
  group_by(Class) %>% summarise_if(is.numeric, sum) %>%
  pivot_longer(cols=-c(Class), names_to="samples", values_to="counts")%>%
  left_join(., samplesMetaData, by="samples") %>% 
  #group_by(Order, Cluster,Season, Layer, Season) %>%
  group_by(Class, Tree_species,Season) %>%
  summarise("reads_count"=sum(counts)) %>% 
  ungroup() %>%
  group_by(Season) %>% mutate("reads_percent"= (reads_count/sum(reads_count))*100 )

# samplesMetaData <- samplesMetaData %>% arrange(Tree_number, Season, Layer, samples)
# t <- samplesMetaData$samples 
# df$samples <- factor(df$samples, levels= t )

mycolors <- mycolorsFullList %>% filter(Name %in% df$Class)
#mycolors <- mycolors %>% arrange(Name)
mycolors <- mycolors %>% arrange(myLevel)
mycolor <- mycolors$color
names(mycolor) <- mycolors$Name

df$Class <- factor(df$Class, levels=mycolors$Name)



pc2 <-  ggplot(data=df,aes(y=reads_percent,x=Season, fill=Class))+
  geom_bar(stat="identity")+
  theme_bw()+
  scale_fill_manual(values=mycolor)+
  theme(axis.text.x=element_text(angle = 55,hjust=1#,size=2
  ), axis.title.x = element_blank()
  )+
  
  #facet_wrap(~ Season*Season, scales = "free_x") 
  facet_wrap(~ Tree_species, scales = "free_x") + 
  guides(colour = guide_legend(title.position = "right"))

#5*7


pc2


############################  Order - grouped 

top10Taxo<- data %>%
  #dplyr::select(-starts_with("NR_"),-starts_with("Myc"),-ID,-ctg_length, -sum_contig_count) %>% 
  dplyr::select(starts_with(amplicon), Order) %>%
  group_by(Order) %>% summarise_if(is.numeric, sum) %>%
  pivot_longer(cols=-c(Order), names_to="samples", values_to="counts")%>%
  left_join(., samplesMetaData, by="samples") %>% 
  group_by(Order, Tree_species,Season) %>%
  summarise("reads_count"=sum(counts)) %>% 
  group_by(Order) %>%
  summarise( max=max(reads_count), mean=mean(reads_count), median=median(reads_count)) %>%
  slice_max(.,order_by = max,n=12)


df <- data %>%
  #dplyr::select(-starts_with("NR_"),-starts_with("Myc"),-ID,-ctg_length, -sum_contig_count) %>% 
  dplyr::select(starts_with(amplicon), Order) %>%
  mutate(Order = case_when(Order %in% top10Taxo$Order ~ Order, TRUE ~ "Others")) %>% 
  mutate(Order = case_when(Order != "unknown order" ~ Order, TRUE ~ "unidentified")) %>% 
  group_by(Order) %>% summarise_if(is.numeric, sum) %>%
  pivot_longer(cols=-c(Order), names_to="samples", values_to="counts")%>%
  left_join(., samplesMetaData, by="samples") %>% 
  #group_by(Order, Cluster,Season, Layer, Season) %>%
  group_by(Order, Tree_species,Season) %>%
  summarise("reads_count"=sum(counts)) %>% 
  ungroup() %>%
  group_by(Season) %>% mutate("reads_percent"= (reads_count/sum(reads_count))*100 )

# samplesMetaData <- samplesMetaData %>% arrange(Tree_number, Season, Layer, samples)
# t <- samplesMetaData$samples 
# df$samples <- factor(df$samples, levels= t )

mycolors <- mycolorsFullList %>% filter(Name %in% df$Order)
mycolors <- mycolors %>% arrange(myLevel)
mycolor <- mycolors$color
names(mycolor) <- mycolors$Name

df$Order <- factor(df$Order, levels=mycolors$Name)


po2 <- ggplot(data=df,aes(y=reads_percent,x=Season, fill=Order))+
  geom_bar(stat="identity")+
  theme_bw()+
  scale_fill_manual(values=mycolor)+
  theme(axis.text.x=element_text(angle = 55,hjust=1#,size=2
  ), axis.title.x = element_blank())+
  #facet_wrap(~ Season*Season, scales = "free_x") 
  facet_wrap(~ Tree_species, scales = "free_x") 

po2 

# 5*6 -samples; Season


############################### Family - grouped 
top10Taxo<- data %>%
  #dplyr::select(-starts_with("NR_"),-starts_with("Myc"),-ID,-ctg_length, -sum_contig_count) %>% 
  dplyr::select(starts_with(amplicon), Family) %>%
  group_by(Family) %>% summarise_if(is.numeric, sum) %>%
  pivot_longer(cols=-c(Family), names_to="samples", values_to="counts")%>%
  left_join(., samplesMetaData, by="samples") %>% 
  group_by(Family, Tree_species,Season) %>%
  summarise("reads_count"=sum(counts)) %>% 
  group_by(Family) %>%
  summarise( max=max(reads_count), mean=mean(reads_count), median=median(reads_count)) %>%
  slice_max(.,order_by = max,n=12)


df <- data %>%
  #dplyr::select(-starts_with("NR_"),-starts_with("Myc"),-ID,-ctg_length, -sum_contig_count) %>% 
  dplyr::select(starts_with(amplicon), Family) %>%
  mutate(Family = case_when(Family %in% top10Taxo$Family ~ Family, TRUE ~ "Others")) %>%
  mutate(Family = case_when(Family != "unknown family" ~ Family, TRUE ~ "unidentified")) %>% 
  group_by(Family) %>% summarise_if(is.numeric, sum) %>%
  pivot_longer(cols=-c(Family), names_to="samples", values_to="counts")%>%
  left_join(., samplesMetaData, by="samples") %>% 
  #group_by(Order, Cluster,Season, Layer, Season) %>%
  group_by(Family, Tree_species,Season) %>%
  summarise("reads_count"=sum(counts)) %>% 
  ungroup() %>%
  group_by(Season) %>% mutate("reads_percent"= (reads_count/sum(reads_count))*100 )

# samplesMetaData <- samplesMetaData %>% arrange(Tree_number, Season, Layer, samples)
# t <- samplesMetaData$samples 
# df$samples <- factor(df$samples, levels= t )

mycolors <- mycolorsFullList %>% filter(Name %in% df$Family)
#mycolors <- mycolors %>% arrange(Name)
mycolors <- mycolors %>% arrange(myLevel)
mycolor <- mycolors$color
names(mycolor) <- mycolors$Name

df$Family <- factor(df$Family, levels=mycolors$Name)



pf2 <- ggplot(data=df,aes(y=reads_percent,x=Season, fill=Family))+
  geom_bar(stat="identity")+
  theme_bw()+
  scale_fill_manual(values=mycolor)+
  theme(axis.text.x=element_text(angle = 55,hjust=1#,size=2
  ), axis.title.x = element_blank()
  )+
  
  #facet_wrap(~ Season*Season, scales = "free_x") 
  facet_wrap(~ Tree_species, scales = "free_x") + 
  guides(colour = guide_legend(title.position = "right"))

#5*7

pf2
############################### Genus - grouped 
top10Taxo<- data %>%
  #dplyr::select(-starts_with("NR_"),-starts_with("Myc"),-ID,-ctg_length, -sum_contig_count) %>% 
  dplyr::select(starts_with(amplicon), Genus) %>%
  group_by(Genus) %>% summarise_if(is.numeric, sum) %>%
  pivot_longer(cols=-c(Genus), names_to="samples", values_to="counts") %>%
  left_join(., samplesMetaData, by="samples") %>% 
  group_by(Genus, Tree_species,Season) %>%
  summarise("reads_count"=sum(counts)) %>% 
  group_by(Genus) %>%
  summarise( max=max(reads_count), mean=mean(reads_count), median=median(reads_count)) %>%
  slice_max(.,order_by = max,n=12)


df <- data %>%
  #dplyr::select(-starts_with("NR_"),-starts_with("Myc"),-ID,-ctg_length, -sum_contig_count) %>% 
  dplyr::select(starts_with(amplicon), Genus) %>%
  mutate(Genus = case_when(Genus %in% top10Taxo$Genus ~ Genus, TRUE ~ "Others")) %>% 
  mutate(Genus = case_when(Genus != "unknown genus" ~ Genus, TRUE ~ "unidentified")) %>% 
  group_by(Genus) %>% summarise_if(is.numeric, sum) %>%
  pivot_longer(cols=-c(Genus), names_to="samples", values_to="counts")%>%
  left_join(., samplesMetaData, by="samples") %>% 
  #group_by(Order, Cluster,Season, Layer, Season) %>%
  group_by(Genus, Tree_species,Season) %>%
  summarise("reads_count"=sum(counts)) %>% 
  ungroup() %>%
  group_by(Season) %>% mutate("reads_percent"= (reads_count/sum(reads_count))*100 )

# samplesMetaData <- samplesMetaData %>% arrange(Tree_number, Season, Layer, samples)
# t <- samplesMetaData$samples 
# df$samples <- factor(df$samples, levels= t )

mycolors <- mycolorsFullList %>% filter(Name %in% df$Genus)
#mycolors <- mycolors %>% arrange(Name)
mycolors <- mycolors %>% arrange(myLevel)
mycolor <- mycolors$color
names(mycolor) <- mycolors$Name

df$Genus <- factor(df$Genus, levels=mycolors$Name)



pg2 <- ggplot(data=df,aes(y=reads_percent,x=Season, fill=Genus))+
  geom_bar(stat="identity")+
  theme_bw()+
  scale_fill_manual(values=mycolor)+
  theme(axis.text.x=element_text(angle = 55,hjust=1#,size=2
  ), axis.title.x = element_blank()
  )+
  
  #facet_wrap(~ Season*Season, scales = "free_x") 
  facet_wrap(~ Tree_species, scales = "free_x") + 
  guides(colour = guide_legend(title.position = "right"))

#5*6
pg2



########################## Layer ###################
########################## Layer ###################
########################## Layer ###################
########################## Layer ###################
########################## Layer ###################
########################## AMF ###################



############################### Phylum -grouped 
top10Taxo<- data %>%
  #dplyr::select(-starts_with("NR_"),-starts_with("Myc"),-ID,-ctg_length, -sum_contig_count) %>% 
  dplyr::select(starts_with(amplicon), Phylum) %>%
  group_by(Phylum) %>% summarise_if(is.numeric, sum) %>%
  pivot_longer(cols=-c(Phylum), names_to="samples", values_to="counts")%>%
  left_join(., samplesMetaData, by="samples") %>% 
  group_by(Phylum, Tree_species,Layer) %>%
  summarise("reads_count"=sum(counts)) %>% 
  group_by(Phylum) %>%
  summarise( max=max(reads_count), mean=mean(reads_count), median=median(reads_count)) %>%
  slice_max(.,order_by = max,n=12)




df <- data %>%
  #dplyr::select(-starts_with("NR_"),-starts_with("Myc"),-ID,-ctg_length, -sum_contig_count) %>% 
  dplyr::select(starts_with(amplicon), Phylum) %>%
  mutate(Phylum = case_when(Phylum %in% top10Taxo$Phylum ~ Phylum, TRUE ~ "Others")) %>% 
  group_by(Phylum) %>% summarise_if(is.numeric, sum) %>%
  pivot_longer(cols=-c(Phylum), names_to="samples", values_to="counts")%>%
  left_join(., samplesMetaData, by="samples") %>% 
  #group_by(Order, Cluster,Season, Layer, Layer) %>%
  group_by(Phylum, Tree_species,Layer) %>%
  summarise("reads_count"=sum(counts)) %>% 
  ungroup() %>%
  group_by(Layer) %>% mutate("reads_percent"= (reads_count/sum(reads_count))*100 )

# samplesMetaData <- samplesMetaData %>% arrange(Tree_number, Season, Layer, samples)
# t <- samplesMetaData$samples 
# df$samples <- factor(df$samples, levels= t )

mycolors <- mycolorsFullList %>% filter(Name %in% df$Phylum)
mycolors <- mycolors %>% arrange(myLevel)
mycolor <- mycolors$color
names(mycolor) <- mycolors$Name
df$Phylum <- factor(df$Phylum, levels=mycolors$Name)


ppa <- ggplot(data=df,aes(y=reads_percent,x=Layer, fill=Phylum))+
  geom_bar(stat="identity")+
  theme_bw()+
  scale_fill_manual(values=mycolor)+
  theme(axis.text.x=element_text(angle = 55,hjust=1#,size=2
  ), axis.title.x = element_blank()
  )+
  
  #facet_wrap(~ Season*Layer, scales = "free_x") 
  facet_wrap(~ Tree_species, scales = "free_x") + 
  guides(colour = guide_legend(title.position = "right"))#+
#theme(legend.position="bottom")

ppa

############################### Class - grouped 
top10Taxo<- data %>%
  #dplyr::select(-starts_with("NR_"),-starts_with("Myc"),-ID,-ctg_length, -sum_contig_count) %>% 
  dplyr::select(starts_with(amplicon), Class) %>%
  group_by(Class) %>% summarise_if(is.numeric, sum) %>%
  pivot_longer(cols=-c(Class), names_to="samples", values_to="counts")%>%
  left_join(., samplesMetaData, by="samples") %>% 
  group_by(Class, Tree_species,Layer) %>%
  summarise("reads_count"=sum(counts)) %>% 
  group_by(Class) %>%
  summarise( max=max(reads_count), mean=mean(reads_count), median=median(reads_count)) %>%
  slice_max(.,order_by = max,n=12)



df <- data %>%
  #dplyr::select(-starts_with("NR_"),-starts_with("Myc"),-ID,-ctg_length, -sum_contig_count) %>% 
  dplyr::select(starts_with(amplicon), Class) %>%
  mutate(Class = case_when(Class %in% top10Taxo$Class ~ Class, TRUE ~ "Others")) %>% 
  mutate(Class = case_when(Class != "unknown class" ~ Class, TRUE ~ "unidentified")) %>%
  group_by(Class) %>% summarise_if(is.numeric, sum) %>%
  pivot_longer(cols=-c(Class), names_to="samples", values_to="counts")%>%
  left_join(., samplesMetaData, by="samples") %>% 
  #group_by(Order, Cluster,Season, Layer, Layer) %>%
  group_by(Class, Tree_species,Layer) %>%
  summarise("reads_count"=sum(counts)) %>% 
  ungroup() %>%
  group_by(Layer) %>% mutate("reads_percent"= (reads_count/sum(reads_count))*100 )

# samplesMetaData <- samplesMetaData %>% arrange(Tree_number, Season, Layer, samples)
# t <- samplesMetaData$samples 
# df$samples <- factor(df$samples, levels= t )

mycolors <- mycolorsFullList %>% filter(Name %in% df$Class)
#mycolors <- mycolors %>% arrange(Name)
mycolors <- mycolors %>% arrange(myLevel)
mycolor <- mycolors$color
names(mycolor) <- mycolors$Name

df$Class <- factor(df$Class, levels=mycolors$Name)



pca <-  ggplot(data=df,aes(y=reads_percent,x=Layer, fill=Class))+
  geom_bar(stat="identity")+
  theme_bw()+
  scale_fill_manual(values=mycolor)+
  theme(axis.text.x=element_text(angle = 55,hjust=1#,size=2
  ), axis.title.x = element_blank()
  )+
  
  #facet_wrap(~ Season*Layer, scales = "free_x") 
  facet_wrap(~ Tree_species, scales = "free_x") + 
  guides(colour = guide_legend(title.position = "right"))

#5*7


pca


############################  Order - grouped 

top10Taxo<- data %>%
  #dplyr::select(-starts_with("NR_"),-starts_with("Myc"),-ID,-ctg_length, -sum_contig_count) %>% 
  dplyr::select(starts_with(amplicon), Order) %>%
  group_by(Order) %>% summarise_if(is.numeric, sum) %>%
  pivot_longer(cols=-c(Order), names_to="samples", values_to="counts")%>%
  left_join(., samplesMetaData, by="samples") %>% 
  group_by(Order, Tree_species,Layer) %>%
  summarise("reads_count"=sum(counts)) %>% 
  group_by(Order) %>%
  summarise( max=max(reads_count), mean=mean(reads_count), median=median(reads_count)) %>%
  slice_max(.,order_by = max,n=12)


df <- data %>%
  #dplyr::select(-starts_with("NR_"),-starts_with("Myc"),-ID,-ctg_length, -sum_contig_count) %>% 
  dplyr::select(starts_with(amplicon), Order) %>%
  mutate(Order = case_when(Order %in% top10Taxo$Order ~ Order, TRUE ~ "Others")) %>% 
  mutate(Order = case_when(Order != "unknown order" ~ Order, TRUE ~ "unidentified")) %>% 
  group_by(Order) %>% summarise_if(is.numeric, sum) %>%
  pivot_longer(cols=-c(Order), names_to="samples", values_to="counts")%>%
  left_join(., samplesMetaData, by="samples") %>% 
  #group_by(Order, Cluster,Season, Layer, Layer) %>%
  group_by(Order, Tree_species,Layer) %>%
  summarise("reads_count"=sum(counts)) %>% 
  ungroup() %>%
  group_by(Layer) %>% mutate("reads_percent"= (reads_count/sum(reads_count))*100 )

# samplesMetaData <- samplesMetaData %>% arrange(Tree_number, Season, Layer, samples)
# t <- samplesMetaData$samples 
# df$samples <- factor(df$samples, levels= t )

mycolors <- mycolorsFullList %>% filter(Name %in% df$Order)
mycolors <- mycolors %>% arrange(myLevel)
mycolor <- mycolors$color
names(mycolor) <- mycolors$Name

df$Order <- factor(df$Order, levels=mycolors$Name)


poa <- ggplot(data=df,aes(y=reads_percent,x=Layer, fill=Order))+
  geom_bar(stat="identity")+
  theme_bw()+
  scale_fill_manual(values=mycolor)+
  theme(axis.text.x=element_text(angle = 55,hjust=1#,size=2
  ), axis.title.x = element_blank())+
  #facet_wrap(~ Season*Layer, scales = "free_x") 
  facet_wrap(~ Tree_species, scales = "free_x") 

poa 

# 5*6 -samples; Layer


############################### Family - grouped 
top10Taxo<- data %>%
  #dplyr::select(-starts_with("NR_"),-starts_with("Myc"),-ID,-ctg_length, -sum_contig_count) %>% 
  dplyr::select(starts_with(amplicon), Family) %>%
  group_by(Family) %>% summarise_if(is.numeric, sum) %>%
  pivot_longer(cols=-c(Family), names_to="samples", values_to="counts")%>%
  left_join(., samplesMetaData, by="samples") %>% 
  group_by(Family, Tree_species,Layer) %>%
  summarise("reads_count"=sum(counts)) %>% 
  group_by(Family) %>%
  summarise( max=max(reads_count), mean=mean(reads_count), median=median(reads_count)) %>%
  slice_max(.,order_by = max,n=12)


df <- data %>%
  #dplyr::select(-starts_with("NR_"),-starts_with("Myc"),-ID,-ctg_length, -sum_contig_count) %>% 
  dplyr::select(starts_with(amplicon), Family) %>%
  mutate(Family = case_when(Family %in% top10Taxo$Family ~ Family, TRUE ~ "Others")) %>% 
  mutate(Family = case_when(Family  != "unknown family" ~ Family , TRUE ~ "unidentified")) %>% 
  group_by(Family) %>% summarise_if(is.numeric, sum) %>%
  pivot_longer(cols=-c(Family), names_to="samples", values_to="counts")%>%
  left_join(., samplesMetaData, by="samples") %>% 
  #group_by(Order, Cluster,Season, Layer, Layer) %>%
  group_by(Family, Tree_species,Layer) %>%
  summarise("reads_count"=sum(counts)) %>% 
  ungroup() %>%
  group_by(Layer) %>% mutate("reads_percent"= (reads_count/sum(reads_count))*100 )

# samplesMetaData <- samplesMetaData %>% arrange(Tree_number, Season, Layer, samples)
# t <- samplesMetaData$samples 
# df$samples <- factor(df$samples, levels= t )

mycolors <- mycolorsFullList %>% filter(Name %in% df$Family)
#mycolors <- mycolors %>% arrange(Name)
mycolors <- mycolors %>% arrange(myLevel)
mycolor <- mycolors$color
names(mycolor) <- mycolors$Name

df$Family <- factor(df$Family, levels=mycolors$Name)



pfa <- ggplot(data=df,aes(y=reads_percent,x=Layer, fill=Family))+
  geom_bar(stat="identity")+
  theme_bw()+
  scale_fill_manual(values=mycolor)+
  theme(axis.text.x=element_text(angle = 55,hjust=1#,size=2
  ), axis.title.x = element_blank()
  )+
  
  #facet_wrap(~ Season*Layer, scales = "free_x") 
  facet_wrap(~ Tree_species, scales = "free_x") + 
  guides(colour = guide_legend(title.position = "right"))

#5*7

pfa
############################### Genus - grouped 
top10Taxo<- data %>%
  #dplyr::select(-starts_with("NR_"),-starts_with("Myc"),-ID,-ctg_length, -sum_contig_count) %>% 
  dplyr::select(starts_with(amplicon), Genus) %>%
  group_by(Genus) %>% summarise_if(is.numeric, sum) %>%
  pivot_longer(cols=-c(Genus), names_to="samples", values_to="counts") %>%
  left_join(., samplesMetaData, by="samples") %>% 
  group_by(Genus, Tree_species,Layer) %>%
  summarise("reads_count"=sum(counts)) %>% 
  group_by(Genus) %>%
  summarise( max=max(reads_count), mean=mean(reads_count), median=median(reads_count)) %>%
  slice_max(.,order_by = max,n=12)


df <- data %>%
  #dplyr::select(-starts_with("NR_"),-starts_with("Myc"),-ID,-ctg_length, -sum_contig_count) %>% 
  dplyr::select(starts_with(amplicon), Genus) %>%
  mutate(Genus = case_when(Genus %in% top10Taxo$Genus ~ Genus, TRUE ~ "Others")) %>% 
  mutate(Genus = case_when(Genus != "unknown genus" ~ Genus, TRUE ~ "unidentified")) %>% 
  group_by(Genus) %>% summarise_if(is.numeric, sum) %>%
  pivot_longer(cols=-c(Genus), names_to="samples", values_to="counts")%>%
  left_join(., samplesMetaData, by="samples") %>% 
  #group_by(Order, Cluster,Season, Layer, Layer) %>%
  group_by(Genus, Tree_species,Layer) %>%
  summarise("reads_count"=sum(counts)) %>% 
  ungroup() %>%
  group_by(Layer) %>% mutate("reads_percent"= (reads_count/sum(reads_count))*100 )

# samplesMetaData <- samplesMetaData %>% arrange(Tree_number, Season, Layer, samples)
# t <- samplesMetaData$samples 
# df$samples <- factor(df$samples, levels= t )

mycolors <- mycolorsFullList %>% filter(Name %in% df$Genus)
#mycolors <- mycolors %>% arrange(Name)
mycolors <- mycolors %>% arrange(myLevel)
mycolor <- mycolors$color
names(mycolor) <- mycolors$Name

df$Genus <- factor(df$Genus, levels=mycolors$Name)



pga <- ggplot(data=df,aes(y=reads_percent,x=Layer, fill=Genus))+
  geom_bar(stat="identity")+
  theme_bw()+
  scale_fill_manual(values=mycolor)+
  theme(axis.text.x=element_text(angle = 55,hjust=1#,size=2
  ), axis.title.x = element_blank()
  )+
  
  #facet_wrap(~ Season*Layer, scales = "free_x") 
  facet_wrap(~ Tree_species, scales = "free_x") + 
  guides(colour = guide_legend(title.position = "right"))

#5*6
pga


library(patchwork)

#(pp1 | pc1 | po1 | pf1 | pg1 | pe1) / (pp2 | pc2 | po2 | pf2 | pg2 | pe2) / (ppa | pca | poa | pfa | pga | pea)
(pp1 | pp2 | ppa) / (pc1 | pc2 | pca) / (po1 | po2 | poa) / (pf1 | pf2 | pfa) /(pg1|pg2|pga)
#40*20