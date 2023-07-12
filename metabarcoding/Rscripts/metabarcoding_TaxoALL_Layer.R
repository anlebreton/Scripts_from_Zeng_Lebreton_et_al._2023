# ## Import functions
library("plyr");library("dplyr");library("ggplot2");library("tidyr");library("tidyverse")
library(pheatmap)
library(RColorBrewer)


samplesMetaData <- read.csv("../0_Metadata/Yulong_sampling_ITS1_metadata.txt", sep='\t')
samplesMetaData <- samplesMetaData %>% mutate(samples= gsub("-", ".", X)) 
data <- read.csv("../1_data/Yulong_ITS1_normalised_abundance_addCols.tsv", sep='\t')
#data[ data == "unidentified" ] <- "Unidentified" 
#data[ data == "Multi-affiliation" ] <- "Unidentified" 

mycolorsFullList <- read.csv("../..//mycolorFile.txt",sep="\t")
mycolorVector <- mycolorsFullList$color
names(mycolorVector) <- mycolorsFullList$Name

#my_list <- list(as.vector(mycolorsFullList$color),as.vector(mycolorsFullList$Name) )
#names(my_list) <- c(as.vector(mycolorsFullList$Name))


#############################

df <- data %>%
  #dplyr::select(-starts_with("NR_"),-starts_with("Myc"),-ID,-ctg_length, -sum_contig_count) %>% 
  dplyr::select(starts_with("ITS1"), Phylum)%>%
  group_by(Phylum) %>% summarise_if(is.numeric, sum) %>%
  pivot_longer(cols=-c(Phylum), names_to="samples", values_to="counts")%>%
  group_by(Phylum) 

df <- merge(df, samplesMetaData, by="samples")
samplesMetaData <- samplesMetaData %>% arrange(Tree_number, Season, Season, samples)
t <- samplesMetaData$samples 
df$samples <- factor(df$samples, levels= t )

mycolors <- mycolorsFullList %>% filter(Name %in% df$Phylum)
mycolor <- mycolors$color
names(mycolor) <- mycolors$Name


ggplot(data=df,aes(y=counts,x=samples, fill=Phylum))+
  geom_bar(stat="identity")+
  theme_bw()+
  scale_fill_manual(values=mycolor)+
  theme(axis.text.x=element_text(angle = 55,hjust=1,size=2))+
  #facet_wrap(~ Season*Tree_species, scales = "free_x") 
  facet_wrap(~ Tree_species, scales = "free_x") 


# 5 * 15 


############################ CLASS ##########################

top10Taxo<- data %>%
  #dplyr::select(-starts_with("NR_"),-starts_with("Myc"),-ID,-ctg_length, -sum_contig_count) %>% 
  dplyr::select(starts_with("ITS1"), Class) %>%
  group_by(Class) %>% summarise_if(is.numeric, sum) %>%
  pivot_longer(cols=-c(Class), names_to="samples", values_to="counts")%>%
  group_by(Class) %>%
  summarise( max=max(counts), mean=mean(counts), median=median(counts)) %>%
  slice_max(.,order_by = max,n=11)

df <- data %>%
  #dplyr::select(-starts_with("NR_"),-starts_with("Myc"),-ID,-ctg_length, -sum_contig_count) %>% 
  dplyr::select(starts_with("ITS1"), Class) %>%
  mutate(Class = case_when(Class %in% top10Taxo$Class ~ Class, TRUE ~ "Others")) %>% 
  group_by(Class) %>% summarise_if(is.numeric, sum) %>%
  pivot_longer(cols=-c(Class), names_to="samples", values_to="counts")%>%
  group_by(Class) 

df <- merge(df, samplesMetaData, by="samples")
samplesMetaData <- samplesMetaData %>% arrange(Tree_number, Season, Season, samples)
t <- samplesMetaData$samples 
df$samples <- factor(df$samples, levels= t )

mycolors <- mycolorsFullList %>% filter(Name %in% df$Class)
mycolor <- mycolors$color
names(mycolor) <- mycolors$Name

ggplot(data=df,aes(y=counts,x=samples, fill=Class))+
  geom_bar(stat="identity")+
  theme_bw()+
  scale_fill_manual(values=mycolor)+
  theme(axis.text.x=element_text(angle = 55,hjust=1,size=2))+
  #facet_wrap(~ Season*Tree_species, scales = "free_x") 
  facet_wrap(~ Tree_species, scales = "free_x") 

# 5 * 15 


############################ FAMILY ##########################

top10Taxo<- data %>%
  #dplyr::select(-starts_with("NR_"),-starts_with("Myc"),-ID,-ctg_length, -sum_contig_count) %>% 
  dplyr::select(starts_with("ITS1"), Family) %>%
  group_by(Family) %>% summarise_if(is.numeric, sum) %>%
  pivot_longer(cols=-c(Family), names_to="samples", values_to="counts")%>%
  group_by(Family) %>%
  summarise( max=max(counts), mean=mean(counts), median=median(counts)) %>%
  slice_max(.,order_by = max,n=11)

df <- data %>%
  #dplyr::select(-starts_with("NR_"),-starts_with("Myc"),-ID,-ctg_length, -sum_contig_count) %>% 
  dplyr::select(starts_with("ITS1"), Family) %>%
  mutate(Family = case_when(Family %in% top10Taxo$Family ~ Family, TRUE ~ "Others")) %>% 
  group_by(Family) %>% summarise_if(is.numeric, sum) %>%
  pivot_longer(cols=-c(Family), names_to="samples", values_to="counts")%>%
  group_by(Family) 

df <- merge(df, samplesMetaData, by="samples")
samplesMetaData <- samplesMetaData %>% arrange(Tree_number, Season, Season, samples)
t <- samplesMetaData$samples 
df$samples <- factor(df$samples, levels= t )

mycolors <- mycolorsFullList %>% filter(Name %in% df$Family)
mycolors <- mycolors %>% arrange(Name)
mycolor <- mycolors$color
names(mycolor) <- mycolors$Name

ggplot(data=df,aes(y=counts,x=samples, fill=Family))+
  geom_bar(stat="identity")+
  theme_bw()+
  scale_fill_manual(values=mycolor)+
  theme(axis.text.x=element_text(angle = 55,hjust=1,size=2))+
  #facet_wrap(~ Season*Tree_species, scales = "free_x") 
  facet_wrap(~ Tree_species, scales = "free_x") 

# 5 * 15 


############################





top10Taxo<- data %>%
  #dplyr::select(-starts_with("NR_"),-starts_with("Myc"),-ID,-ctg_length, -sum_contig_count) %>% 
  dplyr::select(starts_with("ITS1"), Order) %>%
  group_by(Order) %>% summarise_if(is.numeric, sum) %>%
  pivot_longer(cols=-c(Order), names_to="samples", values_to="counts")%>%
  group_by(Order) %>%
  summarise( max=max(counts), mean=mean(counts), median=median(counts)) %>%
  slice_max(.,order_by = max,n=12)

c(rev(brewer.pal(7,"Blues")))
  
df <- data %>%
  #dplyr::select(-starts_with("NR_"),-starts_with("Myc"),-ID,-ctg_length, -sum_contig_count) %>% 
  dplyr::select(starts_with("ITS1"), Order) %>%
  mutate(Order = case_when(Order %in% top10Taxo$Order ~ Order, TRUE ~ "Others")) %>% 
  group_by(Order) %>% summarise_if(is.numeric, sum) %>%
  pivot_longer(cols=-c(Order), names_to="samples", values_to="counts")%>%
  group_by(Order) 


df <- merge(df, samplesMetaData, by="samples")
samplesMetaData <- samplesMetaData %>% arrange(Tree_number, Season, Season, samples)
t <- samplesMetaData$samples 
df$samples <- factor(df$samples, levels= t )

mycolors <- mycolorsFullList %>% filter(Name %in% df$Order)
mycolor <- mycolors$color
names(mycolor) <- mycolors$Name



ggplot(data=df,aes(y=counts,x=samples, fill=Order))+
  geom_bar(stat="identity")+
  theme_bw()+
  scale_fill_manual(values=mycolor)+
  theme(axis.text.x=element_text(angle = 55,hjust=1,size=2))+
  facet_wrap(~ Season*Tree_species, scales = "free_x") 
  #facet_wrap(~ Tree_species, scales = "free_x") 


# 8*10 - samples;  Season*Tree_species
# 5*15 -samples; Tree_species


#?hclust()

# library(cowplot)
# h1 <- ggplot(data=df,aes(y=counts,x=samples, fill=Order))+
#   geom_bar(stat="identity")+
#   theme_bw()+
#   scale_fill_manual(values=mycolor)+
#   theme(axis.text.x=element_text(angle = 55,hjust=1))+
#   facet_wrap(~ Tree_species, scales = "free_x")
# 
# h2 <- ggplot(data=df, aes(x = samples, y = 1, fill = Season))+
#   geom_bar(stat = "identity", width = 1)+
#   theme_void()+
#   theme(panel.spacing.x = unit(1, "mm"))+
#   facet_wrap(.~Tree_species, scales = "free_x")
# 
# h3 <- ggplot(data=df, aes(x = samples, y = 1, fill = Season))+
#   geom_bar(stat = "identity", width = 1)+
#   theme_void()+
#   theme(panel.spacing.x = unit(1, "mm"))+
#   facet_wrap(.~Tree_species, scales = "free_x")
# 
# h4 <- ggplot(data=df, aes(x = samples, y = 1, fill = Season))+
#   geom_bar(stat = "identity", width = 1)+
#   theme_void()+
#   theme(panel.spacing.x = unit(1, "mm"))+
#   facet_wrap(.~Tree_species, scales = "free_x")
# 
# 
# 
# legend <- plot_grid(get_legend(h2), get_legend(h1), ncol = 1)
# h1 <- h1 + theme(legend.position = "none")
# h2 <- h2 + theme(legend.position = "none")
# h3 <- h3 + theme(legend.position = "none")
# h4 <- h4 + theme(legend.position = "none")
# #plot <- plot_grid(h2, h1, align = "v", ncol = 1, axis = "b", rel_heights = c(0.5, 15))
# plot <- plot_grid(h1, h2,h3,h4, align = "v", ncol = 1, axis = "b", rel_heights = c(15, 0.5,0.5,0.5))
# 
# plot_grid(plot, legend, nrow = 1, rel_widths = c(10, 1.5))




# 
# 
# 
# scale_fill_gradient(low="white", high = "#dfc27d") +
#   new_scale_fill() +
#   geom_fruit(
#     data=plotAbundance,
#     geom=geom_bar,
#     mapping = aes(y=SpeciesID,x=nb, fill=plot_name),
#     #colour="grey50",
#     pwidth=0.8, 
#     orientation="y", 
#     stat="identity",
#     offset =0.075  )
# 
# 
# 
# 
# tmp <- df
# df$Order <- NULL
# rownames(df) <- tmp$Order
# 
# pheatmap(df )

############################  Order - grouped 

top10Taxo<- data %>%
  #dplyr::select(-starts_with("NR_"),-starts_with("Myc"),-ID,-ctg_length, -sum_contig_count) %>% 
  dplyr::select(starts_with("ITS1"), Order) %>%
  group_by(Order) %>% summarise_if(is.numeric, sum) %>%
  pivot_longer(cols=-c(Order), names_to="samples", values_to="counts")%>%
  group_by(Order) %>%
  summarise( max=max(counts), mean=mean(counts), median=median(counts)) %>%
  slice_max(.,order_by = max,n=12)


df <- data %>%
  #dplyr::select(-starts_with("NR_"),-starts_with("Myc"),-ID,-ctg_length, -sum_contig_count) %>% 
  dplyr::select(starts_with("ITS1"), Order) %>%
  mutate(Order = case_when(Order %in% top10Taxo$Order ~ Order, TRUE ~ "Others")) %>% 
  group_by(Order) %>% summarise_if(is.numeric, sum) %>%
  pivot_longer(cols=-c(Order), names_to="samples", values_to="counts")%>%
  left_join(., samplesMetaData, by="samples") %>% 
  #group_by(Order, Season,Season, Season, Tree_species) %>%
  group_by(Order, Season,Tree_species) %>%
  summarise("reads_count"=sum(counts)) %>% 
  ungroup() %>%
  group_by(Season,Tree_species) %>% mutate("reads_percent"= (reads_count/sum(reads_count))*100 )

# samplesMetaData <- samplesMetaData %>% arrange(Tree_number, Season, Season, samples)
# t <- samplesMetaData$samples 
# df$samples <- factor(df$samples, levels= t )

mycolors <- mycolorsFullList %>% filter(Name %in% df$Order)
mycolors <- mycolors %>% arrange(myLevel)
mycolor <- mycolors$color
names(mycolor) <- mycolors$Name

df$Order <- factor(df$Order, levels=mycolors$Name)


ggplot(data=df,aes(y=reads_percent,x=Season, fill=Order))+
  geom_bar(stat="identity")+
  theme_bw()+
  scale_fill_manual(values=mycolor)+
  theme(axis.text.x=element_text(angle = 55,hjust=1#,size=2
                                 ), axis.title.x = element_blank())+
  #facet_wrap(~ Season*Tree_species, scales = "free_x") 
facet_wrap(~ Tree_species, scales = "free_x") 


# 5*6 -samples; Tree_species

############################### Phylum -grouped 
top10Taxo<- data %>%
  #dplyr::select(-starts_with("NR_"),-starts_with("Myc"),-ID,-ctg_length, -sum_contig_count) %>% 
  dplyr::select(starts_with("ITS1"), Phylum) %>%
  group_by(Phylum) %>% summarise_if(is.numeric, sum) %>%
  pivot_longer(cols=-c(Phylum), names_to="samples", values_to="counts")%>%
  group_by(Phylum) %>%
  summarise( max=max(counts), mean=mean(counts), median=median(counts)) %>%
  slice_max(.,order_by = max,n=3)


df <- data %>%
  #dplyr::select(-starts_with("NR_"),-starts_with("Myc"),-ID,-ctg_length, -sum_contig_count) %>% 
  dplyr::select(starts_with("ITS1"), Phylum) %>%
  mutate(Phylum = case_when(Phylum %in% top10Taxo$Phylum ~ Phylum, TRUE ~ "Others")) %>% 
  group_by(Phylum) %>% summarise_if(is.numeric, sum) %>%
  pivot_longer(cols=-c(Phylum), names_to="samples", values_to="counts")%>%
  left_join(., samplesMetaData, by="samples") %>% 
  #group_by(Order, Season,Season, Season, Tree_species) %>%
  group_by(Phylum, Season,Tree_species) %>%
  summarise("reads_count"=sum(counts)) %>% 
  ungroup() %>%
  group_by(Season,Tree_species) %>% mutate("reads_percent"= (reads_count/sum(reads_count))*100 )

# samplesMetaData <- samplesMetaData %>% arrange(Tree_number, Season, Season, samples)
# t <- samplesMetaData$samples 
# df$samples <- factor(df$samples, levels= t )

mycolors <- mycolorsFullList %>% filter(Name %in% df$Phylum)
mycolors <- mycolors %>% arrange(Name)
mycolor <- mycolors$color
names(mycolor) <- mycolors$Name


ggplot(data=df,aes(y=reads_percent,x=Season, fill=Phylum))+
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


############################### Class - grouped 
top10Taxo<- data %>%
  #dplyr::select(-starts_with("NR_"),-starts_with("Myc"),-ID,-ctg_length, -sum_contig_count) %>% 
  dplyr::select(starts_with("ITS1"), Class) %>%
  group_by(Class) %>% summarise_if(is.numeric, sum) %>%
  pivot_longer(cols=-c(Class), names_to="samples", values_to="counts")%>%
  group_by(Class) %>%
  summarise( max=max(counts), mean=mean(counts), median=median(counts)) %>%
  slice_max(.,order_by = max,n=6)


df <- data %>%
  #dplyr::select(-starts_with("NR_"),-starts_with("Myc"),-ID,-ctg_length, -sum_contig_count) %>% 
  dplyr::select(starts_with("ITS1"), Class) %>%
  mutate(Class = case_when(Class %in% top10Taxo$Class ~ Class, TRUE ~ "Others")) %>% 
  group_by(Class) %>% summarise_if(is.numeric, sum) %>%
  pivot_longer(cols=-c(Class), names_to="samples", values_to="counts")%>%
  left_join(., samplesMetaData, by="samples") %>% 
  #group_by(Order, Season,Season, Season, Tree_species) %>%
  group_by(Class, Season,Tree_species) %>%
  summarise("reads_count"=sum(counts)) %>% 
  ungroup() %>%
  group_by(Season,Tree_species) %>% mutate("reads_percent"= (reads_count/sum(reads_count))*100 )

# samplesMetaData <- samplesMetaData %>% arrange(Tree_number, Season, Season, samples)
# t <- samplesMetaData$samples 
# df$samples <- factor(df$samples, levels= t )

mycolors <- mycolorsFullList %>% filter(Name %in% df$Class)
#mycolors <- mycolors %>% arrange(Name)
mycolors <- mycolors %>% arrange(myLevel)
mycolor <- mycolors$color
names(mycolor) <- mycolors$Name

df$Class <- factor(df$Class, levels=mycolors$Name)



ggplot(data=df,aes(y=reads_percent,x=Season, fill=Class))+
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




############################### Family - grouped 
top10Taxo<- data %>%
  #dplyr::select(-starts_with("NR_"),-starts_with("Myc"),-ID,-ctg_length, -sum_contig_count) %>% 
  dplyr::select(starts_with("ITS1"), Family) %>%
  group_by(Family) %>% summarise_if(is.numeric, sum) %>%
  pivot_longer(cols=-c(Family), names_to="samples", values_to="counts")%>%
  group_by(Family) %>%
  summarise( max=max(counts), mean=mean(counts), median=median(counts)) %>%
  slice_max(.,order_by = max,n=11)


df <- data %>%
  #dplyr::select(-starts_with("NR_"),-starts_with("Myc"),-ID,-ctg_length, -sum_contig_count) %>% 
  dplyr::select(starts_with("ITS1"), Family) %>%
  mutate(Family = case_when(Family %in% top10Taxo$Family ~ Family, TRUE ~ "Others")) %>% 
  group_by(Family) %>% summarise_if(is.numeric, sum) %>%
  pivot_longer(cols=-c(Family), names_to="samples", values_to="counts")%>%
  left_join(., samplesMetaData, by="samples") %>% 
  #group_by(Order, Season,Season, Season, Tree_species) %>%
  group_by(Family, Season,Tree_species) %>%
  summarise("reads_count"=sum(counts)) %>% 
  ungroup() %>%
  group_by(Season,Tree_species) %>% mutate("reads_percent"= (reads_count/sum(reads_count))*100 )

# samplesMetaData <- samplesMetaData %>% arrange(Tree_number, Season, Season, samples)
# t <- samplesMetaData$samples 
# df$samples <- factor(df$samples, levels= t )

mycolors <- mycolorsFullList %>% filter(Name %in% df$Family)
#mycolors <- mycolors %>% arrange(Name)
mycolors <- mycolors %>% arrange(myLevel)
mycolor <- mycolors$color
names(mycolor) <- mycolors$Name

df$Family <- factor(df$Family, levels=mycolors$Name)



ggplot(data=df,aes(y=reads_percent,x=Season, fill=Family))+
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


############################### Genus - grouped 
top10Taxo<- data %>%
  #dplyr::select(-starts_with("NR_"),-starts_with("Myc"),-ID,-ctg_length, -sum_contig_count) %>% 
  dplyr::select(starts_with("ITS1"), Genus) %>%
  group_by(Genus) %>% summarise_if(is.numeric, sum) %>%
  pivot_longer(cols=-c(Genus), names_to="samples", values_to="counts")%>%
  group_by(Genus) %>%
  summarise( max=max(counts), mean=mean(counts), median=median(counts)) %>%
  slice_max(.,order_by = max,n=12)


df <- data %>%
  #dplyr::select(-starts_with("NR_"),-starts_with("Myc"),-ID,-ctg_length, -sum_contig_count) %>% 
  dplyr::select(starts_with("ITS1"), Genus) %>%
  mutate(Genus = case_when(Genus %in% top10Taxo$Genus ~ Genus, TRUE ~ "Others")) %>% 
  group_by(Genus) %>% summarise_if(is.numeric, sum) %>%
  pivot_longer(cols=-c(Genus), names_to="samples", values_to="counts")%>%
  left_join(., samplesMetaData, by="samples") %>% 
  #group_by(Order, Season,Season, Season, Tree_species) %>%
  group_by(Genus, Season,Tree_species) %>%
  summarise("reads_count"=sum(counts)) %>% 
  ungroup() %>%
  group_by(Season,Tree_species) %>% mutate("reads_percent"= (reads_count/sum(reads_count))*100 )

# samplesMetaData <- samplesMetaData %>% arrange(Tree_number, Season, Season, samples)
# t <- samplesMetaData$samples 
# df$samples <- factor(df$samples, levels= t )

mycolors <- mycolorsFullList %>% filter(Name %in% df$Genus)
#mycolors <- mycolors %>% arrange(Name)
mycolors <- mycolors %>% arrange(myLevel)
mycolor <- mycolors$color
names(mycolor) <- mycolors$Name

df$Genus <- factor(df$Genus, levels=mycolors$Name)



ggplot(data=df,aes(y=reads_percent,x=Season, fill=Genus))+
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


############################## ECOLOGY ##########################

ecology <-read.csv(file="../../ecology_fungaltrait_genus.txt", sep = "\t")
#ecology <- ecology %>% mutate(primary_lifestyle = case_when((!is.na(primary_lifestyle)) ~ primary_lifestyle, FALSE ~ "Unidentified"))

mycolorsFullList <- read.csv("../..//mycolorFile.txt",sep="\t")

top10Taxo<- data %>%
  left_join(.,ecology) %>%
  #dplyr::select(-starts_with("NR_"),-starts_with("Myc"),-ID,-ctg_length, -sum_contig_count) %>% 
  dplyr::select(starts_with("ITS1"), primary_lifestyle) %>%
  group_by(primary_lifestyle) %>% summarise_if(is.numeric, sum) %>%
  pivot_longer(cols=-c(primary_lifestyle), names_to="samples", values_to="counts")%>%
  group_by(primary_lifestyle) %>%
  summarise( max=max(counts), mean=mean(counts), median=median(counts)) %>%
  slice_max(.,order_by = max,n=11)


df <- data %>%
  #dplyr::select(-starts_with("NR_"),-starts_with("Myc"),-ID,-ctg_length, -sum_contig_count) %>% 
  dplyr::select(starts_with("ITS1"), Genus) %>%
  group_by(Genus) %>% summarise_if(is.numeric, sum) %>%
  pivot_longer(cols=-c(Genus), names_to="samples", values_to="counts")%>%
  left_join(., samplesMetaData, by="samples") %>% 
  left_join(., ecology, by="Genus") %>%
  mutate(primary_lifestyle = case_when(primary_lifestyle %in% top10Taxo$primary_lifestyle ~ primary_lifestyle, TRUE ~ "Others")) %>% 
  group_by(samples,primary_lifestyle, Season,Tree_species) %>%
  summarise("reads_count"=sum(counts)) %>% 
  ungroup() %>%
  group_by(samples) %>% mutate("reads_percent"= (reads_count/sum(reads_count))*100 ) 

df$primary_lifestyle[is.na(df$primary_lifestyle)] <- "Unidentified"
#df <- df %>% filter(primary_lifestyle != "Unidentified")

mycolors <- mycolorsFullList %>% filter(Name %in% df$primary_lifestyle)
mycolors <- mycolors %>% arrange(myLevel)
mycolor <- mycolors$color
names(mycolor) <- mycolors$Name

df$primary_lifestyle <- factor(df$primary_lifestyle, levels=mycolors$Name)


ggplot(data=df,aes(y=reads_percent,x=samples, fill=primary_lifestyle))+
  geom_bar(stat="identity")+
  theme_bw()+
  scale_fill_manual(values=mycolor)+
  theme(axis.text.x=element_text(angle = 55,hjust=1,size=2
  ), axis.title.x = element_blank()
  )+
  
  #facet_wrap(~ Season*Tree_species, scales = "free_x") 
  facet_wrap(~ Tree_species, scales = "free_x") + 
  guides(colour = guide_legend(title.position = "right"))


# 5 * 15 

############################## ECOLOGY ##########################

ecology <-read.csv(file="../../ecology_fungaltrait_genus.txt", sep = "\t")
#ecology <- ecology %>% mutate(primary_lifestyle = case_when((!is.na(primary_lifestyle)) ~ primary_lifestyle, FALSE ~ "Unidentified"))

mycolorsFullList <- read.csv("../..//mycolorFile.txt",sep="\t")

top10Taxo<- data %>%
  left_join(.,ecology) %>%
  #dplyr::select(-starts_with("NR_"),-starts_with("Myc"),-ID,-ctg_length, -sum_contig_count) %>% 
  dplyr::select(starts_with("ITS1"), primary_lifestyle) %>%
  group_by(primary_lifestyle) %>% summarise_if(is.numeric, sum) %>%
  pivot_longer(cols=-c(primary_lifestyle), names_to="samples", values_to="counts")%>%
  group_by(primary_lifestyle) %>%
  summarise( max=max(counts), mean=mean(counts), median=median(counts)) %>%
  slice_max(.,order_by = max,n=8)


df <- data %>%
  #dplyr::select(-starts_with("NR_"),-starts_with("Myc"),-ID,-ctg_length, -sum_contig_count) %>% 
  dplyr::select(starts_with("ITS1"), Genus) %>%
  group_by(Genus) %>% summarise_if(is.numeric, sum) %>%
  pivot_longer(cols=-c(Genus), names_to="samples", values_to="counts")%>%
  left_join(., samplesMetaData, by="samples") %>% 
  left_join(., ecology, by="Genus") %>%
  mutate(primary_lifestyle = case_when(primary_lifestyle %in% top10Taxo$primary_lifestyle ~ primary_lifestyle, TRUE ~ "Others")) %>% 
  group_by(primary_lifestyle, Season,Tree_species) %>%
  summarise("reads_count"=sum(counts)) %>% 
  ungroup() %>%
  group_by(Season,Tree_species) %>% mutate("reads_percent"= (reads_count/sum(reads_count))*100 ) 

df$primary_lifestyle[is.na(df$primary_lifestyle)] <- "Unidentified"
#df <- df %>% filter(primary_lifestyle != "Unidentified")

mycolors <- mycolorsFullList %>% filter(Name %in% df$primary_lifestyle)
mycolors <- mycolors %>% arrange(myLevel)
mycolor <- mycolors$color
names(mycolor) <- mycolors$Name

df$primary_lifestyle <- factor(df$primary_lifestyle, levels=mycolors$Name)


ggplot(data=df,aes(y=reads_percent,x=Season, fill=primary_lifestyle))+
  geom_bar(stat="identity")+
  theme_bw()+
  scale_fill_manual(values=mycolor)+
  theme(axis.text.x=element_text(angle = 55,hjust=1#,size=2
  ), axis.title.x = element_blank()
  )+
  
  #facet_wrap(~ Season*Tree_species, scales = "free_x") 
  facet_wrap(~ Tree_species, scales = "free_x") + 
  guides(colour = guide_legend(title.position = "right"))


# 5 * 15 
