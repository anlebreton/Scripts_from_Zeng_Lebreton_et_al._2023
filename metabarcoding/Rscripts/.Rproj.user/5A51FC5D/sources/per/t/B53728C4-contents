#library("readr")
library("plyr");library("dplyr");library("ggplot2");library("tidyr");library("tidyverse")
library("RColorBrewer");library("ggpubr");library("vegan")


abundance <- read.csv("../1_data/Yulong_bacteria_normalised_abundance_addCols.tsv", sep="\t" )
abundance$ID <- paste(abundance$observation_name, abundance$Species, sep="_")
#abundance$ID <- paste(abundance$observation_name, abundance$Order, sep="_")

abundance$ID <- str_replace_all(abundance$ID , "-", ".")
abundance$ID <- str_replace_all(abundance$ID , " ", ".")
abundance$ID <- str_replace_all(abundance$ID , "[(]", ".")
abundance$ID <- str_replace_all(abundance$ID , "[)]", ".")


samples_metadata <- read.table("../0_Metadata/Yulong_sampling_bacteria_metadata.txt")
samples_metadata$Sample <- rownames(samples_metadata)
samples_metadata$Sample <- str_replace_all(samples_metadata$Sample , "-", ".")

ecology <-read.csv(file="../../ecology_fungaltrait_genus.txt", sep = "\t")
amplicon="bacteria"

df <-  abundance #%>% left_join(ecology) %>% filter(grepl(pattern = "sapro", primary_lifestyle)) 

#%>% filter(observation_sum >50)
#df <-  abundance %>% filter(observation_name %in% DE$observation_name)
#df$ID <- paste(df$observation_name, df$Species, sep="_")
#df$ID <- str_replace_all(df$ID , "-", ".")

df <-  abundance %>% select(c("ID", starts_with(amplicon))) %>% 
pivot_longer(cols=-c("ID"), names_to="Sample", values_to="count") %>% 
  left_join(samples_metadata ) %>%
  group_by(ID, Cluster) %>% 
  summarise(median=median(count),.groups = 'keep') %>% 
  pivot_wider(id_cols = ID, names_from =Cluster, values_from =median  ) #%>%
  #filter_all(all_vars(. != 0))

df2 <- df
tmp <- df
df2$ID <- NULL
rownames(df2) <- tmp$ID
nb <- length(colnames(df2))
df2[df2 != 0] <- 1

df2$NbCluster <- rowSums(df2)

summary(df2$NbCluster)
t <-df2 %>% filter(NbCluster >7)

ggplot(df2, aes(x = NbCluster))+
  geom_bar()

ggplot(df2, aes(x = NbCluster)) +
  #fill=JGI_ID)) +
  #geom_histogram(binwidth=2)+
  theme_light()+
  #theme(#strip.text.y = element_text(angle = 0),
  #      strip.text.x = element_text(angle = 90) )+
  theme(legend.position="bottom") +
  guides( color = FALSE) +
  geom_freqpoly(binwidth=1) 





data <-  df %>% select(c("ID", starts_with(amplicon)))
data[data==0] <- NA
data2<-data[complete.cases(data),]

df2 <- df %>% select(c("ID", starts_with(amplicon))) 
tmp <- df2
df2$ID <- NULL
rownames(df2) <- tmp$ID
nb <- length(colnames(df2))
df2[df2 != 0] <- 1

df2$Sum <- (rowSums(df2)/nb)*100

ggplot(df2, aes(x = Sum)) +
  #fill=JGI_ID)) +
  #geom_histogram(binwidth=2)+
  theme_light()+
  #theme(#strip.text.y = element_text(angle = 0),
  #      strip.text.x = element_text(angle = 90) )+
  theme(legend.position="bottom") +
  guides( color = FALSE) +
  geom_freqpoly(binwidth=2) #+ 
  facet_grid(Lifestyle ~ LTR_ID )



df2$ID <- rownames(df2)

tmp <- df2 %>% pivot_longer(cols=-c("ID"), names_to="samples", values_to="count") %>%
group_by(ID,samples) %>% summarise("counts"=sum(count)) 

tmp <- df2 %>% pivot_longer(cols=-c("ID"), names_to="ID", values_to="count") %>%
  group_by(ID,samples) %>% summarise("counts"=sum(count)) 


data$new <- rowSums(data[43:167])



df2 <- as.data.frame(t(df2))
df2$Sample <- rownames(df2)
df2 <- df2 %>% left_join(samples_metadata )

