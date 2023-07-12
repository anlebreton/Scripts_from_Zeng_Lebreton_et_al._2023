#install.packages("devtools")
#devtools::install_github("ropenscilabs/datastorr")
#devtools::install_github("traitecoevo/fungaltraits")
library(fungaltraits)
library(dplyr)
library(tidyr)
library(ggplot2)

fungal_traits<- fungal_traits()
fg <- fungal_traits %>% 
  select(species, guild_fg, trophic_mode_fg) %>%
  filter(!is.na(trophic_mode_fg)) %>%
  unique()
  

abundance <- read.csv("../1_Abundance_tables/Yulong_ITS1_abundance_addCol.tsv", sep="\t" )
samples_metadata <- read.table("../0_Metadata/Yulong_sampling_ITS1_metadata.txt")
filteredFungalDB2 <-read.csv(file="../filteredFungalDB_Genus_Merged.tsv", sep = "\t")


samples_metadata$Sample <- rownames(samples_metadata)

abundance$observation_perc <- (abundance$observation_sum/sum(abundance$observation_sum))*100

number_samples <- length(samples_metadata$Sample)

abundance <- merge(filteredFungalDB2, abundance, by="Genus", all.y=T)

#abundance <- merge(fg, abundance, by.x="species", by.y="Species", all.y=T)
fg2 <- abundance %>% select(observation_name, guild_fg, trophic_mode_fg)


tmp <-abundance %>% select(c(observation_name | starts_with("ITS1")))
tmp2 <- gather(tmp,Sample,observation, -observation_name)
tmp2 <- tmp2 %>% filter(observation != 0)
tmp3 <- merge(tmp2,samples_metadata)
tmp3 <- merge(tmp3, fg2, all.x=T)


df <- tmp2 %>% filter(observation != 0) %>% group_by(observation_name) %>%
  summarise(percSample=n()/number_samples*100) 


colnames(tmp4)
tmp4 <- tmp3 %>% filter(!is.na(guild_fg))

ggplot(tmp3 ,aes(x=plot_name , y=observation, fill=guild_fg)) +
  geom_bar(stat="identity")+
  #scale_fill_manual(values=myColor)#+
theme(axis.text.x = element_text(angle = 45,hjust=1)) +
facet_grid(layer ~ Year_Season)

