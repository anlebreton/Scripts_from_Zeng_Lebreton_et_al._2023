################ CORE COMMUNITIES ############
library(dplyr)
library(tidyr)

#abundance <- read.csv("../1_data/1_Abundance_tables/Yulong_ITS1_abundance_addCol.tsv", sep="\t" )
abundance <- read.csv("../1_data/Yulong_ITS1_normalised_abundance_addCols.tsv", sep="\t" )

samples_metadata <- read.table("../0_Metadata/Yulong_sampling_ITS1_metadata.txt")
samples_metadata$Sample <- rownames(samples_metadata)
filteredFungalDB2 <-read.csv(file="../filteredFungalDB_Genus_Merged.tsv", sep = "\t")

abundance$SpeciesID <- paste(abundance$observation_name, abundance$Species, sep="_")
abundance$observation_perc <- (abundance$observation_sum/sum(abundance$observation_sum))*100
number_samples <- length(samples_metadata$Sample)


tmp <-abundance %>% select(c(SpeciesID | starts_with("ITS1")))
tmp2 <- gather(tmp,Sample,observation, -SpeciesID)
df <- tmp2 %>% filter(observation != 0) %>% group_by(SpeciesID) %>%
  summarise(percSample=n()/number_samples*100) 

df <- merge(df, abundance, by="SpeciesID" )
#DF <- df %>% select(observation_name,percSample, observation_perc)

res <-  df %>%  #filter(percSample > 50) #%>%
  filter(observation_perc > 0.59)

sum(res$observation_perc)

#tmp <- res %>% select(observation_name, Species, Genus)
#write.table(tmp, file="~/Desktop/Ailaoshan_lithocarpus_Figure_tables/Ailaoshan_12nov/newName.tsv", sep="\t")

resSequence <- res %>% select(SpeciesID, seed_sequence)
#write.table(resSequence, file="../1_Abundance_tables/ITS1_dominant-taxa.tsv", row.names = F)
#write.table(res, file="~/Desktop/Ailaoshan_lithocarpus_Figure_tables/ITS1_dominant-taxa_v2.tsv", row.names = F, sep="\t")

rm(number_samples, tmp, tmp2, resSequence)

######## CORE COMMUNITIES : TREE DISPLAY ##########

#edits in vi to switch to fasta
#mafft --auto ITS1_dominant-taxa.fa > ITS1_dominant-taxa.mafft 
#modeltest-ng -i ITS1_dominant-taxa.mafft 
#raxml-ng --msa ITS1_dominant-taxa.mafft --model TIM2+I+G4 --bs-trees autoMRE{100} -threads 1


library(ggtree)
library(RColorBrewer)

tree <- read.tree("..//2_dominant/ITS1_dominant-taxa.mafft.raxml.bestTree")
tree$tip.label
res$SpeciesID

tmp <- res %>% 
  group_by(Order) %>% 
  summarise(i = list(SpeciesID))
myGroup <- setNames(tmp$i,tmp$Order)
rm(tmp)

#ggtree(tree, layout = "circular", branch.length="none")
mycolor=c(rev(brewer.pal(7,"BrBG")))
mycolor=c("0"="black", "Agaricales"="#8c510a","Archaeorhizomycetales"="#bf812e","Atheliales"="#dfc27d", #"Chaetothyriales"="#dfc27d",
          "Eurotiales"="darkblue", #"Geminibasidiales"="blue",
          "Russulales"="#80cdc1", "Sebacinales"="#35978f","Cantharellales"="blue", "Eurotiales"="darkblue",
          "Hypocreales" ="red"
          
         # "Thelephorales"="red", "Trechisporales"="darkred","Tremellales"="black"
         )

tree2 <- groupOTU(tree, myGroup, 
                  group_name = "Order")

p <- ggtree(tree2,#layout="fan", 
            #layout = "circular",
            branch.length="none",
            aes(color=Order)) + 
  geom_tiplab(aes(color=Order), align=TRUE,
              linesize=0,
              #as_ylab =(vjust=0), 
              hjust=0,
              #offset =4,
              size=3
  )+
  #geom_nodelab() +
  #,align=F,size=1) +
  scale_color_manual(values=mycolor) #+
#xlim(NA, 24)
p
#p <- rotate_tree(p, -90)
#p



################ TEST #############################
library(ggtreeExtra)
library(ggtree)
library(treeio)
library(tidytree)
library(ggstar)
library(ggplot2)
library(ggnewscale)

res2 <- res %>% select(starts_with("ITS1"),SpeciesID )
tmp <-res2 %>%  gather(Sample, abundance, -SpeciesID)
res2 <- merge(tmp, samples_metadata)

order_newname <- res %>% select(SpeciesID, Genus)
order_newname <- merge(order_newname, filteredFungalDB2, all.x = T)

res2 <- merge(res2,order_newname, all.x = T)


dat1 <- res2 %>% select(SpeciesID,guild_fg) %>% unique()
rn <- dat1$SpeciesID
dat1$SpeciesID <- NULL
rownames(dat1) <- rn
dat1$guild_fg[is.na(dat1$guild_fg)] <- "Unidentified"


colorEco=c("Ectomycorrhizal"="forestGreen",
           "Ectomycorrhizal-Orchid Mycorrhizal-Root Associated Biotroph"="#8cd9b3",
           "Ectomycorrhizal-Undefined Saprotroph"="blue",
           "Endophyte-Plant Pathogen"="brown",
           "Unidentified"="white" #"Endophyte"="brown",
           )

p2 <- gheatmap(p, dat1, offset = 9, #colnames_offset_y = 0.2, font.size=2,
               width=.1, color=NULL, 
               colnames_position="top", 
               colnames_angle=90,  
               hjust=0)+
  scale_fill_manual(values=c(colorEco ))

p2

plotByLayer <- res2 %>% group_by(layer,SpeciesID) %>%
  summarise(nb=sum(abundance))


plotByLayer$layer <- factor(plotByLayer$layer, levels = c("HF","OS","OM"))


plotBySeason <- res2 %>% group_by(Year_Season,SpeciesID) %>%
  summarise(nb=sum(abundance)) %>% ungroup() %>%
  group_by(SpeciesID) %>% mutate(abundance_by_Name=sum(nb)) %>% ungroup() %>%
  mutate(pcA_per_Name=nb/abundance_by_Name*100)

plotByLayer <- res2 %>% group_by(layer,SpeciesID) %>%
  summarise(nb=sum(abundance))%>% ungroup() %>%
  group_by(SpeciesID) %>% mutate(abundance_by_Name=sum(nb)) %>% ungroup() %>%
  mutate(pcA_per_Name_layer=nb/abundance_by_Name*100)

plotAbundance <- res2 %>% group_by(SpeciesID, plot_name) %>%
  summarise(nb=sum(abundance)) %>% ungroup() %>%
  group_by(SpeciesID) %>% mutate(abundance_by_Name=sum(nb)) %>% ungroup() %>%
  mutate(pcA_per_Name_tree=nb/abundance_by_Name*100)




# p <- ggtree(tree, layout = "circular", branch.length="none") +
#   geom_tiplab(align=TRUE,
#               as_ylab =(vjust=0),
#               hjust=1,
#               offset = 30)
# p


colorLayer=c( "OS"="grey","OM"="beige")
colorSeason=c("2019_wet"="#35978f", "2020_dry"="#dfc27d","2020_wet"="navyblue")
treecolor=c("Abies"="#ffb3b3","Picea"="#c6ecd9","Quercus"="#66ccff")

p2  + new_scale_fill() +
  geom_fruit(data=plotBySeason, geom=geom_tile,
             mapping=aes(y=SpeciesID, x=Year_Season, fill=pcA_per_Name),
             colour = "grey50",
             offset = 0.58,
             size = 0.2)+
  scale_fill_gradient(low="white", high = "red") +
  new_scale_fill() +
  geom_fruit(data=plotByLayer, geom=geom_tile,
             mapping=aes(y=SpeciesID, x=layer, fill=pcA_per_Name_layer),
             colour = "grey50",
             offset = 0.075, #,size = 0.05
             pwidth = 0.1
             )+
  scale_fill_gradient(low="white", high = "blue") +
  new_scale_fill() +
  geom_fruit(
    data=plotAbundance,
    geom=geom_tile,
    mapping = aes(y=SpeciesID,x=plot_name, fill=pcA_per_Name_tree),
    colour = "grey50",
    pwidth=0.2, 
    offset =0.075 ,size = 0.2 ) +
  scale_fill_gradient(low="white", high = "#dfc27d") +
  new_scale_fill() +
  geom_fruit(
    data=plotAbundance,
    geom=geom_bar,
    mapping = aes(y=SpeciesID,x=nb, fill=plot_name),
    #colour="grey50",
    pwidth=0.8, 
    orientation="y", 
    stat="identity",
    offset =0.075  ) +
  scale_fill_discrete(treecolor)
#scale_fill_gradient(low="white", high = "blue") 








