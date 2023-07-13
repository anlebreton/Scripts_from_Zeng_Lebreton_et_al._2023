################ CORE COMMUNITIES  PART 1 ############
library(dplyr)
library(tidyr)
library("seqRFLP")

## output of the BIOM to TSV of FROGS : contingency table with multiple metadata attached (OTU in row, sample in column), 
# fasta sequence and taxonomy metadata are used.
abundance <- read.csv("../1_data/Yulong_ITS2_normalised_abundance_addCols.tsv", sep="\t" )


samples_metadata <- read.table("../0_Metadata/Yulong_sampling_ITS2_metadata.txt")
samples_metadata$Sample <- rownames(samples_metadata)
samples_metadata <- samples_metadata %>% mutate(Sample= gsub("-", ".", Sample)) 

# all samples here have a name starting by the amplicon (column selection purpose) 
amplicon="ITS2"

### ecology obtained on supplementary files of FungalTrait article (https://doi.org/10.1007/s13225-020-00466-2)
ecology <-read.csv(file="../0_Metadata/ecology_fungaltrait_genus.txt", sep = "\t")

abundance$SpeciesID <- paste(abundance$observation_name, abundance$Species, sep="_")
abundance <- abundance %>% mutate(SpeciesID= gsub(" ", "_", SpeciesID)) 


### obtain observation % and number of sample 
abundance$observation_perc <- (abundance$observation_sum/sum(abundance$observation_sum))*100
number_samples <- length(samples_metadata$Sample)

tmp <-abundance %>% select(c(SpeciesID | starts_with(amplicon)))
tmp2 <- gather(tmp,Sample,observation, -SpeciesID)
df <- tmp2 %>% filter(observation != 0) %>% group_by(SpeciesID) %>%
  summarise(percSample=n()/number_samples*100) 

# add these information to the main dataframe 
df <- merge(df, abundance, by="SpeciesID" )

# identify the OTU containing 50% of the abundance 
res <-  df %>%  
  arrange(desc(observation_perc)) %>%
  mutate(cumul_observation_perc=cumsum(observation_perc)) %>%
  filter(cumul_observation_perc <= 50) # 25 for bacteria, 50 for fungi
lenValue <- length(res$SpeciesID)

res <-  df %>%  
  arrange(desc(observation_perc)) %>%
  mutate(cumul_observation_perc=cumsum(observation_perc)) %>%
  slice_head(n=lenValue+1)


# get the sequence in fasta of the top OTU
resSequence <- res %>% select(SpeciesID, seed_sequence)
#write it down to perform the phylogenetic tree outside of R
#write.table(resSequence,file="ITS1_50perc-Abundance.tsv", row.names = F, sep="\t")
#df.fasta = dataframe2fas(resSequence, file="~/Desktop/bacteria_dominant-taxa_order.fa")

rm(number_samples, resSequence) 

######## CORE COMMUNITIES : TREE DISPLAY ##########

# outside R on a cluster...
# mafft/7.471 
# modeltest-ng/0.1.6
# raxml-ng/0.9.0
# mafft --auto ITS1_dominant-taxa.fa > ITS1_dominant-taxa.mafft 
# modeltest-ng -i ITS1_dominant-taxa.mafft 
# raxml-ng --msa ITS1_dominant-taxa.mafft --model < model test res > --bs-trees autoMRE{100} -threads 1
# root the tree on R or outside (e.g. in itol  https://itol.embl.de )

library(ggtree)
library(RColorBrewer)
library(ggplot2)
library(pheatmap)
library(ggpubr)

mycolorsFullList <- read.csv("../../../mycolorFile.txt",sep="\t")

#import tree
tree <- read.tree("/Volumes/Hungry/Yulong/metabarcoding/2_dominant/rooted_ITS2_top50perc.nwk")
tree$tip.label
res$SpeciesID

# check if the ID of tree and ID of dataframe are the same 
setdiff(tree$tip.label, res$SpeciesID)
setdiff( res$SpeciesID, tree$tip.label)

# get the order of OTU ID in the tree
d=fortify(tree)
dd = subset(d, isTip)
tips.order <- dd$label[order(dd$y, decreasing=FALSE)]
rm(d,dd)

# make group for the color of the branch of the tree
tmp <- res %>% 
  group_by(Order) %>% 
  summarise(i = list(SpeciesID))
myGroup <- setNames(tmp$i,tmp$Order)

rm(tmp)

# get colors
mycolors <- mycolorsFullList %>% filter(Name %in% res$Order)
mycolor <- mycolors$color
names(mycolor) <- mycolors$Name


# integrate the grouping in the tree
 tree2 <- groupOTU(tree, myGroup, 
                   group_name = "Order")
ggtree(tree)+  geom_tiplab()

# plot tree 
p <- ggtree(tree2,
            branch.length="none",
            aes(color=Order)) + 
  geom_tiplab(aes(color=Order),
              linesize=0,
              hjust=0,
              size=3
  )+
  theme(legend.position = "bottom")+
  scale_color_manual(values=mycolor)+xlim(NA, 24)
p

#600*450 SVG

#long version of the abundance with all metadata
res2 <- res %>% select(SpeciesID, Genus, Phylum, starts_with(amplicon) ) %>% 
                pivot_longer(cols= -c(SpeciesID,Genus,Phylum), names_to ="Sample", values_to = "abundance" ) %>%
                left_join(samples_metadata, by="Sample") #%>%


#### Layer ####
plotByLayer <- res2 %>% group_by(Layer,SpeciesID) %>%
  summarise(nb=sum(abundance))%>% ungroup() %>%
  group_by(SpeciesID) %>% mutate(abundance_by_Name=sum(nb)) %>% ungroup() %>%
  mutate(pcA_per_Name_layer=nb/abundance_by_Name*100)

mycolors <- mycolorsFullList %>% filter(Name %in% plotByLayer$Layer)
mycolors <- mycolors %>% arrange(myLevel)
mycolor <- mycolors$color
names(mycolor) <- mycolors$Name
plotByLayer$Layer <- factor(plotByLayer$Layer, levels=rev(mycolors$Name))
plotByLayer$SpeciesID <- factor(plotByLayer$SpeciesID, levels=tips.order)

pL <- ggplot(data=plotByLayer,aes(y=SpeciesID,x=pcA_per_Name_layer, fill=Layer))+
      geom_bar(stat="identity")+
      ggtitle("Layer") +
      theme_bw()+
      scale_fill_manual(values=mycolor)+
      theme(axis.title = element_blank(), legend.position = "bottom") 

pL
pL <- pL+theme( axis.text.y=element_blank(),
                axis.ticks.y=element_blank(), legend.title =element_blank() ) 

### Season ###
plotBySeason<- res2 %>% group_by(Season,SpeciesID) %>%
  summarise(nb=sum(abundance))%>% ungroup() %>%
  group_by(SpeciesID) %>% mutate(abundance_by_Name=sum(nb)) %>% ungroup() %>%
  mutate(pcA_per_Name_Season=nb/abundance_by_Name*100)

mycolors <- mycolorsFullList %>% filter(Name %in% plotBySeason$Season)
mycolors <- mycolors %>% arrange(myLevel)
mycolor <- mycolors$color
names(mycolor) <- mycolors$Name
plotBySeason$Season <- factor(plotBySeason$Season, levels=rev(mycolors$Name))
plotBySeason$SpeciesID <- factor(plotBySeason$SpeciesID, levels=tips.order)

pS <- ggplot(data=plotBySeason,aes(y=SpeciesID,x=pcA_per_Name_Season, fill=Season))+
      geom_bar(stat="identity")+
      ggtitle("Season") +
      theme_bw()+
      scale_fill_manual(values=mycolor)+
      theme(axis.title = element_blank(), legend.position = "bottom") 
pS
pS <- pS+theme( axis.text.y=element_blank(),
                axis.ticks.y=element_blank(), legend.title =element_blank() ) 

### Tree_species ###
plotByTree_species<- res2 %>% group_by(Tree_species,SpeciesID) %>%
  summarise(nb=sum(abundance))%>% ungroup() %>%
  group_by(SpeciesID) %>% mutate(abundance_by_Name=sum(nb)) %>% ungroup() %>%
  mutate(pcA_per_Name_Tree_species=nb/abundance_by_Name*100)

mycolors <- mycolorsFullList %>% filter(Name %in% plotByTree_species$Tree_species)
mycolors <- mycolors %>% arrange(myLevel)
mycolor <- mycolors$color
names(mycolor) <- mycolors$Name
mycolor = c("Abies"="#d64933", "Picea"="#91f086", "Quercus"="#5a92ed" )
plotByTree_species$Tree_species <- factor(plotByTree_species$Tree_species, levels=rev(mycolors$Name))
plotByTree_species$SpeciesID <- factor(plotByTree_species$SpeciesID, levels=tips.order)

pT <- ggplot(data=plotByTree_species,aes(y=SpeciesID,x=pcA_per_Name_Tree_species, fill=Tree_species))+
      geom_bar(stat="identity")+
      ggtitle("Tree species") +
      theme_bw()+
      scale_fill_manual(values=mycolor)+
      theme(axis.title = element_blank(), legend.position = "bottom") 

pT
pT <- pT+theme( axis.text.y =element_blank(),
                axis.ticks.y=element_blank(), legend.title =element_blank()  ) 

###### Abundance panel ###
res2 <- res2 %>% left_join(ecology)
plotAbundance<- res2 %>% group_by(SpeciesID,primary_lifestyle) %>%
  summarise(nb=sum(abundance)/1000)

mycolors <- mycolorsFullList %>% filter(Name %in% plotAbundance$primary_lifestyle)
mycolors <- mycolors %>% arrange(myLevel)
mycolor <- mycolors$color
names(mycolor) <- mycolors$Name
plotAbundance$primary_lifestyle <- factor(plotAbundance$primary_lifestyle, levels=rev(mycolors$Name))


plotAbundance$SpeciesID <- factor(plotAbundance$SpeciesID, levels=tips.order)

pA <- ggplot(data=plotAbundance,aes(y=SpeciesID,x=nb, fill=primary_lifestyle
                                    #fill=Oligo_copio_final
                                    ))+
  geom_bar(stat="identity")+
  ggtitle("Overall abundance (thousand reads)") +
  theme_bw()+
  scale_fill_manual(values=mycolor)+
  theme(axis.title = element_blank(), legend.position = "bottom") 

pA
pA <- pA+theme( axis.text.y =element_blank(),
                axis.ticks.y=element_blank(), legend.title =element_blank()  ) 


### multi pannel figure 

 ggarrange(p, pT,pL,pS, pA, 
           align="hv",
           ncol = 5, nrow = 1,
           widths =c(3,1,1,1,2))
 







