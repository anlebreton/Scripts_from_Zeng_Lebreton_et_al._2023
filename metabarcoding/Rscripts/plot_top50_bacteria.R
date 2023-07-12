################ CORE COMMUNITIES ############
library(dplyr)
library(tidyr)


abundance <- read.csv("../1_data/Yulong_bacteria_normalised_abundance_addCols.tsv", sep="\t" )
samples_metadata <- read.table("../0_Metadata/Yulong_sampling_bacteria_metadata.txt")
samples_metadata$Sample <- rownames(samples_metadata)
samples_metadata <- samples_metadata %>% mutate(Sample= gsub("-", ".", Sample)) 

amplicon="bacteria"

abundance$SpeciesID <- paste(abundance$observation_name, abundance$Order, sep="_")
abundance <- abundance %>% mutate(SpeciesID= gsub(" ", "_", SpeciesID)) 


abundance$observation_perc <- (abundance$observation_sum/sum(abundance$observation_sum))*100
number_samples <- length(samples_metadata$Sample)

tmp <-abundance %>% select(c(SpeciesID | starts_with(amplicon)))
tmp2 <- gather(tmp,Sample,observation, -SpeciesID)
df <- tmp2 %>% filter(observation != 0) %>% group_by(SpeciesID) %>%
  summarise(percSample=n()/number_samples*100) 

df <- merge(df, abundance, by="SpeciesID" )
#DF <- df %>% select(observation_name,percSample, observation_perc)

res <-  df %>%  #filter(percSample > 50) #%>%
  arrange(desc(observation_perc)) %>%
   mutate(cumul_observation_perc=cumsum(observation_perc)) %>% 
   filter(cumul_observation_perc <= 25)
 lenValue <- length(res$SpeciesID)
# 
res <-  df %>%  #filter(percSample > 50) #%>%
  arrange(desc(observation_perc)) %>%
  mutate(cumul_observation_perc=cumsum(observation_perc)) %>%
  slice_head(n=lenValue+1)

#tmp <- res %>% select(observation_name, Species, Genus)
#write.table(tmp, file="~/Desktop/Ailaoshan_lithocarpus_Figure_tables/Ailaoshan_12nov/newName.tsv", sep="\t")

library("seqRFLP")
resSequence <- res %>% select(SpeciesID, seed_sequence)
#df.fasta = dataframe2fas(resSequence, file="~/Desktop/bacteria_25perc_order.fa")

rm(number_samples, tmp, tmp2, resSequence)

######## CORE COMMUNITIES : TREE DISPLAY ##########

#scp ITS1_dominant-taxa.fa alebreton@biocomp.nancy.inra.fr:/home/alebreton/Yunnan/top_50perc/
#  module add mafft/7.471 
#module add modeltest-ng/0.1.6
#module add raxml-ng/0.9.0
#mafft --auto ITS1_dominant-taxa.fa > ITS1_dominant-taxa.mafft 
#modeltest-ng -i ITS1_dominant-taxa.mafft 
#raxml-ng --msa ITS1_dominant-taxa.mafft --model TIM2+I+G4 --bs-trees autoMRE{100} -threads 1

library(ggtree)
library(RColorBrewer)
library(ggplot2)

mycolorsFullList <- read.csv("../..//mycolorFile.txt",sep="\t")

#tree <- read.tree("../2_dominant/rooted_ITS2_top50perc.nwk")
tree <- read.tree("../2_dominant/rooted_bacteria_top25perc_order.nwk")
tree$tip.label
res$SpeciesID

setdiff(tree$tip.label, res$SpeciesID)
setdiff( res$SpeciesID, tree$tip.label)

d=fortify(tree)
dd = subset(d, isTip)
tips.order <- dd$label[order(dd$y, decreasing=FALSE)]
rm(d,dd)


tmp <- res %>% 
  #group_by(Order) %>% 
  group_by(Phylum) %>% 
  summarise(i = list(SpeciesID))
#myGroup <- setNames(tmp$i,tmp$Order)
myGroup <- setNames(tmp$i,tmp$Phylum)

rm(tmp)

#ggtree(tree, layout = "circular", branch.length="none")

mycolors <- mycolorsFullList %>% filter(Name %in% res$Phylum)
#mycolors <- mycolorsFullList %>% filter(Name %in% res$Order)
mycolor <- mycolors$color
names(mycolor) <- mycolors$Name

tree2 <- groupOTU(tree, myGroup, 
                  group_name = "Phylum")

# tree2 <- groupOTU(tree, myGroup, 
#                   group_name = "Order")
ggtree(tree)+  geom_tiplab( align=TRUE)

p <- ggtree(tree2,#layout="fan", 
            #layout = "circular",
            branch.length="none",
            aes(color=Phylum)) + 
            #aes(color=Order)) + 
  #geom_tiplab(aes(color=Order), align=TRUE,
              
  geom_tiplab(aes(color=Phylum), align=TRUE,
              linesize=0,
              #as_ylab =(vjust=0), 
              hjust=0,
              #offset =4,
              size=3
  )+
  theme(legend.position = "bottom")+
  #geom_nodelab() +
  #,align=F,size=1) +
  scale_color_manual(values=mycolor)+xlim(NA, 24)
p
#p <- rotate_tree(p, -90)
#p

#600*450 SVG

library(pheatmap)

head(res)

res2 <- res %>% select(SpeciesID, Genus, Phylum, starts_with(amplicon) ) %>% 
                pivot_longer(cols= -c(SpeciesID,Genus,Phylum), names_to ="Sample", values_to = "abundance" ) %>%
                left_join(samples_metadata, by="Sample") #%>%
                #left_join(ecology)

#### Layer ####
plotByLayer <- res2 %>% group_by(Layer,SpeciesID) %>%
  summarise(nb=sum(abundance))%>% ungroup() %>%
  group_by(SpeciesID) %>% mutate(abundance_by_Name=sum(nb)) %>% ungroup() %>%
  mutate(pcA_per_Name_layer=nb/abundance_by_Name*100)

mycolors <- mycolorsFullList %>% filter(Name %in% plotByLayer$Layer)
mycolors <- mycolors %>% arrange(myLevel)
mycolor <- mycolors$color
names(mycolor) <- mycolors$Name
plotByLayer$Layer <- factor(plotByLayer$Layer, levels=mycolors$Name)
plotByLayer$SpeciesID <- factor(plotByLayer$SpeciesID, levels=tips.order)

pL <- ggplot(data=plotByLayer,aes(y=SpeciesID,x=pcA_per_Name_layer, fill=Layer))+
      geom_bar(stat="identity")+
      ggtitle("Layer") +
      theme_bw()+
      scale_fill_manual(values=mycolor)+
      theme(axis.title = element_blank(), legend.position = "bottom")  #axis.text.x=element_text(angle = 55,hjust=1,size=2))#+
      #facet_wrap(~ Season*Tree_species, scales = "free_x") 
      #facet_wrap(~ Tree_species, scales = "free_x") 

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
plotBySeason$Season <- factor(plotBySeason$Season, levels=mycolors$Name)
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

# 
# plotAbundance<- res2 %>% group_by(SpeciesID,primary_lifestyle) %>%
#   summarise(nb=sum(abundance))
# 
# mycolors <- mycolorsFullList %>% filter(Name %in% plotAbundance$primary_lifestyle)
# mycolors <- mycolors %>% arrange(myLevel)
# mycolor <- mycolors$color
# names(mycolor) <- mycolors$Name
# plotAbundance$primary_lifestyle <- factor(plotAbundance$primary_lifestyle, levels=rev(mycolors$Name))

### bacteria
plotAbundance<- res2 %>% group_by(SpeciesID,Phylum) %>%
  summarise(nb=sum(abundance))
mycolors <- mycolorsFullList %>% filter(Name %in% plotAbundance$Phylum)
mycolors <- mycolors %>% arrange(myLevel)
mycolor <- mycolors$color
names(mycolor) <- mycolors$Name
plotAbundance$Phylum <- factor(plotAbundance$Phylum, levels=rev(mycolors$Name))


### bacteria

plotAbundance$SpeciesID <- factor(plotAbundance$SpeciesID, levels=tips.order)

pA <- ggplot(data=plotAbundance,aes(y=SpeciesID,x=nb, #fill=primary_lifestyle
                                    fill=Phylum
                                    ))+
  geom_bar(stat="identity")+
  ggtitle("Overall abundance") +
  theme_bw()+
  scale_fill_manual(values=mycolor)+
  theme(axis.title = element_blank(), legend.position = "bottom") 

pA
pA <- pA+theme( axis.text.y =element_blank(),
                axis.ticks.y=element_blank(), legend.title =element_blank()  ) 



 library(ggpubr)
 ggarrange(p, pT,pL,pS, pA, 
           align="hv",
           ncol = 5, nrow = 1,
           widths =c(3,1,1,1,2))
 







