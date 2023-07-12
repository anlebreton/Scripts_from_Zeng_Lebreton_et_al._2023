################ CORE COMMUNITIES ############
library(dplyr)
library(tidyr)

abundance <- read.csv("../1_data/Yulong_ITS2_normalised_abundance_addCols.tsv", sep="\t" )
samples_metadata <- read.table("../0_Metadata/Yulong_sampling_ITS2_metadata.txt")
samples_metadata$Sample <- rownames(samples_metadata)
samples_metadata <- samples_metadata %>% mutate(Sample= gsub("-", ".", Sample)) 

ecology <-read.csv(file="../../ecology_fungaltrait_genus.txt", sep = "\t")

abundance$SpeciesID <- paste(abundance$observation_name, abundance$Species, sep="_")
abundance$observation_perc <- (abundance$observation_sum/sum(abundance$observation_sum))*100
number_samples <- length(samples_metadata$Sample)


tmp <-abundance %>% select(c(SpeciesID | starts_with("ITS2")))
tmp2 <- gather(tmp,Sample,observation, -SpeciesID)
df <- tmp2 %>% filter(observation != 0) %>% group_by(SpeciesID) %>%
  summarise(percSample=n()/number_samples*100) 

df <- merge(df, abundance, by="SpeciesID" )
#DF <- df %>% select(observation_name,percSample, observation_perc)

res <-  df %>%  #filter(percSample > 50) #%>%
  arrange(desc(observation_perc)) %>%
  mutate(cumul_observation_perc=cumsum(observation_perc)) %>%
  filter(cumul_observation_perc <= 50)
lenValue <- length(res$SpeciesID)

res <-  df %>%  #filter(percSample > 50) #%>%
  arrange(desc(observation_perc)) %>%
  mutate(cumul_observation_perc=cumsum(observation_perc)) %>%
  slice_head(n=lenValue+1)

#tmp <- res %>% select(observation_name, Species, Genus)
#write.table(tmp, file="~/Desktop/Ailaoshan_lithocarpus_Figure_tables/Ailaoshan_12nov/newName.tsv", sep="\t")

library("seqRFLP")
resSequence <- res %>% select(SpeciesID, seed_sequence)
#df.fasta = dataframe2fas(resSequence, file="~/Desktop/ITS1_dominant-taxa.fa")

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

mycolorsFullList <- read.csv("../..//mycolorFile.txt",sep="\t")

tree <- read.tree("../2_dominant/rooted_ITS2_top50perc.nwk")
tree$tip.label
res$SpeciesID

tmp <- res %>% 
  group_by(Order) %>% 
  summarise(i = list(SpeciesID))
myGroup <- setNames(tmp$i,tmp$Order)
rm(tmp)

#ggtree(tree, layout = "circular", branch.length="none")

mycolors <- mycolorsFullList %>% filter(Name %in% res$Order)
mycolor <- mycolors$color
names(mycolor) <- mycolors$Name


tree2 <- groupOTU(tree, myGroup, 
                  group_name = "Order")
ggtree(tree)+  geom_tiplab( )

p <- ggtree(tree2,#layout="fan", 
            #layout = "circular",
            branch.length="none",
            aes(color=Order)) + 
  geom_tiplab(aes(color=Order),
              linesize=0,
              #as_ylab =(vjust=0), 
              hjust=0,
              #offset =4,
              size=3
  )+
  #geom_nodelab() +
  #,align=F,size=1) +
  scale_color_manual(values=mycolor)+xlim(NA, 24)
p
#p <- rotate_tree(p, -90)
#p

#600*450 SVG

################ TEST #############################
library(ggtreeExtra)
library(ggtree)
library(treeio)
library(tidytree)
library(ggstar)
library(ggplot2)
library(ggnewscale)

res2 <- res %>% select(SpeciesID, starts_with("ITS2") )
tmp <-res2 %>%  gather(Sample, abundance, -SpeciesID)
res2 <- merge(tmp, samples_metadata, by="Sample")

order_newname <- res %>% select(SpeciesID, Genus)
order_newname <- merge(order_newname, ecology, all.x = T)

res2 <- merge(res2,order_newname, all.x = T)


dat1 <- res2 %>% select(SpeciesID,primary_lifestyle) %>% unique()
rn <- dat1$SpeciesID
dat1$SpeciesID <- NULL
rownames(dat1) <- rn
dat1$primary_lifestyle[is.na(dat1$primary_lifestyle)] <- "Unidentified"

mycolors<-  mycolorsFullList %>% filter(Name %in% dat1$primary_lifestyle)
colorEco <- mycolors$color
names(colorEco) <- mycolors$Name


p2 <- gheatmap(p, dat1, offset = 7.2, #colnames_offset_y = 0.2, font.size=2,
               width=.1, color=NULL, 
               colnames_position="top", 
               colnames_angle=90,  
               hjust=0)+
  scale_fill_manual(values=c(colorEco ))

p2

library(pheatmap)
pheatmap(dat1)

ggplot(data=df,aes(y=counts,x=samples, fill=Phylum))+
  geom_bar(stat="identity")+
  theme_bw()+
  scale_fill_manual(values=mycolor)+
  theme(axis.text.x=element_text(angle = 55,hjust=1,size=2))+
  #facet_wrap(~ Season*Tree_species, scales = "free_x") 
  facet_wrap(~ Tree_species, scales = "free_x") 




plotBySeason <- res2 %>% group_by(Season,SpeciesID) %>%
  summarise(nb=sum(abundance)) %>% ungroup() %>%
  group_by(SpeciesID) %>% mutate(abundance_by_Name=sum(nb)) %>% ungroup() %>%
  mutate(pcA_per_Name=nb/abundance_by_Name*100)

t <- plotBySeason %>% select(Season,SpeciesID, pcA_per_Name) %>% pivot_wider(names_from = Season, values_from = pcA_per_Name)
t2 <- t
t$SpeciesID <- NULL
rownames(t) <- t2$SpeciesID
library(pheatmap)

d1 <- data.frame(id=plotBySeason$SpeciesID, season=plotBySeason$pcA_per_Name)
p1 <- p %<+% d1 + geom_bar(aes(color="red"))
p1

p2 <- pheatmap(t, cluster_rows = F, cluster_cols = F)
p2

plotByLayer <- res2 %>% group_by(Layer,SpeciesID) %>%
  summarise(nb=sum(abundance))%>% ungroup() %>%
  group_by(SpeciesID) %>% mutate(abundance_by_Name=sum(nb)) %>% ungroup() %>%
  mutate(pcA_per_Name_layer=nb/abundance_by_Name*100)

plotAbundance <- res2 %>% group_by(SpeciesID, Tree_species) %>%
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
             mapping=aes(y=SpeciesID, x=Season, fill=pcA_per_Name),
             colour = "grey50",
             offset = 0.58,
             size = 0.2)+
  scale_fill_gradient(low="white", high = "red") +
  new_scale_fill() +
  geom_fruit(data=plotByLayer, geom=geom_tile,
             mapping=aes(y=SpeciesID, x=Layer, fill=pcA_per_Name_layer),
             colour = "grey50",
             offset = 0.075, #,size = 0.05
             pwidth = 0.1
             )+
  scale_fill_gradient(low="white", high = "blue") +
  new_scale_fill() +
  geom_fruit(
    data=plotAbundance,
    geom=geom_tile,
    mapping = aes(y=SpeciesID,x=Tree_species, fill=pcA_per_Name_tree),
    colour = "grey50",
    pwidth=0.2, 
    offset =0.075 ,size = 0.2 ) +
  scale_fill_gradient(low="white", high = "#dfc27d") +
  new_scale_fill() +
  geom_fruit(
    data=plotAbundance,
    geom=geom_bar,
    mapping = aes(y=SpeciesID,x=nb, fill=Tree_species),
    #colour="grey50",
    pwidth=0.8, 
    orientation="y", 
    stat="identity",
    offset =0.075  ) +
  scale_fill_discrete(treecolor)
#scale_fill_gradient(low="white", high = "blue") 

library(ggpubr)



pSeason <- geom_fruit(data=plotBySeason, geom=geom_tile,
             mapping=aes(y=SpeciesID, x=Season, fill=pcA_per_Name),
             colour = "grey50",
             offset = 0.58,
             size = 0.2)


ggarrange(p, p2, # p4, #+ rremove("x.text"), 
          #labels = c("A", "B", "C","D"),
          #ncol = 2, nrow = 2,
          ncol = 2, nrow = 1,
          common.legend = FALSE, legend="bottom") 



library(complex)





