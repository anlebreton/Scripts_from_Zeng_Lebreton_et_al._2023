


abundance <- read.csv("../1_Abundance_tables/Yulong_ITS1_abundance_addCol.tsv", sep="\t" )
samples_metadata <- read.table("../0_Metadata/Yulong_sampling_ITS1_metadata.txt")
samples_metadata$Sample <- rownames(samples_metadata)
samples_metadata$SampleTag <- paste(samples_metadata$cluster_number, samples_metadata$Year_Season, samples_metadata$layer, sep="_") 

filteredFungalDB2 <-read.csv(file="../filteredFungalDB_Genus_Merged.tsv", sep = "\t")
DE <- read.csv(file="../DE_OTU_LFC5.txt", sep = "\t")
DE <- read.csv(file="../DE_OTU_season_LFC5.txt", sep = "\t")


df <-  abundance %>% filter(observation_name %in% DE$observation_name)
df$ID <- paste(df$observation_name, df$Species, sep="_")

df2 <- df %>% select(c("ID", starts_with("ITS1")))

 restmp2 <-  gather(df2,Sample,observation, -ID)
 restmp <- merge(restmp2, samples_metadata) %>% unique()
 restmp$SampleID <- paste(restmp$cluster_number, restmp$Year_Season, sep="_" )
 restmp2 <- restmp %>% select(ID,SampleID, observation) %>%
   group_by(SampleID, ID) %>%
   summarise(observation_sum=sum(observation)) 
 

 df2 <- spread(restmp2,SampleID,observation_sum)

#df2 <- spread(restmp2,Sample, observation)

m <- df2
m$ID <- NULL
rownames(m) <- df2$ID


df3 <- merge(df,filteredFungalDB2, all.x = T)
df3$guild_fg[is.na(df3$guild_fg)] <- "Unidentified"

ann <- df3 %>% select(guild_fg, Order,Phylum )
rownames(ann) <- df3$ID

topNelement <- df3 %>%  group_by(Order) %>% mutate(nbObs =n() ) %>%
  select(Order, nbObs) %>% unique() %>% ungroup() %>%
  filter(nbObs > 4)
  #slice_max(nbObs, n = 7)

ann$Order[!(ann$Order %in% topNelement$Order) & !(ann$Order == "unidentified")] <- "Other"

#annSample <- samples_metadata %>% select(Year_Season, plot_name, layer)
annSample <- samples_metadata %>% select(Year_Season, plot_name, cluster_number) %>% unique()
rownames(annSample) <-  paste(annSample$cluster_number, annSample$Year_Season, sep="_" )
annSample$cluster_number <- NULL
  
annotation_colors = list("guild_fg"=c(
                                      "Animal Pathogen-Plant Pathogen-Endophyte"="black",
                                      "Ectomycorrhizal"="forestgreen",
                                      "Ectomycorrhizal-Orchid Mycorrhizal-Root Associated Biotroph"="#8cd9b3",
                                      "Ectomycorrhizal-Undefined Saprotroph" = "tan3", 
                                      "Ectomycorrhizal-Wood Saprotroph" = "darkred",
                                      #"Endophyte-Plant Pathogen"="steelblue",
                                      "Endophyte"="yellow",
                                      #"Undefined Saprotroph"="blue"
                                      #"Leaf Saprotroph"  = "darkblue"
                                      #"Other"="grey",
                                      "Unidentified"="white"
                                      ),
                         "Phylum"=c("Basidiomycota"="grey","Ascomycota"="black","Multi-affiliation"="white",
                                         "unidentified"="white"),
                         "Order"=c(#"0"="black", 
                                   "Agaricales"="#8c510a",
                                   "Atheliales"="#dfc27d",
                                   #"Archaeorhizomycetales"="#bf812e", "Chaetothyriales"="#dfc27d",
                                   #"Eurotiales"="darkblue", #"Geminibasidiales"="blue",
                                   "Russulales"="#80cdc1", "Sebacinales"="#35978f",
                                   "Thelephorales"="black", #"Trechisporales"="darkred","Tremellales"="black",
                                   #"Cantharellales"="grey",
                                   #"Pezizales"="grey",
                                   #"Boletales"="grey",
                                   #"Multi-affiliation" ="white",
                                   #"Helotiales"="grey",
                                   #"Sordariales"="grey",   
                                   "Other"="grey",
                                   "unidentified"="white"),
                         "layer"=c("OS"="grey","OM"="beige"),
                         "plot_name"=c("Abies"="#ffb3b3","Picea"="#c6ecd9","Quercus"="#66ccff"),
                         "Year_Season"=c("2019_wet"="#35978f", "2020_dry"="#dfc27d","2020_wet"="navyblue")
                         )





pheatmap(log(m+1),
         colorRampPalette(c("#e6f5ff", "navyblue"), space = "Lab")(100),
         clustering_method = "ward.D2",
         annotation_row = ann,
         annotation_col = annSample,
         annotation_colors = annotation_colors,
         #show_colnames = FALSE,
         fontsize = 8,
angle_col = 45,
cutree_rows = 3,
#cutree_cols =7
#fontsize_row = 3
)


?pheatmap



pheatmap(log(m+1))
