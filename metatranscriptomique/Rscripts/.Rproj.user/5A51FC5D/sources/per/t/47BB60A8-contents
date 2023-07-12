library("readr")
library(rstatix)
library("plyr");library("dplyr");library("ggplot2");library("tidyr");library("tidyverse")
library("RColorBrewer");library("ggpubr");library("vegan")

"%notin%" <- Negate("%in%")


GetMyColor <- function(colDF){
  mycolorsCAZymes <- read.csv("../../mycolorFile.txt",sep="\t")
  mycolors <- mycolorsCAZymes %>% filter(Name %in% colDF)
  mycolors <- mycolors %>% arrange(myLevel)
  mycolor <- mycolors$color
  names(mycolor) <- mycolors$Name
  #colDF <- factor(colDF, levels=mycolors$Name
  return(mycolor)
}

GetMyColorCAZymes <- function(colDF){
  mycolorsCAZymes <- read.csv("../0_metadata/mycolorCAZymes_v2.txt",sep="\t")
  mycolors <- mycolorsCAZymes %>% filter(Name %in% colDF)
  mycolors <- mycolors %>% arrange(myLevel)
  mycolor <- mycolors$color
  names(mycolor) <- mycolors$Name
  #colDF <- factor(colDF, levels=mycolors$Name
  return(mycolor)
}

generate_label_df <- function(tukey.res, variable){
  # Extract labels and factor levels from Tukey post-hoc 
  Tukey.levels <- tukey.res[[variable]][,4]
  Tukey.labels <- data.frame(multcompLetters(Tukey.levels)['Letters'])
  
  #I need to put the labels in the same order as in the boxplot :
  Tukey.labels$treatment=rownames(Tukey.labels)
  Tukey.labels=Tukey.labels[order(Tukey.labels$treatment) , ]
  return(Tukey.labels)
}
# 
# Fungi_ratio <- Fungi_TPM %>%
#   dplyr::select(-starts_with("NR_"),-starts_with("Myc"),-ID,-ctg_length, -sum_contig_count) %>% 
#   group_by(tax_phylum) %>% summarise_if(is.numeric, sum) %>% pivot_longer(cols=-c(tax_phylum), names_to="samples", values_to="tpm_counts") %>%
#   left_join(treatment_metadata) #%>%




############# TEST OS-OM ##################
### t test/wilcox paired #####

# 
# WhereToTest <- "tax_phylum"
# SelectMetadataValue="all"
# columnMetadataValue="all"

get.Fungi_ratio <- function(Fungi_TPM, WhereToTest, SelectMetadataValue, columnMetadataValue ){
  Fungi_ratio <- Fungi_TPM %>%
    dplyr::select(-starts_with("NR_"),-starts_with("Myc"),-ID,-ctg_length, -sum_contig_count) %>% 
    group_by(get(WhereToTest)) %>% summarise_if(is.numeric, sum) %>% pivot_longer(cols=-c("get(WhereToTest)"), names_to="samples", values_to="tpm_counts") %>%
    left_join(treatment_metadata)  %>% 
    filter(!grepl("score", samples))%>% 
    filter(!grepl("evalue", samples))%>% 
    filter(!grepl("coverage", samples))
  colnames(Fungi_ratio)[colnames(Fungi_ratio) == "get(WhereToTest)"] <- WhereToTest
  
  if(SelectMetadataValue != "all"){
    Fungi_ratio <- Fungi_ratio %>% filter( get(columnMetadataValue) == SelectMetadataValue )
  }
  
  return(Fungi_ratio)
}
 
# Fungi_ratio <- get.Fungi_ratio(Fungi_TPM, WhereToTest,SelectMetadataValue="all", columnMetadataValue="all")
# Fungi_ratio <- get.Fungi_ratio(Fungi_TPM, WhereToTest,SelectMetadataValue="Quercus", columnMetadataValue="tree_species")

get.Fungi_ratio.2groupments <- function(Fungi_TPM, WhereToTest,WhereToTest2, SelectMetadataValue, columnMetadataValue ){
  Fungi_ratio <- Fungi_TPM %>%
    dplyr::select(-starts_with("NR_"),-starts_with("Myc"),-ID,-ctg_length, -sum_contig_count) %>% 
    group_by(get(WhereToTest), get(WhereToTest2)) %>% summarise_if(is.numeric, sum) %>% pivot_longer(cols=-c("get(WhereToTest)", "get(WhereToTest2)"), names_to="samples", values_to="tpm_counts") %>%
    left_join(treatment_metadata) %>% 
    filter(!grepl("score", samples))%>% 
    filter(!grepl("evalue", samples))%>% 
    filter(!grepl("coverage", samples))
    
    
    
  
  colnames(Fungi_ratio)[colnames(Fungi_ratio) == "get(WhereToTest)"] <- WhereToTest
  colnames(Fungi_ratio)[colnames(Fungi_ratio) == "get(WhereToTest2)"] <- WhereToTest2
  
  if(SelectMetadataValue != "all"){
    Fungi_ratio <- Fungi_ratio %>% filter( get(columnMetadataValue) == SelectMetadataValue )
  }
  return(Fungi_ratio)
}

#Fungi_ratio <- get.Fungi_ratio.2groupments(Fungi_TPM_filtered_Annot, WhereToTest="KOG_group",WhereToTest2="LifeStyle",SelectMetadataValue="all", columnMetadataValue="all")





# GroupToTest <- unique(Fungi_ratio$tax_phylum) 
# GroupToTest <- GroupToTest[!GroupToTest %in% c("0", NA)]
# 
# byCond= "Layer"
# condValues= c("OS","OM")


collect_TPM_prc_byGroup <- function(Fungi_ratio, GroupToTest, WhereToTest, byCond, condValues ){
  perc_group <- data.frame(group=GroupToTest)
  colnames(perc_group) <- c("group")
  for(l in condValues){
    tag <- paste(l,"perc",sep="_")
    tmp <- Fungi_ratio %>% filter(get(byCond) == l)%>%
      group_by(samples) %>% mutate(somme=sum(tpm_counts)) %>% ungroup() %>%
      mutate(perc=(tpm_counts/somme)*100 ) %>%  # left_join(treatment_metadata)
      group_by(get(WhereToTest)) %>% summarise(median=median(perc))
    colnames(tmp) <- c("group",tag)
    perc_group <- merge(tmp,perc_group, by="group")
    
    tag <- paste(l,"TPM",sep="_")
    tmp <- Fungi_ratio %>% filter(get(byCond) == l)%>%
      group_by(get(WhereToTest))  %>% summarise(median=round(median(tpm_counts), digits = 0), .groups = "keep")
    colnames(tmp) <- c("group",tag)
    perc_group <- merge(tmp,perc_group, by="group")
  }
  return(perc_group)
}
# 
# collect_TPM_prc_byGroup(Fungi_ratio, GroupToTest, WhereToTest="tax_phylum", byCond="Layer", condValues=c("OS","OM"))
# collect_TPM_prc_byGroup(Fungi_ratio, GroupToTest, WhereToTest="tax_phylum", byCond="tree_species", condValues=c("Quercus","Abies","Picea"))

t.test_paired_Layer <- function(Fungi_ratio, GroupToTest, WhereToTest){
  t.res <- data.frame(group=character(), p.value=numeric() )
  for(grp in GroupToTest){
    t <- Fungi_ratio %>% filter(get(WhereToTest) ==grp)%>%# filter(tree_species == "Abies") %>%
      select(tree_number,Layer,tpm_counts) %>%
      pivot_wider(names_from = Layer, values_from = tpm_counts) %>%
      filter(!is.na(OS)) %>%
      mutate(difference=OS-OM) 
    
    
    # test si suit une loi normale
    shaptest.res <-  shapiro.test(t$difference)
    bartlett.res <- bartlett.test(list(t$OM, t$OS))
    if(shaptest.res$p.value > 0.05 & bartlett.res$p.value > 0.05) #si > 0.05 loi normale et homogénéité des variances
    {   res.test <- t.test(t$OS, t$OM, paired=TRUE)
    }
    else{
      res.test <-  wilcox.test(t$OS, t$OM, paired=TRUE)
    }
    tmpz <- data.frame(group=grp, p.value=res.test$p.value )
    t.res <- bind_rows(t.res,tmpz)
    rm(tmpz)
    #print(c(grp, res.test$p.value))
  }
  ## post-hoc correction p.value
  if(length(t.res$group) < 20 ){
    t.res$pajust <- p.adjust(t.res$p.value, method = "bonferroni")  #fdr ; bonferroni
  }else{
    t.res$pajust <- p.adjust(t.res$p.value, method = "fdr") #fdr ; bonferroni
  }
  
  t.res <- t.res %>% select(group, pajust)
  return(t.res)
}

#t.res <- t.test_paired_Layer(Fungi_ratio, GroupToTest, WhereToTest)

getLayerTestResFormated <- function(Fungi_ratio, GroupToTest, WhereToTest ){
  
  perc_group <- collect_TPM_prc_byGroup(Fungi_ratio, GroupToTest, WhereToTest, byCond="Layer", condValues=c("OS","OM"))
  t.res <- t.test_paired_Layer(Fungi_ratio, GroupToTest, WhereToTest)
  
  dataframeRes <- merge(t.res,perc_group)
  
  dataframeRes$Difference_OSvsOM[dataframeRes$OS_perc-dataframeRes$OM_perc > 0 ] <- "OS > OM"
  dataframeRes$Difference_OSvsOM[dataframeRes$OS_perc-dataframeRes$OM_perc < 0 ] <- "OS < OM"
  dataframeRes$Difference_OSvsOM[dataframeRes$pajust > 0.05 ] <- "N.S."
  dataframeRes$FC_OSvsOM <- round(dataframeRes$OS_perc/dataframeRes$OM_perc, digits = 1)
  dataframeRes$FC_OSvsOM[dataframeRes$pajust > 0.05 ] <- "N.S."
  
  return(dataframeRes)
}

# LayerTestResFormated <- getLayerTestResFormated(Fungi_ratio, GroupToTest, WhereToTest)
# 
# 
# 
# tmp <- colnames(Fungi_TPM)
# taxaTEST <- tmp[grepl(pattern = "PD_", tmp)]
# WhereToTest <- "tax_phylum"
# rm(WhereToTest)
# 
# 
# 
# WhereToTest="PD_Order"
# arbre="Picea"
# 
# Taxonomy <- c("PD_Phylum","PD_Class","PD_Order","PD_Family","PD_Genus")
#  KOG <- c("KOG_class")
#  WhereToTest <- "KOG_class"
#   
# for(WhereToTest in Taxonomy ) {
#   Fungi_ratio <- get.Fungi_ratio(Fungi_TPM, WhereToTest,SelectMetadataValue="all", columnMetadataValue="all")
#   GroupToTest= unique(Fungi_ratio[,as.character(WhereToTest)])
#   colnames(GroupToTest) <- "group"
#   GroupToTest <- GroupToTest$group[!GroupToTest$group %in% c("0", NA)]
#   LayerTestResFormatedSummary=   data.frame(group=GroupToTest)
#   colnames(LayerTestResFormatedSummary) <- c("group")
#   
#       for(arbre in c("Quercus","Abies","Picea", "all")){
#             
#             Fungi_ratio <- get.Fungi_ratio(Fungi_TPM, WhereToTest,SelectMetadataValue=arbre, columnMetadataValue="tree_species")
#             
#             GroupToTest= unique(Fungi_ratio[,as.character(WhereToTest)])
#             colnames(GroupToTest) <- "group"
#             GroupToTest <- GroupToTest$group[!GroupToTest$group %in% c("0", NA)]
#             
#             LayerTestResFormated <- getLayerTestResFormated(Fungi_ratio, GroupToTest, WhereToTest)
#             
#             nom=paste(WhereToTest,arbre, "LayerTestResFormated.tsv",sep="__")
#             write.table(LayerTestResFormated, file = nom, row.names = F, sep="\t")
#             
#             nom2 <- paste(arbre,WhereToTest,sep="__")
#             tmp <- LayerTestResFormated %>% select(group, Difference_OSvsOM)
#             LayerTestResFormatedSummary<- merge(LayerTestResFormatedSummary, tmp, all = T)
#             names(LayerTestResFormatedSummary)[length(names(LayerTestResFormatedSummary))]<-nom2
#       }
#   
#   nom3 <- paste(WhereToTest, "LayerTestSummary.tsv",sep="__")
#   percColum <- paste(arbre, "__OS_perc", sep="")
#   tmp <- LayerTestResFormated %>% select(group, OS_perc)
#   LayerTestResFormatedSummary<- merge(LayerTestResFormatedSummary, tmp, all = T)
#   names(LayerTestResFormatedSummary)[length(names(LayerTestResFormatedSummary))]<-percColum
#   
#   LayerTestResFormatedSummary
#   
#   write.table(LayerTestResFormatedSummary, file = nom3, row.names = F, sep="\t")
#   
# }
# 




# ######## By tree_type #######
# #### Anova / kruskall ####
# 
# 
#  WhereToTest="PD_Order"
# # arbre="Picea"
# sol_layer="OS"
# column_layer="Layer"
# 
# group2testArbre="tree_species"
# 
# Fungi_ratio <- get.Fungi_ratio(Fungi_TPM, WhereToTest,SelectMetadataValue=sol_layer, columnMetadataValue=column_layer)
# GroupToTest <- unique(Fungi_ratio$PD_Order) 
# GroupToTest <- GroupToTest[!GroupToTest %in% c("0", NA)]
# 
# 
# group2testPossibilities <- sort(unlist(unique(Fungi_ratio[group2testArbre])))



get_Lists_per_grp <- function(filtered_Fungi_ratio, group2testPossibilities, group2testArbre) {
    N <- length(group2testPossibilities)
    tpm_Lists <- vector("list", N)
    n=0
    for(possibilite in group2testPossibilities){
          n=n+1
          tmp <- filtered_Fungi_ratio %>% filter(get(group2testArbre) == possibilite) 
          a <- tmp$tpm_counts
          tpm_Lists[[n]] <- a
              }
    tpm_Lists <- setNames(tpm_Lists, group2testPossibilities)
    return(tpm_Lists) } ## tpm_Lists <-  list(Oak=c(1,3,3),Abies=c(1,4,5),Spruce=c(1,4,6))


test_normalite_Homoscedasticite <- function(filtered_Fungi_ratio, group2testPossibilities, group2testArbre){
  normalite=T
  tpm_Lists <- get_Lists_per_grp(filtered_Fungi_ratio, group2testPossibilities, group2testArbre)
  for(possibilite in tpm_Lists){
    shaptest.res <- shapiro.test(possibilite)
    if(shaptest.res$p.value < 0.05){normalite=F} 
  }
  bartlett.res <- bartlett.test(tpm_Lists)
  #si > 0.05 loi normale et homogénéité des variances
  if(normalite == T & bartlett.res$p.value > 0.05){normalite_Homoscedasticite=T }
  else{normalite_Homoscedasticite=F }  
  return(normalite_Homoscedasticite)
} ### return T or F 


test_Kruskall_Anova <- function(filtered_Fungi_ratio, GroupToTest,WhereToTest,group2testPossibilities,group2testArbre) {
      t.res <- data.frame(group=character(), 
                                 normal=character(), 
                                 test_method=character()#,p.value=numeric()#,
                                 )
      t.res.tmp <- data.frame(p.value=numeric())
      for(grp in GroupToTest){
        filtered_Fungi_ratio <- Fungi_ratio %>% filter(get(WhereToTest) ==grp)
        normalite_Homoscedasticite <- test_normalite_Homoscedasticite(filtered_Fungi_ratio, group2testPossibilities, group2testArbre)
        tpm_Lists <- get_Lists_per_grp(filtered_Fungi_ratio, group2testPossibilities, group2testArbre)
        kruskal_test_res <- kruskal.test(tpm_Lists)
        t.res <- bind_rows(t.res,c(group=grp,
                                       normal=normalite_Homoscedasticite,
                                       test_method=kruskal_test_res$method#, 
                                       #p.value=kruskal_test_res$p.value #,
                                       ))
                                 
        t.res.tmp <- bind_rows(t.res.tmp,c(p.value=kruskal_test_res$p.value))
        }
       t.res2 <- bind_cols(t.res,t.res.tmp)
      
      if(length(t.res$group) < 20 ){
        t.res$padjust <- p.adjust(t.res2$p.value, method = "bonferroni")  #fdr ; bonferroni
      }else{
        t.res$padjust <- p.adjust(t.res2$p.value, method = "fdr") #fdr ; bonferroni
      }
       t.res$signif <- "ns"      
       t.res$signif[t.res$padjust < 0.05] <- "*"
       t.res$signif[t.res$padjust < 0.01] <- "**"
       t.res$signif[t.res$padjust < 0.001] <- "***"
       
       return(t.res)
}

#test_Kruskall_Anova.res <- test_Kruskall_Anova(filtered_Fungi_ratio, GroupToTest,WhereToTest,group2testPossibilities,group2testArbre)


test_Dunn_postKruskall <- function(filtered_Fungi_ratio,WhereToTest, test_Kruskall_Anova.res ) {
          
        signifDiff <- test_Kruskall_Anova.res %>% filter(padjust < 0.05)
        
        t.res.tmp <- data.frame(group=character(), 
                                Abies_Picea.padjust=numeric(), 
                                Abies_Quercus.padjust=numeric(),
                                Picea_Quercus.padjust=numeric())
        
        for(grp in signifDiff$group){
           filtered_Fungi_ratio <- Fungi_ratio %>% filter(get(WhereToTest) ==grp)
           dunn_test.res <- dunn_test(filtered_Fungi_ratio, tpm_counts ~ tree_species, p.adjust.method = "bonferroni") 
           
           dunn_test.res$group <- grp 
           dunn_test.res$comp <- paste(dunn_test.res$group1, dunn_test.res$group2, sep="_")
           dunn_test.res$signif <- paste(dunn_test.res$comp,".signif", sep="")
           dunn_test.res$comp <- paste(dunn_test.res$comp,".padjust", sep="")
           dunn_test.res2 <- dunn_test.res %>% select(group,comp,  p.adj ) %>%
                                              pivot_wider(names_from = comp, values_from = p.adj)
           
           dunn_test.res3 <- dunn_test.res %>% select(group,signif,  p.adj.signif ) %>%
             pivot_wider(names_from = signif, values_from = p.adj.signif)  
           
           dunn_test.res <- left_join(dunn_test.res2, dunn_test.res3)
           
           t.res.tmp <-  bind_rows(t.res.tmp,dunn_test.res )
           }
        t.res.update <- left_join(test_Kruskall_Anova.res, t.res.tmp)
        return(t.res.update)
}



#t.res.update <- test_Dunn_postKruskall(filtered_Fungi_ratio,WhereToTest, test_Kruskall_Anova.res ) 


#TPM_prc_byGroup <-  collect_TPM_prc_byGroup(Fungi_ratio, GroupToTest, WhereToTest, byCond="tree_species", condValues=c("Quercus","Abies","Picea"))

#t.res.update2 <- left_join(t.res.update, TPM_prc_byGroup)




Split_Annot_In_Multiple_Rows <- function(Fungi_TPM, WhereToTest){
  Fungi_TPM_filtered <- Fungi_TPM %>% #dplyr::select(-starts_with("NR_"),-starts_with("Myc"),-ID,-ctg_length, -sum_contig_count) %>% 
    #dplyr::select("uniqueID",all_of(WhereToTest), starts_with("PD"),starts_with("LifeStyle"), starts_with("O") ) %>%
    #filter(!is.na(get(WhereToTest))) %>%  filter(get(WhereToTest) != "NoHit" ) %>%
    #filter(KOG_group %in% c("CELLULAR PROCESSES AND SIGNALING", "INFORMATION STORAGE AND PROCESSING", "METABOLISM")) %>%
    #mutate(KOG_group=case_when(KOG_group != "Unknown" ~ KOG_group, TRUE ~ "UNDEFINED"  )) %>%
    separate_rows(WhereToTest,sep = "', '") %>%
    mutate(!!WhereToTest:=sub("^'", "", get(WhereToTest))) %>%
    mutate(!!WhereToTest:=sub("'$", "", get(WhereToTest))) %>%
    mutate(!!WhereToTest:=trimws(get(WhereToTest),which = c("both"),whitespace = "[ \t\r\n']" ))
  
  Fungi_TPM_filtered_Annot <- unique(Fungi_TPM_filtered)
  return(Fungi_TPM_filtered_Annot)
}

get_dfstatTPM<- function(Fungi_ratio, WhereToTest){
  dfstatTPM <- Fungi_ratio %>%
    group_by(get(WhereToTest)) %>% summarise(max=max(tpm_counts), min=(min(tpm_counts)), mean=mean(tpm_counts), median=median(tpm_counts))
  colnames(dfstatTPM)[colnames(dfstatTPM) == "get(WhereToTest)"] <- WhereToTest
  return(dfstatTPM )
}

get_dfstatPERC<- function(Fungi_ratio, WhereToTest){
  dfstatPERC <- Fungi_ratio %>%  group_by(get(WhereToTest), samples) %>% summarise(tpm=sum(tpm_counts)) %>% ungroup() %>%
    group_by(samples) %>% mutate(somme=sum(tpm)) %>% ungroup() %>%
    mutate(perc=(tpm/somme)*100 ) 
  colnames(dfstatPERC)[colnames(dfstatPERC) == "get(WhereToTest)"] <- WhereToTest
  dfstatPERC <-dfstatPERC %>% group_by(get(WhereToTest)) %>%
    summarise(max= round(max(perc),digits =1), min=round(min(perc),digits =1), mean=round(mean(perc),digits =1), median=round(median(perc),digits =1))
  colnames(dfstatPERC)[colnames(dfstatPERC) == "get(WhereToTest)"] <- WhereToTest
  return(dfstatPERC)
}

get_dfstatTPM_LifeStyle<- function(Fungi_ratio, WhereToTest){
  dfstatTPM <- Fungi_ratio %>%
    group_by(get(WhereToTest), LifeStyle) %>% 
    summarise(max=max(tpm_counts), min=(min(tpm_counts)), mean=mean(tpm_counts), median=median(tpm_counts), .groups = "keep")
  colnames(dfstatTPM)[colnames(dfstatTPM) == "get(WhereToTest)"] <- WhereToTest
  return(dfstatTPM)
}

get_dfstatPERC_LifeStyle <- function(Fungi_ratio, WhereToTest){
  dfstatPERC <- Fungi_ratio %>%  group_by(get(WhereToTest),LifeStyle, samples) %>% summarise(tpm=sum(tpm_counts)) %>% ungroup() %>%
    group_by(samples) %>% mutate(somme=sum(tpm)) %>% ungroup() %>%
    mutate(perc=(tpm/somme)*100 ) 
  colnames(dfstatPERC)[colnames(dfstatPERC) == "get(WhereToTest)"] <- WhereToTest
  dfstatPERC <-dfstatPERC %>% group_by(get(WhereToTest),LifeStyle) %>%
    summarise(max= round(max(perc),digits =1), min=round(min(perc),digits =1), mean=round(mean(perc),digits =1), median=round(median(perc),digits =1))
  colnames(dfstatPERC)[colnames(dfstatPERC) == "get(WhereToTest)"] <- WhereToTest
  return(dfstatPERC)
}








