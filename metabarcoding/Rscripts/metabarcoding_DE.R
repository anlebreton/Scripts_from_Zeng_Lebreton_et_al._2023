## Import packages
library(phyloseq)
library(ggplot2)
library(ape)

## Settin variables
## The OTU abundance matrix with taxonomy annotation file (biom format)
biomfile <- "../1_data/Yulong_ITS1_normalised_abundance.biom"

## The sample metadata file(TSV format)
 samplefile <- "../0_Metadata/Yulong_sampling_ITS1_metadata.txt"

## (optional) the OTU tree file (nwk format). Write "None" if you do not have any tree
 treefile <- "Yulong_ITS1_tree.nwk"

## The ordered taxonomic levels stored in BIOM. Each level is separated by one space.
## default : "Kingdom Phylum Class Order Family Genus Species"
 ranks <- "Kingdom Phylum Class Order Family Genus Species"

## Do you want to normalise your data ? "True" or "False"
 normalisation <- "False"

## Create input and parameters dataframe
 params <- data.frame( "biomfile" = biomfile, "samplefile" = samplefile, "tree" = treefile, "ranks" = ranks, "normalisation" = normalisation)


## Import data
biomfile <- params$biomfile
data     <- import_biom(biomfile)
sampledata <- read.csv(params$samplefile, sep = "\t", row.names = 1)

# if taxonomy starts with k__ it means that its Greengenes like format
# import need to be done using parse_taxonomy_greengenes function
# in this case user taxonomic rank names are ignored
tax      <- tax_table(data)[[1]]
if ((gregexpr('k__', tax))[[1]][1]>0) { 
  cat("Warning : Taxonomic affiliations come from Greengenes database, user specified ranks names are ignored.")
  data <- import_biom(biomfile, parseFunction = parse_taxonomy_greengenes)
} else {
  ## else, custumize rank name with the user specified ranks variable
  new_rank <- as.list(strsplit(params$ranks, " ")[[1]])
  colnames(tax_table(data)) <- new_rank
}
#Warning : Taxonomic affiliations come from Greengenes database, user specified ranks names are ignored.
## add sample name to metadata, as SampleID variable
sampledata$SampleID <- rownames(sampledata)
sample_data(data) <- sampledata

## add tree metadata if available
if (params$treefile != "None"){
  treefile <- read.tree(params$treefile)
  phy_tree(data) <- treefile
}

## change de sample metadata order as in input samplefile
for ( variable in sample_variables(data)){
  variable.order = as.vector(unique(sampledata[,variable]))
  sample_data(data)[,variable] <- factor(get_variable(data, variable),levels=variable.order)
}

## remove empty samples
empty_samples <- sample_names(data)[which(sample_sums(data)==0)]
sample_to_keep <- sample_names(data)[which(sample_sums(data)>0)]
data <- prune_samples(sample_to_keep, data)

## abundance normalisation
if(params$normalisation){ data <- rarefy_even_depth(data, rngseed = 1121983)}

## save phyloseq object in Rdata file
save(data, file=params$outputRdata)





# Import input files
phyloseq_input <- get(load(params$phyloseq_data))
dds_input     <- get(load(params$dds))
result_names <- DESeq2::resultsNames(dds_input)
doc <- ""
if(length(result_names) == 2 || params$mod1 == "None"){
  results <- DESeq2::results(dds, name = result_names[1])
} else if(length(result_names)>=3){
  results <- DESeq2::results(dds_input, contrast = c(params$var, params$mod1, params$mod2) )
  doc <- paste0("You chose to compare ", params$mod1, " to the reference modality ", params$mod2, ". This implies that a positive log2FoldChange means more abundant in ", params$mod1, " than in ", params$mod2, ".")
  cat(doc)
}

#You chose to compare OM to the reference modality HF. This implies that a positive log2FoldChange means more abundant in OM than in HF.



#Then we extract significant OTUs at the p-value adjusted threshold level (after correction) and enrich results with taxonomic informations and sort taxa by pvalue.

da.otus_table <- results %>% as_tibble(rownames = "OTU") %>% 
  inner_join(tax_table(phyloseq_input) %>% as.data.frame() %>% mutate(OTU = taxa_names(phyloseq_input)), by = "OTU") %>% 
  filter(padj < params$padj) %>%  
  arrange(padj)

#Using customize datatable we can explore the table of differentially expressed OTUs

DT_foldchange <- da.otus_table %>% datatable(filter = "top",                # ==> add filter for each column
                                             extensions = 'Buttons' , 
                                             options = list(dom = '<fBtlip>', 
                                                            scrollX = TRUE, 
                                                            # lengthMenu = list(c(5,10,25,50,100,-1),c(5,10,25,50,100,"All")), # ==> trop couteux en mÃ©moire, mÃªme sur petit jeu
                                                            buttons = list(list(
                                                              extend = 'collection',
                                                              buttons = c('csv', 'excel'),
                                                              text = 'Download')
                                                            )
                                             )) %>% 
  formatSignif(columns = c("baseMean", "log2FoldChange", "lfcSE", "stat", "padj", "pvalue"), digits = 6) %>% 
  formatStyle(columns = "log2FoldChange", color = DT::styleInterval(0, c('#759493', '#C6792B')))

DT_foldchange






da_pie <- data.frame(
  OTUs=c("Differentially Abundant (log-fold change > 0)",
         "Differentially Abundant (log-fold change < 0)",
         "Not Differentiallly Abundant"
  ),
  Number=c(sum(da.otus_table$log2FoldChange > 0), 
           sum(da.otus_table$log2FoldChange < 0), 
           nrow(results) - nrow(da.otus_table)
  )
)

Pie <- ggplot(da_pie, aes(x="", y=Number, fill=OTUs)) + 
  geom_bar(stat = "identity") + 
  coord_polar("y", start=0) + 
  scale_fill_manual(values=c('#8EADAC', '#DE9F73', '#C0C0C0')) + 
  geom_text(aes(label = Number), position = position_stack(vjust = 0.5), color = "white") +
  labs(x = NULL, y = NULL, fill = NULL, title = "Pie chart to view OTUs number of Differential Abundance test") +
  theme_classic() + theme(
    axis.line = element_blank(), 
    axis.text = element_blank(),
    axis.ticks = element_blank())

Pie


## Volcano plot by DESeq2
# Add OTU names
# Add the -log10 pvalue 
# Add the pre-calculated log2 fold change
da_volcano <- data.frame(
  otu      = row.names(results),
  evidence = -log10(results$padj), 
  lfc      = results$log2FoldChange)

# Remove rows that have NA values
da_volcano <- na.omit(da_volcano)
# add a threshol line
y_axix_volcano_line <- -log10(params$padj)

# Modify dataset to add new coloumn of colors
da_volcano <- da_volcano %>%
  mutate(
    color = case_when(
      lfc > 0 & evidence > y_axix_volcano_line ~ "More", 
      lfc < 0 & evidence > y_axix_volcano_line ~ "Less", 
      TRUE                                     ~ "Equal"
    )
  )

# Color corresponds to fold change directionality
volcano_plot <- ggplot(da_volcano, aes(x = lfc, y = evidence)) + 
  geom_point(aes(color = factor(color)), size = 1.75, alpha = 0.8, na.rm = T) + # add gene points
  theme_bw(base_size = 16) +                                                    # clean up theme
  theme(legend.position = "none") +                                             # remove legend 
  ggtitle(label = "Volcano Plot", subtitle = "Colored by effect sign") +        # add title
  xlab(expression(log[2]("FoldChange"))) +                                      # x-axis label
  ylab(expression(-log[10]("adjusted p-value"))) +                              # y-axis label
  geom_vline(xintercept = 0, colour = "grey80", linetype = 2) +                                # add line at 0
  geom_hline(aes(yintercept = y_axix_volcano_line), yintercept = y_axix_volcano_line, colour = "grey80", linetype = 2) +
  annotate(geom = "text", 
           label = paste("padj =", params$padj), 
           x = min(da_volcano$lfc), 
           y = y_axix_volcano_line + 0.25, 
           size = 4,
           colour = "black",
           vjust = 0,
           hjust = 0) + # add pvalue threshold           
  scale_color_manual(values = c("More" = "#C6792B", "Less" = "#759493", "Equal" = "#C0C0C0")) # change colors

# Plot figure
volcano_plot + scale_y_continuous(trans = "log1p")





## Could do with built-in DESeq2 function:
DESeq2::plotMA(results, ylim = c(-10,10),
               ylab = "log2 fold change", alpha=params$padj,
               cex  = 0.8, colNonSig="#C0C0C0", colSig="#A2A32F", colLine="#C0C0C0",
               main = "Post Normalisation DESeq2: MA plot of log2FoldChange",
               legend = TRUE)



da.otus_table <- da.otus_table %>% arrange(log2FoldChange)

if ( length(result_names)==2 || params$mod1 != "None"){
  
  subtitle <- result_names[2]
  if ( params$mod1!= "None" ){
    phyloseq_input <- prune_samples(get_variable(phyloseq_input, params$var) %in% c(params$mod1, params$mod2), phyloseq_input)
    subtitle <- paste(sep="", params$var, "_", params$mod1, "_vs_", params$mod2)
  }
  
  ##Heatmap plot with 2 conditions
  hplot <- plot_heatmap(prune_taxa(da.otus_table$OTU, phyloseq_input), 
                        taxa.label = NULL,
                        taxa.order = da.otus_table$OTU, low = "yellow", high = "red", na.value = "white") + 
    facet_grid(as.formula(paste("~", params$var)), scales = "free_x") + 
    ggtitle("Heatmap plot of DA otus, between 2 conditions", subtitle = subtitle)
  ## Quantitative variable
} else {
  idx <- order(get_variable(phyloseq_input, params$var))
  my.sample.order <- sample_names(phyloseq_input)[idx]
  ##Heatmap plot with all conditions of variable
  hplot <- plot_heatmap(prune_taxa(da.otus_table$OTU, phyloseq_input), 
                        taxa.label = NULL,
                        taxa.order = da.otus_table$OTU, sample.order = my.sample.order, 
                        low = "yellow", high = "red", na.value = "white") + 
    ggtitle(paste0("Heatmap plot of DA otus with samples ordered according to ", params$var))
}
# ggplotly(hplot)
hplot





