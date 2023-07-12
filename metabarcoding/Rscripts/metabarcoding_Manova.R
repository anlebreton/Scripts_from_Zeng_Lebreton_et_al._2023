library(phyloseq)
library(vegan)

## Setting variables
## The Phyloseq object (format rdata)
 phyloseq <- load("../1_data/Yulong_ITS1.Rdata")

## The beta diversity distance matrix file
 distance <- "../4_betaDiv/ITS1/bray.tsv"

## The experiment variable that you want to analyse
varExp <- "Tree_species+Season+Layer+Cluster"

## Create input and parameters dataframe
 params <- data.frame( "phyloseq" = phyloseq, "distance" = distance, "varExp" = varExp)
## Load data
# the phyloseq object, nammed data in FROGSSTAT Phyloseq Import data
#load(params$phyloseq)
# Convert sample_data to data.frame
metadata <- as(sample_data(data), "data.frame") 

# the distance matrix file
A        <- read.table(file=params$distance, row.names=1)
dist     <- as.dist(A)

## Multivariate ANOVA performed with adonis
adonis <- paste('adonis(dist ~ ', params$varExp, ', data = metadata, perm = 9999)')
#adonis(formula = dist ~ Tree_species * Layer, data = metadata,      permutations = 9999) 
eval(parse(text = adonis))


# 
# library(olsrr)
# 
# mod.iris <- lm(cbind(pH, Sepal.Width, Petal.Length, Petal.Width)
#                ~ Species, data=iris)
# 
# 
# # model <- lm(mpg ~ disp + hp + wt + qsec, data = mtcars)
#  ols_test_normality(res)
# # model <- lm(mpg ~ disp + hp + wt + qsec, data = mtcars)
#  ols_test_correlation(res)


#Div.fit1=lmer (shannon/chao1 ~ Tree species*Season * Soil type + (1|block),REML=TRUE,data=sp))
