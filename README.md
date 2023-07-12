# Scripts from Zeng, Lebreton et al. (2023)

Q.Zeng #, A. Lebreton #, L. Auer, X. Man, L. Jia, G. Wang, S. Gong, V. Lombard, M. Buée, G. Wu, Y. Dai, Z. Yang and F.M. Martin

#These authors contributed equally to this work
# Stable functional structure despite high taxonomic variability across fungal communities in soils of old-growth montane forests



# metabarcoding
---- Preprocessing ---
 
The preprocessing of the files were done in Galaxy on UseGalaxy.fr instance. 

Workflows are available here:

ITS2 (and ITS1) primers:  https://usegalaxy.fr/u/a.lebreton/w/its2metabarcodingfrogsprevalence

18S primers: https://usegalaxy.fr/u/a.lebreton/w/amfmetabarcodingfrogsprevalence

Starting from quality checked reads, it leads to the table of OTU abundance.
Then, FROGS Tree and FROGSSTAT Phyloseq Import Data (using normalisation if not previously performed) were used.
These OTU abundance tables and generated Rdata were the starting point of the script referenced here. 


 # metatranscriptomic
 ---- Preprocessing ---
 
A fungal TPM count table was generated as describded in the article. 
