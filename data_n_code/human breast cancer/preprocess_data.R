# ----------------------o0o---------------------- #
#                 Load library
# ----------------------o0o---------------------- #

#library
devtools::install_github("huynguyen250896/oCEM", force = T) #NOTE: Download all dependencies of the tool!

x=c("oCEM", "dplyr", "dynamicTreeCut", "flashClust","Hmisc", "WGCNA", "moments", "fastICA", "tidyr", "fdrtool", 
    "mixOmics", "cluster", "purrr", "parallel")
lapply(x, require, character.only = TRUE)

# ----------------------o0o---------------------- #
#                 Pre-process data
# ----------------------o0o---------------------- #
#load raw data
exp = read.table('data_mRNA_median_Zscores.txt', sep = '\t', check.names = FALSE, header = TRUE, row.names = NULL)
clinical = read.table('data_clinical_patient.txt', sep = '\t', check.names = FALSE, header = TRUE, row.names = 1,fill=TRUE)

#31 identified driver genes
driver=c("MAP2K4", "ARID1A", "PIK3CA", "TBX3", "MAP3K1", "TP53", "AKT1", "GATA3", "CDH1", "RB1", "CDKN1B", "NCOR1", "CDKN2A", "ERBB2", "KRAS", "BRCA2", "BAP1", "PTEN", "CBFB", "KMT2C", "RUNX1", "NF1", "PIK3R1", "ERBB3", "FOXO3", "SMAD4", "GPS2", "AGTR2", "ZFP36L1", "MEN1","SF3B1")
length(driver) #31 driver genes

#only keep the 31 driver genes in exp
exp=exp %>%
  dplyr::filter(.$Hugo_Symbol %in% driver) %>%
  tibble::column_to_rownames('Hugo_Symbol') %>%
  dplyr::select(-Entrez_Gene_Id)

#check dimension
dim(exp) # 31 1904
dim(clinical) # 2509   21

#match patients sharing between exp versus clinical
#exp are the matrix whose rows are samples and columns are genes
exp = exp[,intersect(colnames(exp), rownames(clinical))] %>% t()
clinicalEXP = clinical[intersect(rownames(exp), rownames(clinical)), ]

#preprocess clinicalEXP
clinicalEXP = clinicalEXP %>%
  tibble::rownames_to_column('sample') %>%
  dplyr::select(c(sample, LYMPH_NODES_EXAMINED_POSITIVE, NPI, stage)) %>%
  tibble::column_to_rownames('sample')

colnames(clinicalEXP)[1:3] =  c("lymph", "npi", "stage")
str(clinicalEXP)
# 'data.frame':	1904 obs. of  5 variables:
# $ lymph : int  1 5 8 1 0 1 0 2 0 6 ...
# $ npi   : num  4.04 6.03 6.03 5.04 3.05 ...
# $ stage : int  2 2 3 2 2 2 1 2 2 4 ...

# ----------------------o0o---------------------- #
#                      WGCNA
# ----------------------o0o---------------------- #
wgcna = blockwiseModules(exp, power = 6,
                       TOMType = "unsigned", minModuleSize = 10,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       verbose = 3)

# ----------------------o0o---------------------- #
#                  oCEM
# ----------------------o0o---------------------- #
num_cor <- optimizeCOM(data = exp, cores = 5)
# >> oCEM suggests choosing the optimal number of components is: 9
# >> oCEM also suggests using ICA for your case.

cem <- overlapCEM(data = exp, clinical = clinicalEXP, ncomp = num_cor, cex.text = 1.0)

# ----------------------o0o---------------------- #
#                     iWGNCA
# ----------------------o0o---------------------- #
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
enableWGCNAThreads() ### Allowing parallel execution

#we choose the soft-thresholding of 6 based on our prior work
softPower = 6; 
adjacency = adjacency(exp, power = softPower,
                      type = "signed");

# Turn adjacency into topological overlap
TOM = TOMsimilarity(adjacency, TOMType = "signed");
dissTOM = 1-TOM

# Hierichical clustering
#Seek the optimal agglomeration method
#methods to assess
m <- c( "average", "single", "complete", "ward")
names(m) <- c( "average", "single", "complete", "ward")
# function to compute agglomerative coefficient
set.seed(25896)
ac <- function(x) {
  agnes(exp, method = x)$ac
}

map_dbl(m, ac) # Agglomerative coefficient of each agglomeration method
# average    single  complete      ward 
# 0.4136737 0.3510251 0.4815638 0.5563910  

#assign gene names from adjacency to dissTOM
rownames(dissTOM) = rownames(adjacency) 
colnames(dissTOM) = colnames(adjacency) 
# Call the hierarchical clustering function
geneTree = hclust(as.dist(dissTOM), method = "ward.D2");

#We set the minimum module size at 10:
minModuleSize = 10;

# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            minClusterSize = minModuleSize);
table(dynamicMods)
# dynamicMods
# 1  2 
# 16 15 

# Convert numeric lables into colors
moduleColors = labels2colors(dynamicMods)

#Define numbers of genes and samples
nGenes = ncol(exp);
nSamples = nrow(exp);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(exp, moduleColors)$eigengenes
MEs = orderMEs(MEs0)


# ----------------------o0o---------------------- #
#                     Comparison
# ----------------------o0o---------------------- #
# comparison between WGCNA and overlappingCGM and between iWGCNA and overlappingCGM
# WGCNA versus overlappingCGM
cor(wgcna$MEs, cem$patterns)

# iWGCNA versus overlappingCGM
cor(MEs, cem$patterns)


