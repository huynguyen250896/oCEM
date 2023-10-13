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
exp = read.table('ecoli_data.tsv', sep = '\t', check.names = FALSE, header = TRUE, row.names = NULL) %>% t()
gene = read.table('ecoli_gene_names.tsv', sep = '\t', check.names = FALSE, header = TRUE, row.names = NULL, col.names = NULL)

#assign gene names to gene expression data
colnames(exp) = gene$row.names

#detect outliers
sampleTree = hclust(dist(exp), method = "average")
par(cex = 0.6);par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="",
     cex.lab = 1.5,cex.axis = 1.5, cex.main = 2)

# Plot a line to show the cut
abline(h = 80, col = "red");
# Determine cluster under the line
clust = cutreeStatic(sampleTree, cutHeight = 80, minSize = 10)
table(clust)
# clust 1 contains the samples we want to keep.
keepSamples = (clust==1)
exp = exp[keepSamples, ]

#remove genes whose std in expression levels < 1
exp1 = apply(exp, 2, function(x) sd(x)) %>%
  as.data.frame()
exp1 = exp1 %>% dplyr::filter(. > 1)

exp = exp[, colnames(exp) %in% rownames(exp1)]
exp <- exp[!duplicated(rownames(exp)),]
#ensuring rows are samples and its columns are genes.
dim(exp) # 752 samples 572 genes. 

# ----------------------o0o---------------------- #
#                      WGCNA
# ----------------------o0o---------------------- #
wgcna = blockwiseModules(exp, power = 6,
                         TOMType = "unsigned", minModuleSize = 10,
                         reassignThreshold = 0, mergeCutHeight = 0.25,
                         numericLabels = TRUE, pamRespectsDendro = FALSE,
                         verbose = 3)

# ----------------------o0o---------------------- #
#                     oCEM
# ----------------------o0o---------------------- #
num_PC <- optimizeCOM(data = exp, cores = 5)
# >> oCEM suggests choosing the optimal number of components is: 25 
# >> Both ICA and IPCA-FDR are appropriate for your case. Please use a more stringent approach to make the best decision. 

cem <- overlapCEM(data = exp, ncomp = num_PC,
                  cex.text = 1.0)
#We tried running oCEM with ICA-Zscore (default) with a stricter threshold, Z-score = 1.5

# ----------------------o0o---------------------- #
#                     iWGNCA
# ----------------------o0o---------------------- #
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
enableWGCNAThreads() ### Allowing parallel execution

# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(exp, powerVector = powers, verbose = 5,
                        networkType = "signed")
# Plot the results:
sizeGrWindow(9, 5)
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.81,col="red")
#we choose the soft-thresholding of 14 based on our prior work
softPower = 14; 
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
#  1   2   3   4   5   6   7   8   9 
# 221 120  74  52  25  25  24  20  11

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
# comparison between WGCNA and oCEM and between iWGCNA and oCEM
# WGCNA versus oCEM
cor(wgcna$MEs, cem$patterns)

# iWGCNA versus oCEM
cor(MEs, cem$patterns)


