#library
library(WGCNA)
library(overlappingCGM)

# ----------------------o0o---------------------- #
#                 Pre-process data
# ----------------------o0o---------------------- #
#load raw file
exp = read.table("LiverFemale3600.csv", header = T, check.names = F, sep=",")
cli = read.table("ClinicalTraits.csv", header = T, check.names = F, sep=",", row.names = 1)

#remove missing gene names and duplicated genes in exp
exp = exp[which(exp$gene_symbol != "0"),]

dup = duplicated(exp$gene_symbol)
exp = exp[which(dup == FALSE),]

#turn exp1 into satisfactory format of DrGA
#overlappingCGM requires data whose rows are samples and columns are genes.
exp = exp %>%
  dplyr::select(-c(substanceBXH, LocusLinkID, ProteomeID, cytogeneticLoc,
                   CHROMOSOME, StartPosition, EndPosition)) %>%
  tibble::remove_rownames() %>%
  tibble::column_to_rownames('gene_symbol') %>%
  drop_na() %>% t()

#detect outliers
sampleTree = hclust(dist(exp), method = "average")
par(cex = 0.6);par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="",
     cex.lab = 1.5,cex.axis = 1.5, cex.main = 2)

# Plot a line to show the cut
abline(h = 12.2, col = "red");
# Determine cluster under the line
clust = cutreeStatic(sampleTree, cutHeight = 12.2, minSize = 10)
table(clust)
# clust 1 contains the samples we want to keep.
keepSamples = (clust==1)
exp = exp[keepSamples, ]

#match mouses that share between cli versus exp
cli = cli[cli$Mice %in% rownames(exp),]
cli = cli %>%
  remove_rownames() %>%
  tibble::column_to_rownames('Mice') %>% 
  dplyr::select(-c(Number, sex, Mouse_ID, Strain, DOB, parents, Western_Diet,
                   Sac_Date, comments, Note))

#how the clinical traitsrelate to the sample dendrogram.
# Re-cluster samples
sampleTree2 = hclust(dist(exp), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(cli, signed = FALSE);
# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree2, traitColors,groupLabels = names(cli),
                    main = "Sample dendrogram and trait heatmap")
#white means a low value, red a high value, and grey a missing entry.
#Only keep several clinical features with red color
cli = cli %>%
  dplyr::select(weight_g, length_cm, ab_fat,
                total_fat, UC, FFA, Glucose, 
                LDL_plus_VLDL)

#make sure that mice that share between exp and cli are included at their rows and
#in exactly the same order 
all(rownames(exp) == rownames(cli))
#[1] FALSE
exp = exp[rownames(cli), ]

#check dimension
dim(exp)
#[1] 134 2281
dim(cli)
#[1] 134  8


# ----------------------o0o---------------------- #
#                      WGCNA
# ----------------------o0o---------------------- #
wgnca = blockwiseModules(exp, power = 6,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       verbose = 3)

# ----------------------o0o---------------------- #
#                  overlappingCGM
# ----------------------o0o---------------------- #
daBEST(data = exp)
# >> overlappingCGM suggests choosing the optimal number of components is: 18 
# >> Both ICA and IPCA-FDR are appropriate for your case. Please use another more stringent approach to make the best decision. 

cgm=overlappingCGM(data = exp, clinical = cli, ncomp = 18,
                        method = 'ICA-Zscore')

# ----------------------o0o---------------------- #
#                     iWGNCA
# ----------------------o0o---------------------- #
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
enableWGCNAThreads() ### Allowing parallel execution

#we realize that the soft-thresholding of 5 is best fit for the data
softPower = 5; 
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
# 1   2   3   4   5   6   7   8   9  10  11  12 
# 341 311 272 230 216 211 188 168 112  84  80  68 

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
cor(wgcna$MEs, cgm$patterns)

# iWGCNA versus overlappingCGM
cor(MEs, cgm$patterns)



