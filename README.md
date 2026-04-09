# oCEM: Automatic detection and analysis of overlapping co-expressed gene modules
#### I. Introduction
---
When it comes to the co-expressed gene module detection, its typical challenges consist of overlap between identified modules and local co-expression in a subset of biological samples. Recent studies have reported that the decomposition methods are the most appropriate for solving these challenges. In this study, we represent an R tool, termed Overlapping CoExpressed gene Module (oCEM), which possesses those methods with a wholly automatic analysis framework to help non-technical users to easily perform complicated statistical analyses and then gain robust results. We also develop a novel auxiliary statistical approach to select the optimal number of principal components using a permutation procedure. Three example datasets are used, related to human breast cancer, mouse metabolic syndrome, and E.coli gene expression compendium, to enable the illustration of the straightforward use of the tool. Computational experiment results show that oCEM outperforms state-of-the-art techniques in the ability to detect biologically relevant modules additionally.

#### II. Understanding the tool and Data Structure
---
The following are parameters included in overlapCEM and their role:
- data: a data frame or matrix. `data` has its rows are samples and its columns are genes.

- clinical: a data frame or matrix. Input data serve to perform Pearson's correlations between each identified module and each clinical feature. It includes its rows are samples, and its columns are clinical features of your choice.

- ncomp: positive integer. The optimal number of principal components. It should be >= 1.

- standardize: logical. If your `data` are not standardized, just feed `T` or `TRUE` to this parameter. Default value is `T`.

- method: string. Post-processing methods. Allowed values are `ICA-FDR`, `ICA-Zscore`, or `IPCA-FDR`. Default value is `ICA-Zscore`.

- cex.text: numeric. Change the font size of texts in cells of the heatmap showing correlations between each identified module and each clinical feature. Default value is 0.7.

- P: positive integer. The number of permutations used to estimate the null distribution. Default value is `1000`.

- cores: positive integer. The number of CPU cores used for parallel computation during permutation. Default value is `1`.

- type: character string. The method used to determine the optimal number of principal components. Allowed values are `intersection` and `perm`. Default value is intersection.

- alpha: numeric. The significance threshold used when `type = "perm"`. Default value is 0.05.

- adjust.method: character string. The method used for multiple testing correction of permutation p-values. This argument is passed to p.adjust(). Default value is BH. Note that this parameter is only applied when `type = "perm"`.

- seed: numeric or NULL. A random seed for reproducibility. If NULL, no seed is set. Default value is NULL.

- verbose: logical. Show the running time to complete running the whole pipeline. Default value is `T`.

Please download datasets [data_n_code](https://github.com/huynguyen250896/oCEM/tree/main/data_n_code) and read [Additional File 1](https://github.com/huynguyen250896/oCEM/blob/main/Additional%20File%201.pdf) (highly recommended) as examples to well grasp oCEM's easy-to-meet format and its usage.

#### III. Pipeline
---
![Figure](https://imgur.com/lPoY1UX.png)
**Figure:** Pipeline of the package oCEM.

#### IV. Implementation
---
Use the following command to install directly from GitHub;
```sh
devtools::install_github("huynguyen250896/oCEM", dependencies = T)
```
Call the nescessary libraries;
```sh
x = c("oCEM", "dplyr", "dynamicTreeCut", "flashClust","Hmisc",
  "WGCNA", "moments", "fastICA", "tidyr", "fdrtool", "mixOmics",
  "cluster", "purrr", "parallel")
lapply(x, require, character.only = TRUE)
```
running example:
```sh
# oCEM
num_pc <- optimizeCOM(data = exp, cores = 5)
# >> oCEM suggests choosing the optimal number of components is: 9

cem <-overlapCEM(data = exp, clinical = clinicalEXP, ncomp = num_pc)
```

#### V. What's new
---
- 2026-04-09: We introduced a permutation-based approach (`type = "perm"`) to determine the optimal number of principal components in optimizeCOM. In addition to the original intersection rule introduced in the paper, users can now choose between two strategies:

• `intersection`: selects the number of components based on the crossing point between the observed explained variance curve and the average permutation-based reference curve. This approach is computationally efficient and consistent with the original implementation, but it is heuristic in nature.

• `perm`: selects the number of components using permutation-based p-values. For each component, its explained variance is compared against a null distribution generated from permuted data to assess statistical significance. This approach is more statistically grounded, although it may require more computational time.

We also introduced the `adjust.method` parameter to control multiple testing correction when `type = "perm"` (default: `"BH"`). While this provides stricter statistical control, it may be overly conservative when the number of components is large or the number of permutations `P` is small. For exploratory analysis, users may consider setting `adjust.method = "none"`.

In addition, a `seed` parameter was added to ensure reproducibility of the permutation procedure. Setting a fixed seed allows users to obtain stable and consistent results across runs.

- 2023-10-13: I made a bad decision that required the users to input both the mRNA and clinical data into `overlapCEM` to be able to run the tool successfully. Now, I made the input clinical data optional, meaning that the mRNA data is the only data for the tool to run. Besides, I refactored the codes comprehensively that would make them readable more (and hope that it runs more rapidly also!). Besides, default value to `optimizeCOM`'s `method` parameter set to `ICA-Zscore` makes both life scientists and bioinformatics scientists not confused about what to select and serves to compare the performance of oCEM with that of other tools. This decision was based on the results of a wonderful paper [[1](https://www.nature.com/articles/s41467-018-03424-4)].
- 2023-10-08: Users now can set the number of cores to the `optimizeCOM` algorithm on their own using its new argument `cores`, meaning that they can parallely perform the algorithm and get the optimal number of PCs more rapidly . Unfortunately, this feature is not available to Window users this time!

#### VI. Citation
---
Please kindly cite the following paper (and Star this Github repository if you find this tool of interest) if you use the tool in this repo: </br>
```sh
Reference Type: Journal Article
Author: Nguyen, Quang-Huy
Le, Duc-Hau
Year: 2022
Title: oCEM: Automatic detection and analysis of overlapping co-expressed gene modules
Journal: BMC Genomics
Volume: 23
Issue: 1
Pages: 39
Date: 2022/01/08
ISSN: 1471-2164
DOI: 10.1186/s12864-021-08072-5
```
Feel free to contact [Quang-Huy Nguyen](https://github.com/huynguyen250896) <huynguyen96.dnu AT gmail DOT com> for any questions about the code and results.
