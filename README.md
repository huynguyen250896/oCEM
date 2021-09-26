# oCEM: Automatic detection and analysis of overlapping co-expressed gene modules
#### I. Introduction
---
When it comes to the co-expressed gene module detection, its typical challenges consist of overlap between identified modules and local co-expression in a subset of biological samples. Recent studies have reported that the decomposition methods are the most appropriate for solving these challenges. In this study, we represent an R tool, termed Overlapping CoExpressed gene Module (oCEM), which possesses those methods with a wholly automatic analysis framework to help non-technical users to easily perform complicated statistical analyses and then gain robust results. We also develop a novel auxiliary statistical approach to select the optimal number of principal components using a permutation procedure. Three example datasets are used, related to human breast cancer, mouse metabolic syndrome, and E.coli gene expression compendium, to enable the illustration of the straightforward use of the tool. Computational experiment results show that oCEM outperforms state-of-the-art techniques in the ability to detect biologically relevant modules additionally .

#### II. Understanding the tool and Data Structure
---
The following are parameters included in overlapCEM and their role:
- data: a data frame or matrix. `data` has its rows are samples and its columns are genes.

- clinical: a data frame or matrix. Input data serve to perform Pearson's correlations between each identified module and each clinical feature. It includes its rows are samples, and its columns are clinical features of your choice.

- ncomp: positive integer. The optimal number of principal components. It should be >= 2.

- standardize: logical. If your `data` are not standardized, just feed `T` or `TRUE` to this parameter. Default value is `T`.

- method: Post-processing methods. Allowed values are `ICA-FDR`, `ICA-Zscore`, or `IPCA-FDR`.

- cex.text: numeric. Change the font size of texts in cells of the heatmap showing correlations between each identified module and each clinical feature. Default value is 0.7.

Please download datasets [data_n_code](https://github.com/huynguyen250896/oCEM/tree/main/data_n_code) and read [Additional File 1](https://github.com/huynguyen250896/oCEM/blob/main/Additional%20File%201.pdf) (highly recommended) as examples to well grasp oCEM's easy-to-meet format and its usage.

#### III. Pipeline
---
![Figure](https://imgur.com/lPoY1UX.png)
**Figure:** Pipeline of the package oCEM.

#### IV. Implementation
---
Use the following command to install directly from GitHub;
```sh
devtools::install_github("huynguyen250896/oCEM")
```
Call the nescessary libraries;
```sh
x = c("oCEM", "dplyr", "dynamicTreeCut", "flashClust","Hmisc",
  "WGCNA", "moments", "fastICA", "tidyr", "fdrtool", "mixOmics",
  "cluster", "purrr")
lapply(x, require, character.only = TRUE)
```
running example:
```sh
# oCEM
optimizeCOM(data = exp)
# >> oCEM suggests choosing the optimal number of components is: 9
# >> oCEM also suggests using ICA for your case. 

cem=overlapCEM(data = exp, clinical = clinicalEXP, ncomp = 9,
                        method = 'ICA-Zscore', cex.text = 1.0)
```

#### V. Citation
---
Please kindly cite the following paper (and Star this Github repository if you find this tool of interest) if you use the tool in this repo: </br>
```sh
Reference Type: Journal Article
Author: Nguyen, Quang-Huy
Le, Duc-Hau
Year: 2021
Title: oCEM: Automatic detection and analysis of overlapping co-expressed gene modules
Journal: bioRxiv
Date: 2021/03/15
DOI: 10.1101/2021.03.15.435373
```
Feel free to contact [Quang-Huy Nguyen](https://github.com/huynguyen250896) <huynguyen96.dnu AT gmail DOT com> for any questions about the code and results.
