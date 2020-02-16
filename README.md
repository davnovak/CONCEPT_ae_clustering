# PROOF OF CONCEPT: dimensionality reduction via autoencoder as pre-processing for cytometry data clustering

David Novak

Faculty of Science, Charles University

16 Feb 2020

---

## Overview

This repository contains scripts provided as proof-of-concept for reducing dimensionality of cytometry data using a lower-dimensional projection in latent space of an auto-encoder (AE). Previously published work [1, 2, 3] showcases the use of auto-encoders with objective functions designed for clustering of events. Moreover, the use of AEs for dimensionality reduction of single-cell data (scRNA-seq, cytometry) is featured in some recent publications [4, 5].

On the other hand, the approach shown here is extremely simplistic. It serves as an initial attempt to evaluate the change in quality of clustering with using the dimensionality reduction step. A portion of annotated data provided with [6] is used for building auto-encoder models with various latent space dimensionality (45, 40, 35, 30, 25, 20, 15, 10, 5, 2). Afterward, both the original data and each lower-dimensional projection thereof are clustered using a 30-by-30 FlowSOM grid (and subsequent meta-clustering) [7] and *k*-means clustering. The quality of each set of cluster assignments is then evaluated using pairwise F1-scores:

```
N ........ event count
L ........ vector of length N; true label assignments per each event
C ........ vector of length N; cluster assignments per each event
n_iter ... number of sampling iterations

P, Q <- empty numeric vectors

For i from 1 to n_iter:
	Sample p, q from N, such that p != q
	Append p to P
	Append q to Q

Remove duplicated pairs (p, q) from P, Q

Y_true, Y_predicted <- empty logical vector

For i in 1 to Length(P):
	y_true      <- (L[P[i]] == L[Q[i]])
	y_predicted <- (C[P[i]] == C[Q[i]])
	Append y_true to Y_true
	Append y_predicted to Y_predicted

TP, TN, FP, FN <- 0 # true positives, true negatives, false positives & false negatives

For i in 1 to Length(P):
	TP += 1 if Y_true[i] AND Y_pred[i]
	TN += 1 if (NOT Y_true[i]) AND (NOT Y_pred[i])
	FP += 1 if (NOT Y_true[i]) AND Y_pred[i]
	FN += 1 if Y_true[i] AND (NOT Y_pred[i])

precision <- TP / (TP + FP)
recall    <- TP / (TP + FN)

F1 <- 2 * precision * recall / (precision + recall)

Output F1

```

We also generate plots showing that even an (overly) simplistic AE-generated 2-dimensional embedding does a somewhat fine job of separating out annotated cell populations (we show a comparison with PCA, t-SNE [8] and UMAP [9]).

## Results

Will add these soon.

## How to run the example yourself

To run the example, you will need to install `R` and the packages `flowCore`, `FlowSOM`, `tidyverse`, `Rtsne`, `uwot` and `gridExtra`.

In addition, the example uses the `Keras` API for `R`. I recommend to set up Keras to work with the `TensorFlow` backend. Installation instructions can be found at *https://keras.rstudio.com/*.

Pull this GitHub repository and download input FCS files from *https://web.stanford.edu/~samusik/Panorama%20BM%201-10.zip*. Place the data into the repository folder. For the analysis, you need a directory named `Panorama` with the FCS files and the `population_assignments` TXT file in it.

Then, run the script `run_analysis.R` to produce the above images. The script is heavily commented, allowing for modifications (including use of different input data).

## References

[1] Chunfeng Song, Feng Liu, Yongzhen Huang, Liang Wang, and Tieniu Tan. Auto-encoder ba- sed data clustering. In José Ruiz-Shulcloper and Gabriella Sanniti di Baja, editors, Progress in Pattern Recognition, Image Analysis, Computer Vision, and Applications, pages 117–124, Berlin, Heidelberg, 2013. Springer Berlin Heidelberg.

[2] Kai Tian, Shuigeng Zhou, and Jihong Guan. Deepcluster: A general clustering framework based on deep learning. In Michelangelo Ceci, Jaakko Hollmén, Ljupčo Todorovski, Celine Vens, and Sašo Džeroski, editors, Machine Learning and Knowledge Discovery in Databases, pages 809–825, Cham, 2017. Springer International Publishing.

[3] Junyuan Xie, Ross B. Girshick, and Ali Farhadi. Unsupervised deep embedding for clustering analysis. CoRR, abs/1511.06335, 2015.

[4] Hu, Qiwen, and Casey S Greene. “Parameter tuning is a key part of dimensionality reduction via deep variational autoencoders for single cell RNA transcriptomics.” Pacific Symposium on Biocomputing. Pacific Symposium on Biocomputing vol. 24 (2019): 362-373.

[5] Szubert, B., Cole, J.E., Monaco, C. et al. Structure-preserving visualisation of high dimensional single-cell datasets. Sci Rep 9, 8914 (2019). https://doi.org/10.1038/s41598-019-45301-0

[6] Samusik, Nikolay et al. “Automated mapping of phenotype space with single-cell data.” Nature methods vol. 13,6 (2016): 493-6. doi:10.1038/nmeth.3863

[7] Van Gassen, S., Callebaut, B., Van Helden, M.J., Lambrecht, B.N., Demeester, P., Dhaene, T. and Saeys, Y. (2015), FlowSOM: Using self‐organizing maps for visualization and interpretation of cytometry data. Cytometry, 87: 636-645. doi:10.1002/cyto.a.22625

[8] Van Der Maaten LJP, Hinton GE. Visualizing high-dimensional data using t-SNE. J Mach Learn Res 2008;92579–2605.

[9] 1. McInnes L, Healy J, Melville J. UMAP : Uniform Manifold Approximation and Projection for Dimension Reduction [Epub ahead of print].
