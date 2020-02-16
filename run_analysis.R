#############################################################################################
## PROOF OF CONCEPT: autoencoder use as pre-processing step for cytometry data clustering  ##
## David Novak (davidnovakcz@hotmail.com), Faculty of Science, Charles University          ##
## 16 Feb 2020                                                                             ##
## !!! Required packages: flowCore, keras, FlowSOM, tidyverse, Rtsne, uwot, gridExtra !!!  ##
## (Consider pulling a rocker image with tidyverse, TensorFlow & Keras pre-installed to    ##
## if you don't want to clutter your machine.)                                             ##
## !!! Required input data: https://web.stanford.edu/~samusik/Panorama%20BM%201-10.zip !!! ##
## (Samusik, Nikolay et al. “Automated mapping of phenotype space with single-cell data.”  ##
## Nature methods vol. 13,6 (2016): 493-6. doi:10.1038/nmeth.3863)                         ##
#############################################################################################

## Load packages

pkgs    <- rownames(installed.packages())
reqs    <- c('flowCore', 'keras', 'FlowSOM', 'tidyverse', 'Rtsne', 'uwot', 'gridExtra')
if (length(reqs[!reqs %in% pkgs] -> missing) > 0) stop(paste0('Missing packages: ', paste(missing, collapse = ', '), '.'))
for (pkg in reqs) library(pkg, character.only = T)

## Create sample dataset

set.seed(1)

input_fcs    <- 'Panorama/BM2_cct_normalized_01_non-Neutrophils.fcs'
ff           <- flowCore::read.FCS(input_fcs)
all_labs     <- read.csv('Panorama/population_assignments.txt', sep = '\t', header = FALSE)

all_labs$V1  <- as.character(all_labs$V1)

# Get indices of rows in the table corresponding to labelled events
idcs_in_csv  <- grep(gsub('[.]fcs', '', basename(input_fcs)), all_labs$V1)

# Retrieve the corresponding indices in the FCS file rows
idcs_in_fcs  <- as.numeric(sapply(strsplit(all_labs$V1[idcs_in_csv], ' '), function(x) tail(x, 1))) + 1
labs         <- all_labs$V2[idcs_in_csv]

# Generate expression matrix of labelled events
efcs         <- ff@exprs[idcs_in_fcs, ]
efcs         <- efcs[, !colnames(efcs) %in% c('Time', 'Cell_length', 'beadDist')] # omit some parameters

# Transform the data: cofactor 1/5 is good for CyTOF data (around 1/120 for flow)
efcs <- asinh(efcs/5)
nc   <- ncol(efcs)

rm(ff)

## Define function for building and compiling an AE model

ae_projection <- function(latent_dim = 20L, input_data) {
  
  # Divide training and test data
  n_train      <- min(10000, floor(nrow(input_data) * .85))
  idcs_train   <- sort(sample(1:nrow(input_data), n_train, replace = FALSE))
  idcs_test    <- sample((1:nrow(input_data))[-idcs_train], 3000, replace = FALSE)
  
  labs.u       <- unique(as.character(labs))
  legend       <- data.frame(label = labs.u, code = 1:length(labs.u), stringsAsFactors = FALSE)
  
  train_data   <- efcs[idcs_train, ]
  test_data    <- efcs[idcs_test, ]
  train_labels <- labs[idcs_train]
  test_labels  <- labs[idcs_test]
  
  model <- keras_model_sequential()
  
  model %>%
    # Encoder
    layer_dropout(0.2, input_shape = nc) %>% # dropout layer for regularisation
    layer_dense(units = 128, activation = 'selu') %>%
    layer_dense(units = 256, activation = 'selu') %>%
    layer_dense(units = 512, activation = 'selu') %>%
    layer_dense(units = 64,  activation = 'selu') %>%
    layer_dense(units = 32,  activation = 'selu') %>%
    layer_dense(units = latent_dim, activation = 'selu', name = 'latent_variable') %>%
    # Decoder
    layer_dense(units = 32,  activation = 'selu') %>%
    layer_dense(units = 64,  activation = 'selu') %>%
    layer_dense(units = 512, activation = 'selu') %>%
    layer_dense(units = 256, activation = 'selu') %>%
    layer_dense(units = 128, activation = 'selu') %>%
    layer_dense(units = nc)
  
  model %>% compile(
    loss      = 'mean_squared_error',
    optimizer = 'adam'
  )
  
  model %>% fit(
    x = train_data,
    y = train_data,
    epochs = 20
  )
  
  latent <- keras_model(
    inputs  = model$input,
    outputs = get_layer(model, 'latent_variable')$output
  )
  
  proj <- data.frame(predict(latent, input_data))
}

## Create projections with different latent space dimensions

projections <- list(proj.45 = ae_projection(45L, efcs),
                    proj.40 = ae_projection(40L, efcs),
                    proj.35 = ae_projection(35L, efcs),
                    proj.30 = ae_projection(30L, efcs),
                    proj.25 = ae_projection(25L, efcs),
                    proj.20 = ae_projection(20L, efcs),
                    proj.15 = ae_projection(15L, efcs),
                    proj.10 = ae_projection(10L, efcs),
                    proj.5  = ae_projection(5L, efcs),
                    proj.2  = ae_projection(2L, efcs))

## For original data & for all projections:
#### train a SOM with 30 x 30 Kohonen layer & apply hierarchical clustering for meta-cluster generation

n_pops             <- length(unique(na.omit(labs[idcs_in_fcs])))

efcs.clustered.som <- FlowSOM::BuildSOM(
                        FlowSOM::ReadInput(
                          flowCore::flowFrame(efcs), scale = TRUE),
                        colsToUse = 1:nc, xdim = 30, ydim = 30)

efcs.metaclusters  <- FlowSOM::metaClustering_consensus(
                        efcs.clustered.som$map$codes, k = n_pops)

proj.clustered.som <- lapply(projections, function(x)
                        FlowSOM::BuildSOM(
                          FlowSOM::ReadInput(
                            flowCore::flowFrame(as.matrix(x)), scale = TRUE),
                          colsToUse = 1:ncol(x), xdim = 30, ydim = 30))

proj.metaclusters  <- lapply(proj.clustered.som, function(x)
                          FlowSOM::metaClustering_consensus(x$map$codes, k = n_pops))

## For original data & for all projections:
#### apply k-means clustering

efcs.clustered.kmeans <- kmeans(efcs, iter.max = 50L, centers = n_pops)
proj.clustered.kmeans <- lapply(projections, function(x) kmeans(x, iter.max = 50L, centers = n_pops))

## Define function for computing F1-score

F_score <- function(y_true, y_pred) {
  TP <- sum( y_true &  y_pred)
  TN <- sum(!y_true & !y_pred)
  FP <- sum(!y_true &  y_pred)
  FN <- sum( y_true & !y_pred)
  
  precision <- TP / (TP + FP)
  recall    <- TP / (TP + FN)
  
  2 * precision * recall / (precision + recall)
}

## Define function to evaluate clustering using pairwise F1-scoring

eval.clustering <- function(n_sample, y_true, y_pred1, ...) {
  params <- list(y_true, y_pred1, ...)
  
  # Sample pairs of events
  s1     <- sample(x = 1:length(y_true), size = n_sample, replace = TRUE)
  s2     <- sample(x = 1:length(y_true), size = n_sample, replace = TRUE)
  pairs  <- cbind(s1, s2)
  
  # Remove duplicates
  hashes <- apply(pairs, 1, function(x) paste0(x, collapse = ' '))
  dupes  <- which(duplicated(hashes))
  if (length(dupes) > 0) pairs <- pairs[-dupes, ]
  
  # Evaluate differences in event annotation (same-population vs. different-population) and clustering (same-cluster vs. different-cluster)
  diffs <- lapply(params, function(x) apply(pairs, 1, function(y) x[y[1]] == x[y[2]]))
  
  lapply(2:length(diffs), function(idx) F_score(diffs[[1]], diffs[[idx]]))
}

## Retrieve SOM meta-cluster assignments per event

cl.som     <- lapply(1:length(proj.clustered.som),
                     function(idx) proj.metaclusters[[idx]][proj.clustered.som[[idx]]$map$mapping[, 1]])
cl_all.som <- c(list(y_true  = as.numeric(labs),
                     y_pred1 = efcs.metaclusters[efcs.clustered.som$map$mapping[, 1]]),
                     cl.som)

## Compute pairwise F1-scoring for SOM of original data & latent variable representations

params     <- c(n_sample = list(1000000), cl_all.som)
e1         <- do.call(eval.clustering, params)
e2         <- do.call(eval.clustering, params)
e3         <- do.call(eval.clustering, params)
e4         <- do.call(eval.clustering, params)
e5         <- do.call(eval.clustering, params)

f1.som <- rbind(data.frame(latent_dim = c(48, 45, 40, 35, 30, 25, 20, 15, 10, 5, 2), f_score = unlist(e1)),
                data.frame(latent_dim = c(48, 45, 40, 35, 30, 25, 20, 15, 10, 5, 2), f_score = unlist(e2)),
                data.frame(latent_dim = c(48, 45, 40, 35, 30, 25, 20, 15, 10, 5, 2), f_score = unlist(e3)),
                data.frame(latent_dim = c(48, 45, 40, 35, 30, 25, 20, 15, 10, 5, 2), f_score = unlist(e4)),
                data.frame(latent_dim = c(48, 45, 40, 35, 30, 25, 20, 15, 10, 5, 2), f_score = unlist(e5)))

## Retrieve k-means cluster assignments per event

cl.kmeans     <- lapply(proj.clustered.kmeans, function(x) x$cluster)
cl_all.kmeans <- c(list(y_true  = as.numeric(labs),
                        y_pred1 = efcs.clustered.kmeans$cluster),
                        cl.kmeans)

## Compute pairwise F1-scoring for k-means clustering of original data & latent variable representations

params <- c(n_sample = list(1000000), cl_all.kmeans)
g1     <- do.call(eval.clustering, params)
g2     <- do.call(eval.clustering, params)
g3     <- do.call(eval.clustering, params)
g4     <- do.call(eval.clustering, params)
g5     <- do.call(eval.clustering, params)

f1.kmeans <- rbind(data.frame(latent_dim = c(48, 45, 40, 35, 30, 25, 20, 15, 10, 5, 2), f_score = unlist(g1)),
                   data.frame(latent_dim = c(48, 45, 40, 35, 30, 25, 20, 15, 10, 5, 2), f_score = unlist(g2)),
                   data.frame(latent_dim = c(48, 45, 40, 35, 30, 25, 20, 15, 10, 5, 2), f_score = unlist(g3)),
                   data.frame(latent_dim = c(48, 45, 40, 35, 30, 25, 20, 15, 10, 5, 2), f_score = unlist(g4)),
                   data.frame(latent_dim = c(48, 45, 40, 35, 30, 25, 20, 15, 10, 5, 2), f_score = unlist(g5)))

f1.som    <- cbind(method = rep('SOM',     nrow(f1.som)),    f1.som)
f1.kmeans <- cbind(method = rep('k-means', nrow(f1.kmeans)), f1.kmeans)

f1.scores        <- rbind(f1.kmeans, f1.som)
f1.scores$method <- factor(f1.scores$method)

## Save F1-scores object

saveRDS(f1.scores, 'f1_scores.RDS')

## Generate plot of F1-scores 

p.f1 <- ggplot(f1.scores, aes(x = factor(latent_dim), y = f_score)) + geom_point(aes(colour = method), size = 10) +
        theme_light() + ggtitle('Clustering pairwise F1-scores', subtitle = paste0('Original data dimensionality: ', nc)) +
        labs(y = 'F1-score', x = 'data dimension')
  
## Generate plot of 2-dimensional AE projection

df   <- mutate(as.data.frame(projections$proj.2), class = as.factor(labs))
p.ae <- ggplot(df, aes(x = X1, y = X2, colour = class)) + geom_point(size = .5) +
        ggtitle('Latent variable projection') +  labs(y = 'component 2', x = 'component 1') +
        guides(colour = guide_legend(override.aes = list(size = 10)), legend.text = element_text(size = 3))

## Generate plot of 2-dimensional t-SNE projection

proj.tsne <- Rtsne::Rtsne(efcs, dims = 2)
df.tsne   <- mutate(as.data.frame(proj.tsne$Y), class = as.factor(labs))
p.tsne    <- ggplot(df.tsne, aes(x = V1, y = V2, colour = class)) + geom_point(size = .5) +
             ggtitle('t-SNE projection') +  labs(y = 'component 2', x = 'component 1') +
             guides(colour = guide_legend(override.aes = list(size = 10)), legend.text = element_text(size = 3))

## Generate plot of 2-dimensional UMAP projection

proj.umap <- uwot::umap(efcs)
df.umap   <- mutate(as.data.frame(proj.umap), class = as.factor(labs))
p.umap    <- ggplot(df.umap, aes(x = V1, y = V2, colour = class)) + geom_point(size = .5) +
             ggtitle('UMAP projection') +  labs(y = 'component 2', x = 'component 1') +
             guides(colour = guide_legend(override.aes = list(size = 10)), legend.text = element_text(size = 3))

## Generate plot of 2-dimensional PCA projection

proj.pca  <- prcomp(efcs, rank. = 2)
df.pca    <- mutate(as.data.frame(proj.pca$x), class = as.factor(labs))
p.pca     <- ggplot(df.pca, aes(x = PC1, y = PC2, colour = class)) + geom_point(size = .5) +
             ggtitle('PCA projection') +  labs(y = 'component 2', x = 'component 1') +
             guides(colour = guide_legend(override.aes = list(size = 10)), legend.text = element_text(size = 3))

## Generate a summary plot comparing 2-dimensional projections

p.comparison <- grid.arrange(p.ae   + theme_linedraw() + theme(legend.position = 'none'),
                             p.pca  + theme_linedraw() + theme(legend.position = 'none'),
                             p.tsne + theme_linedraw() + theme(legend.position = 'none'),
                             p.umap + theme_linedraw() + theme(legend.position = 'none'))

## Save workspace data

save.image('workspace.RData')
