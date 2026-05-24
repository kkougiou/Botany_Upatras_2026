label_dist <- function(d, labels) {
  attr(d, "Labels") <- as.character(labels)
  d
}



# --- 1a. Convert your dist object to bioregion pairwise format ---------------
convert_dist_to_bioregion <- function(dist_obj, metric_name = "Simpson",
                                      nb_species = NA) {
  if (!inherits(dist_obj, "dist")) stop("Object must be of class 'dist'.")
  labels <- attr(dist_obj, "Labels")
  mat    <- as.matrix(dist_obj)
  df <- data.frame(
    Site1 = rep(labels, times = length(labels)),
    Site2 = rep(labels, each  = length(labels)),
    value = as.vector(mat),
    stringsAsFactors = FALSE
  )
  df <- df[df$Site1 < df$Site2, ]
  colnames(df)[3] <- metric_name
  class(df) <- c("bioregion.pairwise.metric", "data.frame")
  attr(df, "type")       <- "dissimilarity"
  attr(df, "nb_sites")   <- attr(dist_obj, "Size")
  attr(df, "nb_species") <- nb_species
  return(df)
}


bioregion_metrics_fixed <- function (bioregionalization, comat, map = NULL, col_bioregion = NULL) 
{
  if (inherits(bioregionalization, "bioregion.clusters")) {
    if (inherits(bioregionalization$clusters, "data.frame")) {
      has.clusters <- TRUE
      clusters <- bioregionalization$clusters
      if (ncol(clusters) > 2) {
        stop(paste0("This function is designed to be applied on a single ", 
                    "bioregionalization."), call. = FALSE)
      }
    }
    else {
      if (bioregionalization$name == "hclu_hierarclust") {
        stop(paste0("No clusters have been generated for your hierarchical ", 
                    "tree, please extract clusters from the tree before using ", 
                    "bioregionalization_metrics().\n", "See ?hclu_hierarclust or ?cut_tree."), 
             call. = FALSE)
      }
      else {
        stop(paste0("bioregionalization does not have the expected type of ", 
                    "'clusters' slot."), call. = FALSE)
      }
    }
  }
  else {
    stop(paste0("This function is designed to work on bioregion.clusters ", 
                "objects and on a site x species matrix."), call. = FALSE)
  }
  bioregion:::controls(args = NULL, data = comat, type = "input_matrix")
  if (!is.null(map)) {
    if (!("sf" %in% class(map))) {
      stop(paste0("map must be a 'sf' spatial data.frame with bioregions ", 
                  "and sites."), call. = FALSE)
    }
    if (!("data.frame" %in% class(map))) {
      stop(paste0("map must be a 'sf' spatial data.frame with bioregions ", 
                  "and sites."), call. = FALSE)
    }
    if (ncol(map) < 3) {
      stop(paste0("map must have at least 3 columns: sites, bioregions and ", 
                  "geometry."), call. = FALSE)
    }
    if (is.null(col_bioregion)) {
      stop(paste0("col_bioregion must be defined ", "it is the column position ", 
                  "of the bioregion."), call. = FALSE)
    }
  }
  if (!is.null(col_bioregion)) {
    if (is.null(map)) {
      warning(paste0("col_bioregion is defined but is not considered since ", 
                     "map is set to NULL."))
    }
    else {
      bioregion:::controls(args = col_bioregion, data = NULL, type = "strict_positive_integer")
      map_test <- map
      sf::st_geometry(map_test) <- NULL
      if (inherits(map_test[, col_bioregion], "logical")) {
        stop("There is no bioregion in the Bioregion column.", 
             call. = FALSE)
      }
      rm(map_test)
    }
  }
  bioregion_df <- geometry <- NULL
  if (bioregionalization$inputs$bipartite == TRUE) {
    clusters <- clusters[which(attributes(clusters)$node_type == 
                                 "site"), ]
  }
  bioregion_df <- data.frame()
  for (j in 1:bioregionalization$cluster_info$n_clust) {
    focal_j <- unique(clusters[, 2])[j]
    sites_j <- clusters[which(clusters[, 2] == focal_j), 
                        "ID"]
    comat_j <- comat[sites_j, colSums(comat[sites_j, , drop = FALSE]) > 
                       0, drop = FALSE]
    comat_not_j <- comat[-which(rownames(comat) %in% sites_j), 
                         , drop = FALSE]
    comat_not_j <- comat_not_j[, colSums(comat_not_j) > 
                                 0, drop = FALSE]
    endemic_j <- sum(!(colnames(comat_j) %in% colnames(comat_not_j)))
    bioregion_df <- rbind(bioregion_df, data.frame(Bioregion = focal_j, 
                                                   Site_number = length(sites_j), Species_number = ncol(comat_j), 
                                                   Endemics = endemic_j, Percentage_Endemic = 100 * 
                                                     endemic_j/ncol(comat_j)))
  }
  if (length(unique(bioregion_df$Bioregion)) != bioregionalization$cluster_info$n_clust) {
    warning("Not all bioregions are in the output.")
  }
  if (!is.null(map)) {
    colnames(map)[col_bioregion] <- "Bioregion"
    if (dplyr::n_distinct(map$Bioregion) == 1) {
      bioregion_df <- data.frame(bioregion_df, Coherence = 100)
    }
    else if (dplyr::n_distinct(map$Bioregion) > 1) {
      spatial_coherence_df <- data.frame()
      for (j in 1:dplyr::n_distinct(map$Bioregion)) {
        bioregion_j <- unique(map[which(map$Bioregion == 
                                          unique(map$Bioregion)[j]), ]$Bioregion)
        map_j <- map[which(map$Bioregion == bioregion_j), 
        ]
        map_j <- dplyr::summarise(map_j, geometry = sf::st_union(geometry))
        map_j <- sf::st_cast(map_j, "POLYGON")
        map_j <- dplyr::mutate(map_j, ID = dplyr::row_number())
        map_j$area <- as.numeric(sf::st_area(map_j))
        sp_coherence_j <- 100 * max(map_j$area, na.rm = TRUE)/sum(map_j$area, 
                                                                  na.rm = TRUE)
        spatial_coherence_df <- rbind(spatial_coherence_df, 
                                      data.frame(Bioregion = bioregion_j, Coherence = sp_coherence_j))
      }
      bioregion_df <- dplyr::left_join(bioregion_df, spatial_coherence_df, 
                                       by = "Bioregion")
    }
    else {
      stop("There is no bioregion in the Bioregion column.", 
           call. = FALSE)
    }
  }
  return(bioregion_df)
}


calc_affiliation_evoreg_fixed <- function(phylo.comp.dist, groups) {
  if (inherits(phylo.comp.dist, "dist") == FALSE) {
    stop("phylo.comp.dist might be from class dist")
  }
  if (length(groups) != nrow(as.matrix(phylo.comp.dist))) {
    stop("Phylogenetic Distance Matrix and group vectors might have the same sites")
  }
  n.groups <- length(as.numeric(levels(groups)))
  comm.groups <- groups
  Gs <- lapply(1:n.groups, function(x) {
    which(comm.groups == x)
  })
  names(Gs) <- paste("G", 1:n.groups, sep = "")
  dist.matrix <- as.matrix(phylo.comp.dist)
  
  # FIX: add drop = FALSE to preserve matrix when group has 1 site
  PGall <- lapply(Gs, function(x) {
    dist.matrix[x, x, drop = FALSE]
  })
  
  PGall_similarity <- lapply(PGall, function(x) 1 - x)
  
  afilliation_by_grp <- lapply(PGall_similarity, function(x) {
    afilliation_by_grp <- matrix(NA, nrow(x), 2, dimnames = list(rownames(x),
                                                                 c("afilliation", "group")))
    for (z in 1:nrow(x)) {
      # FIX: handle single-site groups explicitly
      if (nrow(x) == 1) {
        afilliation_by_grp[z, 1] <- 0
      } else {
        dis <- as.data.frame(x[z, ])[-z, ]
        afilliation_by_grp[z, 1] <- mean(dis)
      }
    }
    return(afilliation_by_grp)
  })
  
  list_afilliation_by_grp <- vector(mode = "list", length = length(afilliation_by_grp))
  for (l in 1:length(afilliation_by_grp)) {
    # FIX: skip rescaling for single-site groups to avoid NaN from zero-range
    if (nrow(afilliation_by_grp[[l]]) == 1) {
      afilliation_by_grp_pad <- afilliation_by_grp[[l]]
      afilliation_by_grp_pad[, 1] <- 0
    } else {
      afilliation_by_grp_pad <- scales::rescale(afilliation_by_grp[[l]],
                                                c(0, 1))
    }
    afilliation_by_grp_pad[, 2] <- l
    list_afilliation_by_grp[[l]] <- afilliation_by_grp_pad
  }
  
  matrix_afilliation <- do.call(rbind, list_afilliation_by_grp)
  matrix_afilliation_org <- matrix_afilliation[match(rownames(dist.matrix),
                                                     rownames(matrix_afilliation)), ]
  return(matrix_afilliation_org)
}



# Function to get optimal K via KGS
get_optimal_k <- function(dist_obj) {
  methods <- c("ward.D", "ward.D2", "single", "complete", "average", "mcquitty")
  # Find best method by cophenetic correlation
  best_method <- methods[which.max(map_dbl(methods, ~cor(cophenetic(hclust(dist_obj, method = .x)), dist_obj)))]
  hc <- hclust(dist_obj, method = best_method)
  kgs_res <- maptree::kgs(hc, dist_obj)
  k <- as.numeric(names(which.min(kgs_res)))
  return(list(k = k, method = best_method))
}

# Function to fix rownames in community matrix
fix_mtn_names <- function(mat, fixes) {
  for (old in names(fixes)) {
    rownames(mat)[rownames(mat) == old] <- fixes[old]
  }
  return(mat)
}

library(vegan)        # For NMDS, PERMANOVA, betadisper
library(fpc)          # For clusterboot
library(indicspecies) # For Indicator Species Analysis (ISA)

##============================================================================##
## Validate Bioregions (NMDS, PERMANOVA, Clusterboot) ----
##============================================================================##
validate_bioregion <- function(dist_mat, clusters, hclust_method = "average", k) {
  
  # Failsafe: If your bio_result didn't output a method (it is NULL), default to "average"
  if (is.null(hclust_method) || length(hclust_method) == 0) {
    hclust_method <- "average"
  }
  
  # Ensure the input is strictly a distance object
  if(!inherits(dist_mat, "dist")) {
    dist_mat <- as.dist(dist_mat)
  }
  
  # 1. NMDS (Non-metric Multidimensional Scaling)
  nmds <- metaMDS(dist_mat, k = 2, trace = 0, trymax = 100)
  stress_val <- nmds$stress
  
  # Prepare NMDS dataframe for plotting
  nmds_df <- data.frame(
    NMDS1 = nmds$points[, 1],
    NMDS2 = nmds$points[, 2],
    Cluster = as.factor(clusters),
    Mountain = names(clusters)
  )
  
  # 2. Homogeneity of Dispersions
  # suppressWarnings hides the expected "negative squared distances" non-Euclidean warning
  suppressWarnings({
    disp <- betadisper(dist_mat, as.factor(clusters))
  })
  perm_disp <- permutest(disp, permutations = 999)
  
  # 3. PERMANOVA
  permanova <- adonis2(dist_mat ~ as.factor(clusters), permutations = 999)
  
  # 4. Cluster Stability (Bootstrap Resampling)
  set.seed(123)
  cboot <- fpc::clusterboot(
    data = dist_mat, 
    distances = TRUE, 
    B = 999, 
    bootmethod = "boot",
    clustermethod = fpc::disthclustCBI,
    method = hclust_method,  # Safely passes the clustering algorithm (e.g., "average")
    k = k,
    showplots = FALSE
  )
  
  stability_scores <- cboot$bootmean
  
  # Return everything nicely
  list(
    nmds_df      = nmds_df,
    nmds_stress  = stress_val,
    dispersion   = perm_disp,
    permanova    = permanova,
    stability    = stability_scores
  )
}
##============================================================================##


##============================================================================##
## Indicator Species Analysis (ISA) using Pearson's Phi ----
##============================================================================##
run_isa <- function(comm_mat, clusters) {
  
  # func = "r.g" calculates the group-size corrected Pearson's phi coefficient 
  # exactly as requested in the paper methodology.
  isa_result <- multipatt(
    x = comm_mat, 
    cluster = as.factor(clusters), 
    func = "r.g", 
    control = how(nperm = 999)
  )
  
  # Extract the summary of significant indicator species (p < 0.05)
  return(isa_result)
}
##============================================================================##


##============================================================================##
## Nearest distance excluding overlapping polygons ----
##============================================================================##
nearest_dist_no_overlap <- function(query_sf, target_sf) {
  # For each feature in query_sf, find the distance to the nearest
  
  # feature in target_sf that does NOT intersect it.
  out <- numeric(nrow(query_sf))
  
  for (i in seq_len(nrow(query_sf))) {
    overlaps    <- st_intersects(query_sf[i, ], target_sf, sparse = FALSE)[1, ]
    candidates  <- target_sf[!overlaps, ]
    
    if (nrow(candidates) > 0) {
      idx    <- st_nearest_feature(query_sf[i, ], candidates)
      out[i] <- as.numeric(st_distance(query_sf[i, ], candidates[idx, ]))
    } else {
      out[i] <- NA
    }
  }
  return(out)
}
##============================================================================##


##============================================================================##
## Establish the AIC function ----
##============================================================================##
AICFxn <- function(model){
  mod <- glm((1 - model$observed) ~ model$ecological,
           family = binomial(link = log))
  
  k <- length(which(model$coefficients > 0)) + 1 
  
  # Number of coefficients, plus 1 for the model intercept
  AIC <- (2*k) - (2*logLik(mod))
  
  dev<-((mod$null.deviance-mod$deviance)/mod$null.deviance)*100
  
  return(list(AIC, dev))
  } # end AICFxn
##============================================================================##


##============================================================================##
## gdm.varImp_v2 ----
##============================================================================##


gdm.varImp_v2 <- function (spTable, geo, splines = NULL, knots = NULL, predSelect = FALSE, 
                           nPerm = 50, pValue = 0.05, parallel = FALSE, cores = 2, 
                           sampleSites = 1, sampleSitePairs = 1, outFile = NULL) 
{
  # 1. Ensure necessary packages are loaded
  require(gdm)
  require(methods)
  require(pbapply)
  if (parallel) {
    require(parallel)
    require(doParallel)
    require(foreach)
  }
  
  # 2. Extract hidden 'gdm' functions safely into the local environment
  permutateSitePair <- utils::getFromNamespace("permutateSitePair", "gdm")
  subsample.sitepair <- utils::getFromNamespace("subsample.sitepair", "gdm")
  
  k <- NULL
  if (!is(spTable, "gdmData")) {
    warning("The spTable object is not of class 'gdmData'. See the formatsitepair function for help.")
  }
  if (!(is(spTable, "gdmData") | is(spTable, "matrix") | is(spTable, "data.frame"))) {
    stop("spTable argument needs to be of class 'gdmData', 'matrix', or 'data frame'")
  }
  if (ncol(spTable) < 6) {
    stop("spTable object requires at least 6 columns: distance, weights, s1.xCoord, s1.yCoord, s2.xCoord, s2.yCoord")
  }
  if (nrow(spTable) < 1) {
    stop("The spTable object contains zero rows of data.")
  }
  if (!(geo == TRUE | geo == FALSE)) {
    stop("The geo argument must be either TRUE or FALSE.")
  }
  if (is.null(splines) == FALSE & !is(splines, "numeric")) {
    stop("The splines argument needs to be a numeric data type.")
  }
  if (is.null(knots) == FALSE & !is(knots, "numeric")) {
    stop("The knots argument needs to be a numeric data type.")
  }
  if (!(predSelect == TRUE | predSelect == FALSE)) {
    stop("The predSelect argument must be either TRUE or FALSE.")
  }
  if ((is.null(nPerm) == FALSE & is.numeric(nPerm) == FALSE) | nPerm < 1) {
    stop("The nPerm argument needs to be a positive integer.")
  }
  if (!(parallel == TRUE | parallel == FALSE)) {
    stop("The parallel argument must be either TRUE or FALSE.")
  }
  if (parallel == TRUE & is.null(cores) == TRUE) {
    stop("If parallel==TRUE, the number of cores must be specified.")
  }
  if ((is.null(cores) == FALSE & is.numeric(cores) == FALSE) | cores < 1) {
    stop("The cores argument needs to be a positive integer.")
  }
  if (is.numeric(sampleSites) == FALSE | sampleSites < 0 | sampleSites > 1) {
    stop("The sampleSites argument needs to be a positive number between 0 and 1.")
  }
  if (is.numeric(sampleSitePairs) == FALSE | sampleSitePairs < 0 | sampleSitePairs > 1) {
    stop("The sampleSitePairs argument needs to be a positive number between 0 and 1.")
  }
  if (sampleSites == 0) {
    stop("A sampleSites value of 0 will remove all sites from the analysis.")
  }
  if (sampleSitePairs == 0) {
    stop("A sampleSitePairs value of 0 will remove all sites from the analysis.")
  }
  if (is.null(outFile) == FALSE) {
    if (is.character(outFile) == FALSE) {
      stop("The outFile argument needs to be a character string of the directory and file name you wish the tables to be written to")
    }
    outFileChar <- nchar(outFile)
    if (substr(outFile, outFileChar - 5, outFileChar) != ".RData") {
      outFile <- paste(outFile, ".RData", sep = "")
    }
    if (length(strsplit(outFile, "/")[[1]]) > 1) {
      splitOutFile <- strsplit(outFile, "/")[[1]][-length(strsplit(outFile, "/")[[1]])]
      dir.create(paste(splitOutFile, collapse = "/"))
    } else {
      outFile <- paste("./", outFile, sep = "")
    }
  }
  nPerm <- as.integer(nPerm)
  cores <- as.integer(cores)
  if (sampleSites < 1) {
    spTable <- subsample.sitepair(spTable, sampleSites = sampleSites)
    if (sampleSitePairs < 1) {
      warning("You have selected to randomly remove sites and/or site-pairs.")
    }
  }
  if (sampleSitePairs < 1) {
    numRm <- sample(1:nrow(spTable), round(nrow(spTable) * (1 - sampleSitePairs)))
    spTable <- spTable[-c(numRm), ]
  }
  rtmp <- spTable[, 1]
  if (length(rtmp[rtmp < 0]) > 0) {
    stop("The spTable contains negative distance values. Must be between 0 - 1.")
  }
  if (length(rtmp[rtmp > 1]) > 0) {
    stop("The spTable contains distance values greater than 1. Must be between 0 - 1.")
  }
  nVars <- (ncol(spTable) - 6)/2
  varNames <- colnames(spTable[c(7:(6 + nVars))])
  varNames <- sapply(strsplit(varNames, "s1\\."), "[[", 2)
  if (geo == TRUE) {
    nVars <- nVars + 1
    varNames <- c("Geographic", varNames)
  }
  if (nVars < 2) {
    stop("Function requires at least two predictor variables.")
  }
  message(paste0("Fitting initial model with all ", nVars, " predictors..."))
  Sys.sleep(0.5)
  fullGDM <- gdm(spTable, geo = geo, splines = splines, knots = knots)
  thiscoeff <- 1
  thisquant <- 1
  sumCoeff <- NULL
  for (i in 1:length(fullGDM$predictors)) {
    numsplines <- fullGDM$splines[[i]]
    holdCoeff <- NULL
    for (j in 1:numsplines) {
      holdCoeff[j] <- fullGDM$coefficients[[thiscoeff]]
      thiscoeff <- thiscoeff + 1
    }
    sumCoeff[i] <- sum(holdCoeff)
  }
  zeroSum <- fullGDM$predictors[which(sumCoeff == 0)]
  if (length(zeroSum) > 0) {
    for (p in 1:length(zeroSum)) {
      message(paste0("Sum of I-spline coefficients for predictor ", zeroSum[p], " = 0"))
      Sys.sleep(0.5)
    }
    for (p in 1:length(zeroSum)) {
      if (zeroSum[p] == "Geographic") {
        message("Setting Geo=FALSE and proceeding with permutation testing...")
      } else {
        message(paste0("Removing ", zeroSum[p], " and proceeding with permutation testing..."))
      }
      Sys.sleep(0.5)
    }
    if (length(grep("Geographic", zeroSum)) == 1) {
      geo <- FALSE
      zeroSum <- zeroSum[-grep("Geographic", zeroSum)]
    }
    for (z in zeroSum) {
      testVarCols1 <- grep(paste("^s1\\.", z, "$", sep = ""), colnames(spTable))
      testVarCols2 <- grep(paste("^s2\\.", z, "$", sep = ""), colnames(spTable))
      spTable <- spTable[, -c(testVarCols1, testVarCols2)]
    }
  }
  nVars <- (ncol(spTable) - 6)/2
  varNames <- colnames(spTable[c(7:(6 + nVars))])
  varNames <- sapply(strsplit(varNames, "s1\\."), "[[", 2)
  if (geo == TRUE) {
    nVars <- nVars + 1
    varNames <- c("Geographic", varNames)
  }
  if (nVars < 2) {
    stop("Function requires at least two predictor variables.")
  }
  if (cores > nVars) {
    cores <- nVars
  }
  splines <- rep(unique(splines), nVars)
  sortMatX <- sapply(1:nrow(spTable), function(i, spTab) { c(spTab[i, 3], spTab[i, 5]) }, spTab = spTable)
  sortMatY <- sapply(1:nrow(spTable), function(i, spTab) { c(spTab[i, 4], spTab[i, 6]) }, spTab = spTable)
  sortMatNum <- sapply(1:nrow(spTable), function(i) { c(1, 2) })
  sortMatRow <- sapply(1:nrow(spTable), function(i) { c(i, i) })
  fullSortMat <- cbind(as.vector(sortMatX), as.vector(sortMatY), as.vector(sortMatNum), as.vector(sortMatRow), rep(NA, length(sortMatX)))
  siteByCoords <- as.data.frame(unique(fullSortMat[, 1:2]))
  numSites <- nrow(siteByCoords)
  for (i in 1:numSites) {
    fullSortMat[which(fullSortMat[, 1] == siteByCoords[i, 1] & fullSortMat[, 2] == siteByCoords[i, 2]), 5] <- i
  }
  indexTab <- matrix(NA, nrow(spTable), 2)
  for (iRow in 1:nrow(fullSortMat)) {
    indexTab[fullSortMat[iRow, 4], fullSortMat[iRow, 3]] <- fullSortMat[iRow, 5]
  }
  rm(fullSortMat)
  rm(sortMatX)
  rm(sortMatY)
  rm(sortMatNum)
  rm(sortMatRow)
  rm(siteByCoords)
  exBySite <- lapply(1:numSites, function(i, index, tab) {
    rowSites <- which(index[, 1] %in% i)
    if (length(rowSites) < 1) {
      rowSites <- which(index[, 2] %in% i)
    }
    exSiteData <- tab[rowSites[1], ]
    return(exSiteData)
  }, index = indexTab, tab = spTable)
  outSite <- which(!(1:numSites %in% indexTab[, 1]))
  for (i in 1:length(exBySite)) {
    siteRow <- exBySite[[i]]
    if (i %in% outSite) {
      siteRow <- siteRow[grep("s2\\.", colnames(siteRow))]
      colnames(siteRow) <- sapply(strsplit(colnames(siteRow), "s2\\."), "[[", 2)
    } else {
      siteRow <- siteRow[grep("s1\\.", colnames(siteRow))]
      colnames(siteRow) <- sapply(strsplit(colnames(siteRow), "s1\\."), "[[", 2)
    }
    exBySite[[i]] <- siteRow
  }
  siteData <- do.call("rbind", exBySite)
  modelTestValues <- matrix(NA, 4, nVars * 5, dimnames = list(c("Model deviance", "Percent deviance explained", "Model p-value", "Fitted permutations"), c("All predictors", "1-removed", paste(seq(2, nVars * 5 - 1), "-removed", sep = ""))))
  varImpTable <- matrix(NA, nVars, nVars * 5 - 1)
  rownames(varImpTable) <- varNames
  colnames(varImpTable) <- c("All predictors", "1-removed", paste(seq(2, nVars * 5 - 2), "-removed", sep = ""))
  pValues <- nModsConverge <- varImpTable
  currSitePair <- spTable
  nullGDMFullFit <- 0
  message(paste0("Creating ", nPerm, " permuted site-pair tables..."))
  
  if (parallel == F | nPerm <= 25) {
    permSpt <- pbreplicate(nPerm, list(permutateSitePair(currSitePair, siteData, indexTab, varNames)))
  }
  
  if (parallel == T & nPerm > 25) {
    cl <- makeCluster(cores)
    registerDoParallel(cl)
    # Notice we removed "gdm:::" entirely from the variables and added exact explicit exports
    permSpt <- foreach(k = 1:nPerm, .verbose = F, .packages = c("gdm"), 
                       .export = c("permutateSitePair", "currSitePair", "siteData", "indexTab", "varNames")) %dopar% {
                         permutateSitePair(currSitePair, siteData, indexTab, varNames)
                       }
    stopCluster(cl)
  }
  varNames.x <- varNames
  message("Starting model assessment...")
  for (v in 1:length(varNames)) {
    
    if (length(varNames.x) < 2) {
      message("Only one predictor remains...variable assessment stopped.")
      message("Returning results evaluated up to this point.")
      v <- v - 1 
      break      
    }
    
    if (v > 1) {
      message(paste0("Removing ", names(elimVar), " and proceeding with the next round of permutations."))
    }
    if (is.numeric(splines)) {
      splines <- rep(unique(splines), length(varNames.x))
    }
    fullGDM <- gdm(currSitePair, geo = geo, splines = splines, knots = knots)
    message(paste0("Percent deviance explained by the full model =  ", round(fullGDM$explained, 3)))
    if (is.null(fullGDM) == TRUE) {
      warning(paste("The model did not converge when testing variable: ", varNames.x[v], ". Terminating analysis and returning output completed up to this point.", sep = ""))
      v <- v - 1 
      break
    }
    message("Fitting GDMs to the permuted site-pair tables...")
    permGDM <- lapply(permSpt, function(x) {
      suppressWarnings(try(gdm(x, geo = geo, splines = splines, knots = knots), silent = TRUE))
    })
    
    permModelDev <- sapply(permGDM, function(mod) {
      if(inherits(mod, "try-error") || is.null(mod)) return(NULL)
      mod$gdmdeviance
    })
    modPerms <- length(which(sapply(permModelDev, is.null) == TRUE))
    if (modPerms > 0) {
      permModelDev <- unlist(permModelDev[-(which(sapply(permModelDev, is.null) == T))])
    }
    
    if (v > 1) {
      colnames(modelTestValues)[v] <- colnames(pValues)[v] <- colnames(varImpTable)[v] <- colnames(nModsConverge)[v] <- paste0(names(elimVar), "-removed")
    }
    modelTestValues[1, v] <- round(fullGDM$gdmdeviance, 3)
    modelTestValues[2, v] <- round(fullGDM$explained, 3)
    modelTestValues[3, v] <- round(sum(permModelDev <= fullGDM$gdmdeviance)/(nPerm - modPerms), 3)
    modelTestValues[4, v] <- round(nPerm - modPerms, 0)
    
    if (parallel == TRUE) {
      if (length(varNames.x) < cores) {
        cores <- length(varNames.x)
      }
      cl <- makeCluster(cores)
      registerDoParallel(cl)
      # Adding completely explicit exports for parallel workers
      permVarDev <- foreach(k = 1:length(varNames.x), .verbose = F, .packages = c("gdm"), 
                            .export = c("currSitePair", "permSpt", "varNames.x", "geo", "splines", "knots")) %dopar% {
                              if (varNames.x[k] != "Geographic") {
                                lll <- lapply(permSpt, function(x, spt = currSitePair) {
                                  idx1 <- grep(paste("^s1\\.", varNames.x[k], "$", sep = ""), colnames(x))
                                  idx2 <- grep(paste("^s2\\.", varNames.x[k], "$", sep = ""), colnames(x))
                                  spt[, c(idx1, idx2)] <- x[, c(idx1, idx2)]
                                  return(spt)
                                })
                              }
                              if (varNames.x[k] == "Geographic") {
                                lll <- lapply(permSpt, function(x, spt = currSitePair) {
                                  s1 <- sample(1:nrow(spt), nrow(spt))
                                  s2 <- sample(1:nrow(spt), nrow(spt))
                                  s3 <- sample(1:nrow(spt), nrow(spt))
                                  s4 <- sample(1:nrow(spt), nrow(spt))
                                  spt[, 3] <- spt[s1, 3]
                                  spt[, 4] <- spt[s2, 4]
                                  spt[, 5] <- spt[s3, 5]
                                  spt[, 6] <- spt[s4, 6]
                                  return(spt)
                                })
                              }
                              gdmPermVar <- lapply(lll, function(x) {
                                suppressWarnings(try(gdm(x, geo = geo, splines = splines, knots = knots), silent=TRUE))
                              })
                              permModelDev <- sapply(gdmPermVar, function(mod) {
                                if(inherits(mod, "try-error") || is.null(mod)) return(NULL)
                                mod$gdmdeviance
                              })
                              return(permModelDev)
                            }
      stopCluster(cl)
    }
    
    if (parallel == FALSE) {
      permVarDev <- list()
      for (k in 1:length(varNames.x)) {
        if (varNames.x[k] != "Geographic") {
          message(paste0("Assessing importance of ", varNames.x[k], "..."))
          lll <- lapply(permSpt, function(x, spt = currSitePair) {
            idx1 <- grep(paste("^s1\\.", varNames.x[k], "$", sep = ""), colnames(x))
            idx2 <- grep(paste("^s2\\.", varNames.x[k], "$", sep = ""), colnames(x))
            spt[, c(idx1, idx2)] <- x[, c(idx1, idx2)]
            return(spt)
          })
        }
        if (varNames.x[k] == "Geographic") {
          message("Assessing importance of geographic distance...")
          lll <- lapply(permSpt, function(x, spt = currSitePair) {
            s1 <- sample(1:nrow(spt), nrow(spt))
            s2 <- sample(1:nrow(spt), nrow(spt))
            s3 <- sample(1:nrow(spt), nrow(spt))
            s4 <- sample(1:nrow(spt), nrow(spt))
            spt[, 3] <- spt[s1, 3]
            spt[, 4] <- spt[s2, 4]
            spt[, 5] <- spt[s3, 5]
            spt[, 6] <- spt[s4, 6]
            return(spt)
          })
        }
        gdmPermVar <- lapply(lll, function(x) {
          suppressWarnings(try(gdm(x, geo = geo, splines = splines, knots = knots), silent=TRUE))
        })
        permVarDev[[k]] <- unlist(sapply(gdmPermVar, function(mod) {
          if(inherits(mod, "try-error") || is.null(mod)) return(NULL)
          mod$gdmdeviance
        }))
      }
    }
    names(permVarDev) <- varNames.x
    nullDev <- fullGDM$nulldeviance
    for (var in varNames.x) {
      grepper <- grep(paste0("^", var, "$"), names(permVarDev))
      varDevTab <- unlist(permVarDev[[grepper]])
      nConv <- length(varDevTab)
      nModsConverge[which(rownames(varImpTable) == var), v] <- nConv
      varDevExplained <- 100 * (nullDev - varDevTab)/nullDev
      varImpTable[which(rownames(varImpTable) == var), v] <- median(100 * abs((varDevExplained - fullGDM$explained)/fullGDM$explained), na.rm=TRUE)
      if (var != "Geographic") {
        testVarCols1 <- grep(paste("^s1\\.", var, "$", sep = ""), colnames(currSitePair))
        testVarCols2 <- grep(paste("^s2\\.", var, "$", sep = ""), colnames(currSitePair))
        testSitePair <- currSitePair[, -c(testVarCols1, testVarCols2)]
        noVarGDM <- gdm(testSitePair, geo = geo, splines = splines[-1], knots = knots)
      } else {
        noVarGDM <- gdm(currSitePair, geo = F, splines = splines[-1], knots = knots)
      }
      permDevReduct <- noVarGDM$gdmdeviance - varDevTab
      pValues[which(rownames(pValues) == var), v] <- sum(permDevReduct >= (varDevTab - fullGDM$gdmdeviance))/(nConv)
    }
    if (max(na.omit(pValues[, v])) < pValue) {
      message("All remaining predictors are significant, ceasing assessment.")
      message(paste0("Percent deviance explained by final model = ", round(fullGDM$explained, 3)))
      message("Final set of predictors returned: ")
      for (vvv in 1:length(fullGDM$predictors)) {
        message(fullGDM$predictors[vvv])
      }
      break
    }
    if (predSelect == T) {
      elimVar <- which.max(pValues[, v])
      if (names(elimVar) != "Geographic") {
        remVarCols1 <- grep(paste("^s1\\.", names(elimVar), "$", sep = ""), colnames(currSitePair))
        remVarCols2 <- grep(paste("^s2\\.", names(elimVar), "$", sep = ""), colnames(currSitePair))
        currSitePair <- currSitePair[, -c(remVarCols1, remVarCols2)]
        permSpt <- lapply(permSpt, function(x) {
          x[, -c(remVarCols1, remVarCols2)]
        })
      } else {
        geo <- F
      }
      varNames.x <- varNames.x[-which(varNames.x == names(elimVar))]
    }
    if (v == 1 & predSelect == F) {
      message("Backwards elimination not selected by user (predSelect=F). Ceasing assessment.")
      message(paste0("Percent deviance explained by final model = ", round(fullGDM$explained, 3)))
      message("Final set of predictors returned: ")
      for (vvv in 1:length(fullGDM$predictors)) {
        message(fullGDM$predictors[vvv])
      }
      break
    }
  }
  
  if (v == 1 & predSelect == F) {
    modelTestVals <- data.frame(matrix(round(modelTestValues[, 1], 3), ncol = 1))
    rownames(modelTestVals) <- rownames(modelTestValues)
    colnames(modelTestVals) <- "All predictors"
    varImpTab <- data.frame(matrix(round(varImpTable[, 1], 3), ncol = 1))
    rownames(varImpTab) <- rownames(varImpTable)
    colnames(varImpTab) <- "All predictors"
    pVals <- varImpTab
    pVals[, 1] <- round(pValues[, 1], 3)
    nModsConv <- varImpTab
    nModsConv[, 1] <- round(nModsConverge[, 1], 3)
    outObject <- list(modelTestVals, varImpTab, pVals, nModsConv)
    names(outObject) <- c("Model assessment", "Predictor Importance", "Predictor p-values", "Model Convergence")
  } else {
    outObject <- list(round(modelTestValues[, 1:v, drop = FALSE], 3), 
                      round(varImpTable[, 1:v, drop = FALSE], 3), 
                      round(pValues[, 1:v, drop = FALSE], 3), 
                      nModsConverge[, 1:v, drop = FALSE])
    names(outObject) <- c("Model assessment", "Predictor Importance", "Predictor p-values", "Model Convergence")
  }
  
  if (is.null(outFile) == FALSE) {
    save(outObject, file = outFile)
  }
  return(outObject)
}
##============================================================================##



##============================================================================##
## 1. COLOUR PALETTE — SINGLE SWITCH POINT                                    ##
##                                                                            ##
## To restyle every panel at once, change PALETTE_SCHEME below.               ##
## Supported: "npg", "lancet", "jama", "nejm", "custom_muted"                ##
## To add MetBrewer / khroma / scico, uncomment the relevant library above    ##
## and add a case to the switch() inside get_palette().                       ##
##============================================================================##

get_palette <- function(k,
                        scheme = c("npg", "lancet", "jama",
                                   "nejm", "custom_muted")) {
  
  scheme <- match.arg(scheme)
  
  pool <- switch(scheme,
                 npg          = ggsci::pal_npg("nrc")(10),
                 lancet       = ggsci::pal_lancet("lanonc")(9),
                 jama         = ggsci::pal_jama("default")(7),
                 nejm         = ggsci::pal_nejm("default")(8),
                 custom_muted = c("#3B4992", "#EE0000", "#008B45", "#631879",
                                  "#008280", "#BB5500", "#5F559B", "#A20056")
                 ## ── add more schemes here ──
                 # met_hiroshige = MetBrewer::met.brewer("Hiroshige", k),
                 # khroma_bright = khroma::colour("bright")(k),
                 # scico_batlow  = scico::scico(k, palette = "batlow"),
  )
  
  if (k > length(pool)) {
    warning("k (", k, ") exceeds palette size (", length(pool),
            "); colours will recycle.", call. = FALSE)
    pool <- rep_len(pool, k)
  }
  
  pal <- pool[seq_len(k)]
  names(pal) <- as.character(seq_len(k))
  pal
}

# ── SET THE ACTIVE SCHEME ────────────────────────────────────────────────────
PALETTE_SCHEME <- "npg"
##============================================================================##


##============================================================================##
## 2. JOURNAL THEME                                                           ##
##============================================================================##

theme_journal <- function(base_size = 8.5) {
  theme_classic(base_size = base_size) %+replace%
    theme(
      text             = element_text(colour = "grey15"),
      axis.text        = element_text(colour = "grey25", size = rel(0.9)),
      axis.title       = element_text(size = rel(1.05)),
      axis.ticks       = element_line(colour = "grey50", linewidth = 0.25),
      axis.line        = element_line(colour = "grey50", linewidth = 0.3),
      plot.tag          = element_text(face = "bold", size = rel(1.4)),
      plot.margin       = margin(4, 6, 4, 4),
      legend.title      = element_text(face = "bold", size = rel(0.95)),
      legend.text       = element_text(size = rel(0.85)),
      legend.key.size   = unit(0.4, "cm"),
      panel.grid        = element_blank(),
      strip.background  = element_blank(),
      strip.text        = element_text(face = "bold", size = rel(0.95))
    )
}
##============================================================================##


##============================================================================##
## 3. ANALYSIS HELPERS                                                        ##
##============================================================================##

# ── 3a. Cophenetic evaluation across every standard linkage method ───────────

evaluate_linkage <- function(dist_mat,
                             methods = c("ward.D", "ward.D2", "single",
                                         "complete", "average", "mcquitty",
                                         "median", "centroid")) {
  tibble(method = methods) %>%
    mutate(
      cophenetic_r = map_dbl(method, function(m) {
        hc <- hclust(dist_mat, method = m)
        cor(cophenetic(hc), dist_mat)
      })
    ) %>%
    arrange(desc(cophenetic_r))
}
##============================================================================##


##============================================================================##
# ── 3b. Optimal K via KGS penalty ───────────────────────────────────────────
##============================================================================##

find_optimal_k <- function(hc, dist_mat) {
  kgs_vals <- maptree::kgs(hc, dist_mat)
  kgs_df   <- tibble(k       = as.integer(names(kgs_vals)),
                     penalty = as.numeric(kgs_vals))
  optimal  <- kgs_df %>% slice_min(penalty, n = 1, with_ties = FALSE)
  
  list(kgs_df = kgs_df, optimal_k = optimal$k)
}
##============================================================================##


##============================================================================##
# ── 3c. Build dendrogram aesthetics aligned to cluster colours ───────────────
##============================================================================##

build_dend_data <- function(hc, k, pal) {
  dend <- as.dendrogram(hc) %>%
    color_branches(k = k, col = pal)
  
  ggd        <- as.ggdend(dend)
  leaf_cols  <- get_leaves_branches_col(dend)
  leaf_names <- labels(dend)
  leaf_clust <- cutree(hc, k = k)[leaf_names]
  
  pal_aligned <- tapply(leaf_cols, leaf_clust, `[`, 1)
  
  ggd$labels$cluster <- factor(
    cutree(hc, k = k)[as.character(ggd$labels$label)],
    levels = names(pal_aligned)
  )
  
  list(dend = dend, ggd = ggd, pal_aligned = pal_aligned,
       leaf_order = leaf_names)
}
##============================================================================##


##============================================================================##
# ── 3d. Transform SES matrix into a valid (non-negative) dissimilarity ──────
##============================================================================##

ses_to_dist <- function(ses_dist) {
  mat <- as.matrix(ses_dist)
  shifted <- mat - min(mat, na.rm = TRUE)
  shifted <- (shifted + t(shifted)) / 2
  diag(shifted) <- 0
  as.dist(shifted)
}
##============================================================================##


##============================================================================##
# ── 3e. Full clustering pipeline (returns everything needed for plotting) ────
##============================================================================##

run_bioregion_pipeline <- function(dist_mat, shape_sf,
                                   id_col         = "Mountain",
                                   link_method    = "average",
                                   palette_scheme = PALETTE_SCHEME) {
  
  ids <- as.character(shape_sf[[id_col]])
  attr(dist_mat, "Labels") <- ids
  
  # Linkage evaluation
  linkage_eval <- evaluate_linkage(dist_mat)
  message("── Cophenetic correlations ──")
  print(linkage_eval)
  
  # Clustering
  hc     <- hclust(dist_mat, method = link_method)
  r_coph <- cor(dist_mat, cophenetic(hc))
  
  # Optimal K
  kgs_out   <- find_optimal_k(hc, dist_mat)
  optimal_k <- kgs_out$optimal_k
  kgs_df    <- kgs_out$kgs_df
  message("Optimal K = ", optimal_k)
  
  # Cut tree & silhouette
  clusters <- cutree(hc, k = optimal_k)
  sil      <- silhouette(clusters, dist_mat)
  
  # Palette & dendrogram data
  pal      <- get_palette(optimal_k, scheme = palette_scheme)
  dend_dat <- build_dend_data(hc, optimal_k, pal)
  
  # Spatial join
  group_df <- tibble(
    idcell  = names(clusters),
    cluster = factor(clusters),
    color   = dend_dat$pal_aligned[as.character(clusters)]
  )
  bioregion_sf <- shape_sf %>%
    mutate(idcell = as.character(.data[[id_col]])) %>%
    left_join(group_df, by = "idcell")
  
  list(
    hc           = hc,
    r_coph       = r_coph,
    optimal_k    = optimal_k,
    clusters     = clusters,
    sil          = sil,
    kgs_df       = kgs_df,
    pal          = pal,
    pal_aligned  = dend_dat$pal_aligned,
    dend         = dend_dat$dend,
    ggd          = dend_dat$ggd,
    leaf_order   = dend_dat$leaf_order,
    bioregion_sf = bioregion_sf,
    linkage_eval = linkage_eval
  )
}
##============================================================================##


##============================================================================##
## 4. PLOTTING FUNCTIONS                                                      ##
##============================================================================##

# ── 4a. Dendrogram ──────────────────────────────────────────────────────────

plot_dendrogram <- function(res,
                            y_label = expression(
                              beta[sim] ~ "dissimilarity (phylogenetic turnover)"
                            ),
                            tag = "(a)") {
  
  ggd         <- res$ggd
  pal_aligned <- res$pal_aligned
  max_h       <- max(ggd$segments$y, na.rm = TRUE)
  
  # Nudge labels slightly to the left of 0
  label_nudge <- -max_h * 0.02
  
  ggplot() +
    geom_segment(
      data      = ggd$segments,
      aes(x = y, y = x, xend = yend, yend = xend),
      colour    = ggd$segments$col,
      linewidth = 0.50,
      lineend   = "square"
    ) +
    geom_point(
      data = ggd$labels,
      aes(x = 0, y = x, colour = cluster),
      size = 1.3
    ) +
    geom_text(
      data = ggd$labels,
      aes(x = label_nudge, y = x, label = label, colour = cluster),
      angle       = 0,
      hjust       = 1,
      vjust       = 0.35,
      size        = 4, # Adjusted for better legibility
      fontface    = "italic",
      show.legend = FALSE
    ) +
    scale_colour_manual(
      values = pal_aligned,
      labels = paste("Region", names(pal_aligned)),
      name   = "Biogeographical\nregion"
    ) +
    # FIX: Increase 'mult' values to add space on left (0.45) and right (0.1)
    scale_x_continuous(
      name   = y_label,
      expand = expansion(mult = c(0.45, 0.1)) 
    ) +
    scale_y_continuous(expand = expansion(add = 0.6)) +
    annotate(
      "text", x = Inf, y = Inf,
      label = paste0("italic(r)[coph]==",
                     formatC(res$r_coph, format = "f", digits = 3)),
      parse = TRUE, hjust = 1.05, vjust = 1.5,
      size = 3, colour = "grey35"
    ) +
    guides(colour = guide_legend(override.aes = list(size = 2.5))) +
    # FIX: Ensure everything is visible even if it crosses the axis line
    coord_cartesian(clip = "off") + 
    theme_journal() +
    theme(
      axis.text.y     = element_blank(),
      axis.ticks.y    = element_blank(),
      axis.title.y    = element_blank(),
      axis.line.y     = element_blank(),
      # FIX: Add left margin padding to prevent labels from touching the edge
      plot.margin     = margin(10, 15, 10, 25), 
      legend.position = "none"
    ) +
    labs(tag = tag)
}
##============================================================================##


##============================================================================##
# ── 4b. KGS penalty ─────────────────────────────────────────────────────────
##============================================================================##

plot_kgs <- function(res, tag = "(b)") {
  
  kgs_df    <- res$kgs_df
  optimal_k <- res$optimal_k
  kgs_df$optimal <- kgs_df$k == optimal_k
  
  accent <- unname(res$pal_aligned[[1]])
  
  ggplot(kgs_df, aes(k, penalty)) +
    geom_line(colour = "grey50", linewidth = 0.50) +
    geom_point(
      aes(fill = optimal),
      shape = 21, size = 2.5, stroke = 0.35, colour = "grey30"
    ) +
    scale_fill_manual(
      values = c("TRUE" = accent, "FALSE" = "white"),
      guide  = "none"
    ) +
    scale_x_continuous(breaks = kgs_df$k) +
    annotate(
      "text",
      x     = optimal_k + 0.4,
      y     = min(kgs_df$penalty),
      label = paste0("italic(k)==", optimal_k),
      parse = TRUE, size = 2.5, colour = accent, hjust = 0
    ) +
    labs(
      x   = expression(italic(k) ~ "(clusters)"),
      y   = "KGS penalty",
      tag = tag
    ) +
    theme_journal()
}
##============================================================================##


##============================================================================##
# ── 4c. Silhouette ──────────────────────────────────────────────────────────
##============================================================================##

plot_silhouette <- function(res, tag = "(c)") {
  
  sil_mat <- res$sil[, 1:3]
  sil_df  <- tibble(
    cluster   = factor(sil_mat[, 1]),
    sil_width = sil_mat[, 3]
  ) %>%
    arrange(cluster, desc(sil_width)) %>%
    mutate(rank = row_number())
  
  mean_s <- mean(sil_df$sil_width)
  
  ggplot(sil_df, aes(rank, sil_width, fill = cluster)) +
    geom_col(width = 0.8, colour = NA) +
    geom_hline(
      yintercept = mean_s, linetype = "dashed",
      colour = "grey40", linewidth = 0.40
    ) +
    annotate(
      "text",
      x     = nrow(sil_df) * 0.75,
      y     = mean_s + 0.04,
      label = paste0("bar(s)==",
                     formatC(mean_s, format = "f", digits = 2)),
      parse = TRUE, size = 2.5, colour = "grey30"
    ) +
    scale_fill_manual(values = res$pal_aligned, guide = "none") +
    scale_y_continuous(
      limits = c(min(0, min(sil_df$sil_width) - 0.05), 1)
    ) +
    labs(
      x   = "Mountains (ordered by region)",
      y   = "Silhouette width",
      tag = tag
    ) +
    theme_journal() +
    theme(axis.text.x  = element_blank(),
          axis.ticks.x = element_blank())
}
##============================================================================##


##============================================================================##
# ── 4d. Spatial map ─────────────────────────────────────────────────────────
##============================================================================##

plot_map <- function(res, tag = "(d)") {
  
  ggplot() +
    geom_sf(
      data = res$bioregion_sf,
      aes(fill = color),
      colour = NA, linewidth = 0.1
    ) +
    scale_fill_identity() +
    theme_void() +
    theme(plot.margin = margin(10, 10, 10, 10)) +
    labs(tag = tag)
}
##============================================================================##


##============================================================================##
# ── 4e. Pairwise significance data + heatmap ────────────────────────────────
##============================================================================##

build_pairwise_df <- function(obs_dist, ses_dist, p_dist,
                              shape_sf, id_col = "Mountain",
                              alpha = 0.05) {
  
  mtn_names <- as.character(shape_sf[[id_col]])
  attr(obs_dist, "Labels") <- mtn_names
  attr(ses_dist, "Labels") <- mtn_names
  attr(p_dist,   "Labels") <- mtn_names
  
  obs_mat <- as.matrix(obs_dist)
  ses_mat <- as.matrix(ses_dist)
  p_mat   <- as.matrix(p_dist)
  
  # Build lower-triangle pairs only (avoid rowwise for speed)
  idx <- which(lower.tri(obs_mat), arr.ind = TRUE)
  tibble(
    Site1        = mtn_names[idx[, 1]],
    Site2        = mtn_names[idx[, 2]],
    obs_turnover = obs_mat[idx],
    SES          = ses_mat[idx],
    p_value      = p_mat[idx]
  ) %>%
    mutate(
      significant = p_value < alpha,
      direction   = case_when(
        SES > 0 & significant ~ "Higher than expected",
        SES < 0 & significant ~ "Lower than expected",
        TRUE                  ~ "Not significant"
      )
    )
}

plot_heatmap <- function(pairwise_df, dend_order, tag = "(e)") {
  
  ggplot(
    pairwise_df,
    aes(x = factor(Site1, levels = dend_order),
        y = factor(Site2, levels = dend_order))
  ) +
    geom_tile(aes(fill = obs_turnover),
              colour = "grey90", linewidth = 0.2) +
    geom_point(
      data  = filter(pairwise_df, significant),
      aes(shape = direction),
      size = 1.5, colour = "red"
    ) +
    scale_fill_viridis_c(
      name      = expression(beta[sim]),
      option    = "mako",
      direction = -1
    ) +
    scale_shape_manual(
      values = c("Higher than expected" = 24,
                 "Lower than expected"  = 25),
      name   = "Null model\n(p < 0.05)"
    ) +
    labs(x = NULL, y = NULL, tag = tag) +
    theme_journal() +
    theme(
      axis.text.x     = element_text(angle = 45, hjust = 1,
                                     face = "italic"),
      axis.text.y     = element_text(face = "italic"),
      legend.position = "right"
    ) +
    coord_fixed()
}
##============================================================================##


##============================================================================##
## 5. ASSEMBLE & EXPORT HELPER                                                ##
##============================================================================##
export_individual_panels <- function(panels, prefix = "results/") {
  walk(panels, function(spec) {
    ggsave(
      filename = paste0(prefix, spec$file),
      plot     = spec$plot,
      width    = spec$w,
      height   = spec$h,
      units    = "mm"
    )
  })
}
##============================================================================##


##============================================================================##
## MASTER ANALYSIS FUNCTION                                                   ##
##                                                                            ##
## Runs the full pipeline for one distance matrix:                            ##
##   clustering → dendrogram → KGS → silhouette → map                        ##
##   + optional SES clustering & pairwise heatmap                             ##
##   + composite figure assembly → export                                     ##
##============================================================================##

run_full_analysis <- function(config,
                              shape_sf       = mtn_alt_cropped,
                              palette_scheme = PALETTE_SCHEME,
                              output_dir     = "results/",
                              rds_dir        = "RDS/") {
  
  label <- config$label
  message("\n", strrep("=", 60))
  message("  RUNNING: ", label)
  message(strrep("=", 60))
  
  # ── 1. Observed clustering pipeline ──────────────────────────────────────
  res <- run_bioregion_pipeline(
    dist_mat       = config$dist_mat,
    shape_sf       = shape_sf,
    id_col         = "Mountain",
    link_method    = config$link_method %||% "average",
    palette_scheme = palette_scheme
  )
  
  saveRDS(res$bioregion_sf, file.path(rds_dir, config$rds_file))
  
  # ── 2. Core panels (always produced) ────────────────────────────────────
  p_dend <- plot_dendrogram(res, y_label = config$y_label, tag = "(a)")
  p_kgs  <- plot_kgs(res,       tag = "(b)")
  p_sil  <- plot_silhouette(res, tag = "(c)")
  p_map  <- plot_map(res,        tag = "(d)")
  
  # ── 3. Optional: SES clustering + pairwise heatmap ─────────────────────
  has_pairwise <- !is.null(config$ses_dist) && !is.null(config$p_dist)
  res_ses     <- NULL
  pairwise_df <- NULL
  p_heatmap   <- NULL
  
  if (has_pairwise) {
    
    ses_valid <- ses_to_dist(config$ses_dist)
    res_ses   <- run_bioregion_pipeline(
      dist_mat       = ses_valid,
      shape_sf       = shape_sf,
      id_col         = "Mountain",
      link_method    = config$link_method %||% "average",
      palette_scheme = palette_scheme
    )
    
    pairwise_df <- build_pairwise_df(
      obs_dist = config$dist_mat,
      ses_dist = config$ses_dist,
      p_dist   = config$p_dist,
      shape_sf = shape_sf,
      id_col   = "Mountain"
    )
    
    cat(sprintf("\n=== PAIRWISE SIGNIFICANCE: %s ===\n", label),
        "Total pairs:            ", nrow(pairwise_df), "\n",
        "Significant (p < 0.05): ", sum(pairwise_df$significant), "\n",
        "Higher than expected:   ",
        sum(pairwise_df$direction == "Higher than expected"), "\n",
        "Lower than expected:    ",
        sum(pairwise_df$direction == "Lower than expected"), "\n")
    
    p_heatmap <- plot_heatmap(pairwise_df,
                              dend_order = res$leaf_order,
                              tag = "(e)")
  }
  
  # ── 4. Composite figures ────────────────────────────────────────────────
  prefix <- config$file_prefix
  
  fig_4panel <- (p_dend | p_map) /
    (p_kgs  | p_sil) +
    plot_layout(heights = c(2, 1))
  
  ggsave(
    file.path(output_dir, paste0("Figure_Bioregionalisation_", prefix, ".png")),
    plot   = fig_4panel,
    width  = 170, height = 190, units = "mm"
  )
  
  if (has_pairwise) {
    fig_5panel <- (p_dend | p_map) /
      (p_kgs  | p_sil) /
      p_heatmap +
      plot_layout(heights = c(2, 1, 2))
    
    ggsave(
      file.path(output_dir, paste0("Figure_Bioregionalisation_Full_", prefix, ".png")),
      plot   = fig_5panel,
      width  = 170, height = 280, units = "mm"
    )
  }
  
  # ── 5. Individual panels ───────────────────────────────────────────────
  panels <- list(
    list(plot = p_dend, file = paste0("p_dend_", prefix, ".png"), w = 150, h = 180),
    list(plot = p_map,  file = paste0("p_map_",  prefix, ".png"), w = 85,  h = 112),
    list(plot = p_kgs,  file = paste0("p_kgs_",  prefix, ".png"), w = 85,  h = 56),
    list(plot = p_sil,  file = paste0("p_sil_",  prefix, ".png"), w = 85,  h = 56)
  )
  
  if (has_pairwise) {
    panels <- c(panels, list(
      list(plot = p_heatmap, file = paste0("p_heatmap_", prefix, ".png"),
           w = 170, h = 112)
    ))
  }
  
  export_individual_panels(panels, prefix = output_dir)
  
  message("Done: ", label, "\n")
  
  # Return everything invisibly for downstream inspection
  invisible(list(
    res_obs     = res,
    res_ses     = res_ses,
    pairwise_df = pairwise_df,
    panels      = list(dend = p_dend, kgs = p_kgs,
                       sil = p_sil, map = p_map,
                       heatmap = p_heatmap)
  ))
}
##============================================================================##


##============================================================================##
## Master function: cluster a distance matrix and return an sf with groups ----
##============================================================================##
build_bioregion <- function(dist_mat, shape_sf, method = "average") {
  
  ## 1. Hierarchical clustering
  hc <- hclust(dist_mat, method = method)
  
  ## 2. Find optimal k via KGS penalty
  kgs_res <- maptree::kgs(hc, dist_mat)
  optimal_k <- kgs_res %>%
    enframe() %>%
    arrange(value) %>%
    slice(1) %>%
    pull(name) %>%
    as.numeric()
  
  ## 3. Cut the tree
  clusters <- cutree(hc, k = optimal_k)
  
  ## 4. Cophenetic correlation
  r_coph <- cor(cophenetic(hc), dist_mat)
  
  ## 5. Silhouette
  sil <- silhouette(clusters, dist_mat)
  mean_sil <- mean(sil[, "sil_width"])
  
  ## 6. Assign clusters back to the spatial object
  group_df <- data.frame(
    Island = names(clusters),
    group    = as.integer(clusters),
    stringsAsFactors = FALSE
  )
  
  bio_sf <- shape_sf %>%
    left_join(group_df, by = "Island")
  
  ## 7. Return everything useful
  list(
    sf        = bio_sf,
    hclust    = hc,
    clusters  = clusters,
    optimal_k = optimal_k,
    r_coph    = r_coph,
    mean_sil  = mean_sil,
    sil       = sil,
    kgs_res   = kgs_res
  )
}
##============================================================================##


run_null_vmeasure <- function(sf_A, sf_B, x_name, y_name,
                              n_permutations = 999, B = 1,
                              precision = NULL, n_cores = 1) {
  
  # ── Required packages ─────────────────────────────────────────────────
  require(sf)
  require(sabre)
  require(future)
  require(future.apply)
  require(progressr)
  require(progress)        # backend for handler_progress (ETA support)
  
  # ── Step 1: Observed V-measure via vmeasure_calc ──────────────────────
  message("Step 1/4 \u2014 Computing observed V-measure...")
  obs_result <- do.call(sabre::vmeasure_calc, list(
    x = sf_A, y = sf_B,
    x_name = as.name(x_name), y_name = as.name(y_name),
    B = B, precision = precision
  ))
  obs_vm <- obs_result$v_measure
  message("  Observed V-measure: ", round(obs_vm, 4))
  
  # ── Step 2: Pre-compute spatial intersection (ONCE) ───────────────────
  message("Step 2/4 \u2014 Pre-computing spatial intersection...")
  
  sf_A_prep <- sf_A
  sf_A_prep$.a_label <- as.character(sf_A[[x_name]])
  
  sf_B_prep <- sf_B
  sf_B_prep$.b_idx <- seq_len(nrow(sf_B))
  
  if (!is.null(precision)) {
    sf_A_prep <- sf::st_set_precision(sf_A_prep, precision)
    sf_B_prep <- sf::st_set_precision(sf_B_prep, precision)
  }
  
  int_sf <- sf::st_intersection(
    sf_A_prep[, ".a_label"],
    sf_B_prep[, ".b_idx"]
  )
  int_sf <- sf::st_collection_extract(int_sf, type = "POLYGON")
  int_sf <- int_sf[!sf::st_is_empty(int_sf), ]
  
  int_areas <- as.numeric(sf::st_area(int_sf))
  int_data  <- sf::st_drop_geometry(int_sf)
  
  a_labels  <- int_data$.a_label
  b_indices <- int_data$.b_idx
  
  orig_b_labels <- as.character(sf_B[[y_name]])
  a_levels <- sort(unique(a_labels))
  b_levels <- sort(unique(orig_b_labels))
  
  message("  ", length(int_areas), " intersection pieces extracted.")
  
  # ── Free heavy spatial objects from memory ────────────────────────────
  rm(sf_A_prep, sf_B_prep, int_sf, int_data)
  sf_A <- NULL
  sf_B <- NULL
  gc(verbose = FALSE)
  
  # ── Step 3: Verify pre-computed approach matches vmeasure_calc ────────
  message("Step 3/4 \u2014 Verifying pre-computed approach...")
  verify_ct <- xtabs(w ~ a + b, data = data.frame(
    a = factor(a_labels, levels = a_levels),
    b = factor(orig_b_labels[b_indices], levels = b_levels),
    w = int_areas
  ))
  verify_vm <- sabre:::vmeasure(
    x = as.table(rowSums(verify_ct)),
    y = as.table(colSums(verify_ct)),
    z = verify_ct, B = B
  )$v_measure
  
  if (abs(verify_vm - obs_vm) > 0.01) {
    warning(
      "Pre-computed V-measure (", round(verify_vm, 4),
      ") differs from vmeasure_calc (", round(obs_vm, 4), ")."
    )
  } else {
    message("  Verification passed (pre-computed: ", round(verify_vm, 4),
            " vs vmeasure_calc: ", round(obs_vm, 4), ")")
  }
  rm(verify_ct)
  gc(verbose = FALSE)
  
  # ── Step 4: Run null permutations ─────────────────────────────────────
  # Memory estimate
  worker_bytes <- object.size(a_labels) + object.size(b_indices) +
    object.size(int_areas) + object.size(orig_b_labels) +
    object.size(a_levels) + object.size(b_levels)
  worker_mb <- as.numeric(worker_bytes) / 1024^2
  
  message("Step 4/4 \u2014 Running ", n_permutations, " permutations on ",
          n_cores, " core(s)")
  message(sprintf("  ~%.1f MB per worker \u00d7 %d cores = ~%.0f MB total",
                  worker_mb, n_cores, worker_mb * n_cores))
  
  # Set up future parallel backend
  if (n_cores > 1) {
    plan(multisession, workers = n_cores)
  } else {
    plan(sequential)
  }
  on.exit(plan(sequential), add = TRUE)
  
  # Configure progress bar with ETA
  handlers(handler_progress(
    format   = ":spin :current/:total [:bar] :percent | Elapsed: :elapsed | ETA: :eta",
    width    = 60,
    complete = "=",
    clear    = FALSE
  ))
  
  # Chunk size: aim for ~200 progress updates for smooth reporting
  chunk_n <- max(1L, n_permutations %/% 200L)
  
  null_vms <- with_progress({
    p <- progressor(steps = n_permutations)
    
    future_sapply(seq_len(n_permutations), function(i) {
      shuffled_b <- sample(orig_b_labels)
      b_perm     <- shuffled_b[b_indices]
      
      ct <- xtabs(w ~ a + b, data = data.frame(
        a = factor(a_labels, levels = a_levels),
        b = factor(b_perm, levels = b_levels),
        w = int_areas
      ))
      
      vm <- sabre:::vmeasure(
        x = as.table(rowSums(ct)),
        y = as.table(colSums(ct)),
        z = ct, B = B
      )$v_measure
      
      p()
      vm
    },
    future.seed       = TRUE,
    future.chunk.size = chunk_n,
    future.packages   = "sabre")
  })
  
  # ── Summary ───────────────────────────────────────────────────────────
  null_mean <- mean(null_vms)
  null_sd   <- sd(null_vms)
  ses       <- (obs_vm - null_mean) / null_sd
  p_val     <- (sum(null_vms >= obs_vm) + 1) / (n_permutations + 1)
  
  message("\nDone! P-value = ", format(p_val, digits = 4),
          " | SES = ", round(ses, 2))
  
  list(
    Observed          = obs_vm,
    Null_Mean         = null_mean,
    Null_SD           = null_sd,
    SES               = ses,
    P_value           = p_val,
    Null_Distribution = null_vms
  )
}

plot_null_vmeasure <- function(result, title = "V-measure Null Model") {
  require(ggplot2)
  
  df <- data.frame(vm = result$Null_Distribution)
  obs <- result$Observed
  
  # Place annotation on the less crowded side
  obs_side <- ifelse(obs > median(df$vm), "left", "right")
  hjust_val <- ifelse(obs_side == "left", 1.1, -0.1)
  
  ggplot(df, aes(x = vm)) +
    geom_histogram(aes(y = after_stat(density)),
                   bins = 40, fill = "#4393C3", color = "white", alpha = 0.8) +
    geom_density(linewidth = 0.8, color = "#2166AC") +
    geom_vline(xintercept = obs,
               color = "#D6604D", linewidth = 1.2, linetype = "dashed") +
    annotate("text", x = obs, y = Inf,
             label = paste0(
               "Observed = ", round(obs, 3),
               "\nSES = ", round(result$SES, 2),
               "\np = ", format(result$P_value, digits = 3)
             ),
             vjust = 1.5, hjust = hjust_val,
             color = "#B2182B", fontface = "bold", size = 4.2) +
    labs(
      x = "V-measure", y = "Density", title = title,
      subtitle = paste0(length(result$Null_Distribution), " permutations")
    ) +
    theme_minimal(base_size = 13) +
    theme(
      plot.title    = element_text(face = "bold"),
      plot.subtitle = element_text(color = "gray40"),
      panel.grid.minor = element_blank()
    )
}




## ── Convert a dist object to a sppTab for bioFormat = 3 ──────────────────
## The gdm package wants a square site-by-site dissimilarity data.frame
## with a numeric "Site" column.
make_sppTab <- function(dist_obj) {
  dist_obj %>%
    as.matrix() %>%                                            # dist -> matrix
    as.data.frame() %>%                                        # matrix -> df
    mutate(across(everything(), ~ replace(.x, is.nan(.x), 0))) %>%  # NaN -> 0
    mutate(Site = row_number()) %>%                            # add Site IDs
    select(Site, everything())
}


## ── Tidy the four pieces of a gdm.varImp result into 2 small tibbles ─────
summarise_varImp <- function(vi, label) {
  
  ## (a) Model-level metrics ("All predictors" column of "Model assessment")
  ma   <- vi[["Model assessment"]]
  full <- ma[, "All predictors"]
  
  model_tbl <- tibble(
    Matrix                 = label,
    Model_Deviance         = full["Model deviance"],
    Pct_Deviance_Explained = full["Percent deviance explained"],
    Model_P_Value          = full["Model p-value"],
    Fitted_Permutations    = full["Fitted permutations"]
  )
  
  ## (b) Predictor-level importance & p-values
  pi       <- vi[["Predictor Importance"]]
  pp       <- vi[["Predictor p-values"]]
  last_col <- colnames(pp)[ncol(pp)]              # last backward-elim. step
  
  pred_tbl <- tibble(
    Matrix                = label,
    Predictor             = rownames(pp),
    Importance_Full_Model = as.numeric(pi[, "All predictors"]),
    PValue_Full_Model     = as.numeric(pp[, "All predictors"]),
    PValue_Last_Elim_Step = as.numeric(pp[, last_col]),
    Last_Elimination_Step = last_col
  )
  
  list(model = model_tbl, predictors = pred_tbl)
}


## ── Identify predictors with p ≤ alpha at the LAST elimination step ──────
## Returns a list with:
##   $geo : TRUE/FALSE  -> was Geographic distance significant?
##   $env : character() -> environmental predictors that survived
##   $all : character() -> full set (geo + env)
get_sig_preds <- function(vi, alpha = 0.05) {
  pp       <- vi[["Predictor p-values"]]
  last_col <- pp[, ncol(pp)]
  names(last_col) <- rownames(pp)
  last_col <- last_col[!is.na(last_col)]            # drop eliminated preds
  sig      <- names(last_col[last_col <= alpha])
  geo_sig  <- "Geographic" %in% sig
  env_sig  <- setdiff(sig, "Geographic")
  list(geo = geo_sig, env = env_sig, all = sig)
}


## ── Compress cross-validation scalar metrics into a single-row tibble ────
compile_cv <- function(cv, label) {
  tibble(
    Matrix                   = label,
    Train_Deviance_Explained = cv$Train.Deviance.Explained,
    Test_Deviance_Explained  = cv$Test.Deviance.Explained,
    Mean_Error               = cv$Mean.Error,
    Mean_Absolute_Error      = cv$Mean.Absolute.Error,
    RMSE                     = cv$Root.Mean.Square.Error,
    Obs_Pred_Correlation     = cv$Obs.Pred.Correlation,
    Equalized_RMSE           = cv$Equalized.RMSE
  )
}


## ── Reshape isplineExtract() output to long format for ggplot ────────────
spline_to_long <- function(spl, label) {
  preds <- colnames(spl$x)
  map_dfr(preds, ~ tibble(
    Matrix    = label,
    Predictor = .x,
    x_value   = as.vector(spl$x[, .x]),
    y_value   = as.vector(spl$y[, .x])
  ))
}



