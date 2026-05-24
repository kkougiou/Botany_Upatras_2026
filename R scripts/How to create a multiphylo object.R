##============================================================================##
## Load the libraries ----
##============================================================================##
library(rtrees)
library(betapart)
library(ape)
library(ggdendro)
library(dendextend)
library(vegan)
library(igraph)
library(bioregion)
library(rnetcarto)
library(maptree)
library(cluster)
library(tidyverse)
library(zoo) 
library(patchwork)
library(cluster)
library(phyloregion)
library(furrr)
library(future.apply)
library(progressr)
##============================================================================##


##============================================================================##
## Load the biotic data ----
##============================================================================##
plants <- readxl::read_excel("Excel/Kiklades.xlsx")
##============================================================================##


##============================================================================##
## The unique taxa for which to create the 100 phylogenies ----
##============================================================================##
taxa <- plants %>% 
  pivot_longer(
    cols = -Island,
    names_to = "Taxon", 
    values_to = "presence"
  ) %>% 
  distinct(Taxon)
##============================================================================##


##============================================================================##
## Pre-load the megatree ONCE ----
##============================================================================##
mega_tree <- megatrees::tree_plant_otl
##============================================================================##


##============================================================================##
## Set up parallel backend (furrr, which rtrees already depends on) ----
##============================================================================##
n_workers <- parallelly::availableCores() - 1
plan(multisession, workers = n_workers)
##============================================================================##


##============================================================================##
## Progress handler with ETA ----
##============================================================================##
handlers(handler_progress(
  format = ":spin [:bar] :current/:total (:percent) | Elapsed: :elapsed | ETA: :eta"
))

n_trees <- 100
##============================================================================##


##============================================================================##
## Generate trees in parallel with progress ----
##============================================================================##
trees_list <- with_progress({
  p <- progressor(steps = n_trees)
  
  future_map(seq_len(n_trees), function(i) {
    tree <- rtrees::get_one_tree(
      sp_list   = taxa,
      tree      = mega_tree,
      taxon     = "plant",
      scenario  = "random_below_basal",
      show_grafted = FALSE,
      .progress = "none"
    )
    p()
    tree
  }, .options = furrr_options(seed = TRUE))
})

names(trees_list) <- paste0("tree", seq_len(n_trees))
class(trees_list) <- "multiPhylo"
##============================================================================##


##============================================================================##
## Clean up ----
##============================================================================##
plan(sequential)
##============================================================================##


##============================================================================##
## Save the trees to disk ----
##============================================================================##
saveRDS(trees_list, "RDS/trees_list_in_islands.rds")
##============================================================================##