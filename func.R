#### ------------------------------------- ####
##     Functions
#### ------------------------------------- ####

#| @title PRRD
#| @description Calculate the PRRD
#| @param t1 Temperature 1
#| @param t2 Temperature 2
#| @param t1_c Mean for temperature 1 for constant group
#| @param t2_c Mean for temperature 2 for constant group
#| @param t1_f Mean for temperature 1 for fluctuating group
#| @param t2_f Mean for temperature 2 for fluctuating group
#| @param sd_t1_c Standard deviation for temperature 1 for constant group
#| @param sd_t2_c Standard deviation for temperature 2 for constant group
#| @param sd_t1_f Standard deviation for temperature 1 for fluctuating group
#| @param sd_t2_f Standard deviation for temperature 2 for fluctuating group
#| @param n_t1_c Sample size for temperature 1 for constant group
#| @param n_t2_c Sample size for temperature 2 for constant group
#| @param n_t1_f Sample size for temperature 1 for fluctuating group		
#| @param n_t2_f Sample size for temperature 2 for fluctuating group
#| @param type Type of output: 'ef' for effect size, 'v' for sampling variance
#| @return PRRD or its sampling variance
#| @examples
#| PRRD(t1 = 20, t2 = 30, t1_c = 5, t2_c = 10, t1_f = 5, t2_f = 10, 0.1, 0.2, 0.1, 0.2, 5, 10, 5, 10, type = 'ef') # Should be 0
#| PRRD(t1 = 20, t2 = 30, t1_c = 10, t2_c = 5, t1_f = 5, t2_f = 10, 0.1, 0.2, 0.1, 0.2, 5, 10, 5, 10, type = 'ef') # Should be negative because c is smaller than f
#| PRRD(t1 = 20, t2 = 30, t1_c = 5, t2_c = 10, t1_f = 10, t2_f = 5, 0.1, 0.2, 0.1, 0.2, 5, 10, 5, 10, type = 'ef') # Should be negative because c is bigger than f
#| PRRD(20, 30, 5, 10, 5, 10, 0.1, 0.2, 0.1, 0.2, 5, 10, 5, 10, type = 'v')
#| PRRD(20, 30, 5, 10, 5, 10, 3, 6, 2, 8, 2, 3, 2, 3, type = 'v')  # Should be bigger than l57
#|  check with single row, need to check direction etc
#| PRRD(t2 = 23.7, t1 = 17.7, 
#|   	 t1_c = 33.23927393, 
#| 	 t2_c = 17.81023102, 
#| 	 t1_f = 26.63861386, 
#| 	 t2_f = 16.82013201, 
#| 	 sd_t1_c = 2.25026431, 
#| 	 sd_t2_c = 1.666842873, 
#| 	 sd_t1_f = 2.833685746, 
#| 	 sd_t2_f = 2.833685746, 
#| 	 n_t1_c = 50, n_t2_c = 50, n_t1_f = 50, n_t2_f = 50, type = 'v')  

PRRD <- function(t1, t2, t1_c, t2_c, t1_f, t2_f, sd_t1_c, sd_t2_c, sd_t1_f, sd_t2_f, 
                 n_t1_c, n_t2_c, n_t1_f, n_t2_f, type = c('ef', 'v')) {
  
  type <- match.arg(type)
  
  # Check for invalid values
  if (any(c(t1_c, t2_c, t1_f, t2_f) <= 0)) {
    stop("All mean values must be positive to compute log response ratios.")
  }
  
  inv_temp_diff2 <- 1 / (t2 - t1)^2  # Precompute for efficiency
  
  # Effect size calculation
  if (type == 'ef') { 
    prrd <- ((log(t2_f) - log(t1_f)) - (log(t2_c) - log(t1_c))) / (t2 - t1)
    return(prrd)
  }
  
  # Sampling variance calculation
  if (type == 'v') { 
    v_prrd <- inv_temp_diff2 * (
      (sd_t1_c^2 / (n_t1_c * t1_c^2)) + 
      (sd_t2_c^2 / (n_t2_c * t2_c^2)) + 
      (sd_t1_f^2 / (n_t1_f * t1_f^2)) + 
      (sd_t2_f^2 / (n_t2_f * t2_f^2))
    )
    return(v_prrd)
  }
}

round_df <- function(df, digits = 2) {
  df[] <- lapply(df, function(col) {
    if (is.numeric(col)) round(col, digits) else col
  })
  return(df)
}


# Build table with results for supp
table_results <- function(model, group = NULL, study_name = "Study_ID", species_name = "Scientific_Name"){
  
  dat = model$data
  
  if(is.null(group)){
    tab <- dat %>% summarise(study = length(unique(!!sym(study_name))),
                             species = length(unique(!!sym(species_name))),
                             k = n())
  } else {
    tab <- dat %>% group_by(!!sym(group)) %>% summarise(study = length(unique(!!sym(study_name))),
                                                     species = length(unique(!!sym(species_name))),
                                                     k = n())
  }
  
  modelvals <- round_df(data.frame(
    "Estimate" = as.vector(model$b),
    "CI Low" = as.vector(model$ci.lb), 
    "CI High" = as.vector(model$ci.ub), 
    "df" = as.vector(model$ddf), 
    "p-value" = as.vector(model$pval), 
    check.names = FALSE, row.names = NULL), 3)
  
  return(cbind(tab, modelvals))
}

#' @title Plotting function for heterogeneity
#' @description Creates a ggplot object for visualizing heterogeneity in data.
#' @param data Data frame containing the data to be plotted.
#' @param x Column name for the x-axis.
#' @param y Column name for the y-axis.
#' @param type Type of plot (e.g., "Heterogeneity").
#' @param col Color for the bars in the plot.
#' @return A ggplot object.
#' 
het_plot <- function(data, x, y, type, col = wes_palette('GrandBudapest1', 4, type = 'discrete')[4]) {
  ggplot(data, aes(x = x, y = y)) +
        geom_col(alpha = 1, color = col, fill = col) +
        labs(y = expression(paste(italic(M)[])), x = "Level" , title = type) +
        scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) + 
        theme_bw() +
        theme(legend.background = element_blank(),
              axis.text = element_text(size = 12, color = "black"),
              axis.title = element_text(size = 12, color = "black"))
}


#' @title Tree checks and pruning
#' @description Checks consistency between a dataset and a phylogenetic tree,
#'              or prunes the tree to match the data.
#'
#' @param data A data frame containing species information.
#' @param tree A phylogenetic tree object (class `phylo`).
#' @param dataCol Column name or index in `data` that contains species names.
#' @param type Whether to run checks ("checks") or prune the tree ("prune").
#'
#' @details 
#' - "checks": returns counts of species in data vs tree, plus mismatched taxa.  
#' - "prune": drops species from the tree that are not in the data. Will error if 
#'            data contains species missing from the tree.
#'
#' @return
#' If `type = "checks"`, returns a list with:
#' - `SpeciesNumbers`: number of unique species in data and tree
#' - `Species_InTree_But_NotData`: taxa found in tree but missing in data
#' - `Species_InData_But_NotTree`: taxa found in data but missing in tree
#' 
#' If `type = "prune"`, returns a pruned `phylo` tree.
#'
tree_checks <- function(data, tree, dataCol, type = c("checks", "prune")){
  type = match.arg(type)
  # How many unique species exist in data and tree
  Numbers <- matrix(nrow = 2, ncol = 1)
  Numbers[1,1] <- length(unique(data[,dataCol])) 
  Numbers[2,1] <- length(tree$tip.label) 
  rownames(Numbers)<- c("Species in data:", "Species in tree:")
  # Missing species or species not spelt correct      
  species_list1= setdiff(sort(tree$tip.label), sort(unique(data[,dataCol])))
  species_list2= setdiff(sort(unique(data[,dataCol])), sort(tree$tip.label) )
  if(type == "checks"){
    return(list(SpeciesNumbers = data.frame(Numbers), 
                Species_InTree_But_NotData=species_list1, 
                Species_InData_But_NotTree=species_list2))
  }
  if(type == "prune"){
    if(length(species_list2) >=1) stop("Sorry, you can only prune a tree when you have no taxa existing in the data that are not in the tree")
    return(ape::drop.tip(tree, species_list1))
  }
}

#' @title Arcsine-square root transformation with SD adjustment
#' @description Applies an arcsine-square root transformation to a mean 
#'              (e.g. proportions) and adjusts the SD using the delta method.
#'
#' @param mean Mean value (typically a proportion between 0 and 1).
#' @param sd Standard deviation on the original scale.
#'
#' @return A list with:
#' - `mean`: transformed mean
#' - `sd`: transformed standard deviation
#'
asine_transform_with_sd <- function(mean, sd) {
  # Transform mean
  mean_t <- asin(sqrt(mean))
  # Delta-method SD on transformed scale. Equivalent to equation in paper
  sd_t   <- sd / (2 * sqrt(mean) * sqrt(1 - mean))
  list(mean = mean_t, sd = sd_t)
}

#' @title Back-transformed mean
#' @description Computes back-transformed mean values for arcsine-square root 
#'              or log-transformed data.
#'
#' @param m Mean value on the transformed scale.
#' @param s Standard deviation on the transformed scale.
#' @param per_transform Whether a proportion (arcsine) transformation was applied ("Yes"/"No").
#' @param ln_transform Whether a log transformation was applied ("Yes"/"No").
#'
#' @return A numeric value of the back-transformed mean.
#'
back_mean <- function(m, s, per_transform = "No", ln_transform = "No") {
  ifelse(per_transform == "Yes",
         asine_transform_with_sd(mean = m, sd = s)$mean,
         ifelse(ln_transform == "Yes",
                exp(m + (s^2)/2),
                m))
}

#' @title Back-transformed SD
#' @description Computes back-transformed standard deviations for arcsine-square root 
#'              or log-transformed data.
#'
#' @param m Mean value on the transformed scale.
#' @param s Standard deviation on the transformed scale.
#' @param per_transform Whether a proportion (arcsine) transformation was applied ("Yes"/"No").
#' @param ln_transform Whether a log transformation was applied ("Yes"/"No").
#'
#' @return A numeric value of the back-transformed standard deviation.
#'
back_sd <- function(m, s, per_transform = "No", ln_transform = "No") {
  ifelse(per_transform == "Yes",
         asine_transform_with_sd(mean = m, sd = s)$sd,
         ifelse(ln_transform == "Yes",
                sqrt((exp(s^2)-1) * exp(2*m + s^2)),
                s))
}
