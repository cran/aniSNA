#' Function to obtain permuted networks from raw datastream
#'
#' @param species_raw A dataframe consisting of raw GPS observations
#' @param temporal_thresh Temporal threshold in minutes
#' @param spatial_thresh Spatial threshold
#' @param n_permutations Number of permuted versions to obtain
#' @param n_cores Number of cores for parallel processing with default 1
#'
#' @return A list of size n_permutations where each element is a network of class igraph obtained by permuting raw datastream 
#' @export
#'
#' @examples
#' \donttest{
#' data(elk_data_2010)
#' permuted_versions <- obtain_permuted_network_versions(elk_data_2010, 
#' temporal_thresh = 7, spatial_thresh = 15, n_permutations = 10, n_cores = 2)
#' }
obtain_permuted_network_versions <- function(species_raw, temporal_thresh, spatial_thresh, n_permutations, n_cores = 1){
  
  #Obtain randomized order of dates of observation
  species_permutations <- obtain_permutation_dates(species_raw, n_permutations, n_cores)
  
  permuted_networks_list <- vector(mode = "list", length = length(species_permutations))
  
  pb <- utils::txtProgressBar(min = 0,      # Minimum value of the progress bar
                       max = length(species_permutations), # Maximum value of the progress bar
                       style = 3,    # Progress bar style (also available style = 1 and style = 2)
                       width = 50,   # Progress bar width. Defaults to getOption("width")
                       char = "=")   # Character used to create the bar
  
  for(i in 1:length(species_permutations)){
    species_raw$datetime <- species_permutations[[i]]
    species_raw$datetime <- as.POSIXct(species_raw$datetime, format = "%Y-%m-%d %H:%M")
    interactions <- get_interactions(species_raw = species_raw, temporal_thresh = temporal_thresh, spatial_thresh = spatial_thresh, n_cores = n_cores)
    permuted_networks_list[[i]] <- network_from_interactions(species_raw = species_raw, interactions = interactions, n_cores = n_cores)
    utils::setTxtProgressBar(pb, i)
  }
  close(pb)
  
  class(permuted_networks_list) <- "list_permuted_networks"
  
  return(permuted_networks_list)
}

#' Function to plot the network metrics distribution of permuted networks
#'
#' @param x A list of igraph objects obtained obtained using the function obtain_permuted_network_versions
#' @param species_original_network An igraph object which is the original network
#' @param network_metrics A vector depicting names of global network metrics. This should be supplied as a character vector and the values 
#' should be chosen from "mean_strength", "density", "diameter", "transitivity". (default = c("mean_strength", "density", "diameter", "transitivity")).
#'
#' @param ... Further arguments are ignored.
#'
#' @return No return value, called for side effects.
#' 
#' @method plot list_permuted_networks
#' 
#' @export
#'
#' @examples
#' \donttest{
#' data(elk_data_2010, elk_network_2010)
#' permuted_versions <- obtain_permuted_network_versions(elk_data_2010, 
#' temporal_thresh = 7, spatial_thresh = 15, n_permutations = 10, n_cores = 2)
#' plot(permuted_versions, elk_network_2010)
#' }
plot.list_permuted_networks <- function(x, 
                                     species_original_network,
                                     network_metrics = c("density", "mean_strength", "diameter", "transitivity"),
                                     ...){
  networks_list <- x
  
  if(!inherits(networks_list, "list_permuted_networks")){
    stop("List passed is not of class 'list_permuted_networks'")
  }
  
  if("density" %in% network_metrics){

    permutations_density <- unlist(lapply(networks_list, function(i) igraph::edge_density(i)))
    den_den <- stats::density(permutations_density, na.rm = TRUE)
    plot(den_den, 
         main = paste("Density Distribution of permuted versions", sep = ""))
    graphics::mtext(paste("Observed Value - ", 
                round(igraph::edge_density(species_original_network),5), sep = ""), 
          side = 3 , col = "red", cex = 0.8)
  }
  
  if("mean_strength" %in% network_metrics){
    
    permutations_mean_strength <- unlist(lapply(networks_list, function(i) mean(igraph::strength(i))))
    den_ms <- stats::density(permutations_mean_strength, na.rm = TRUE)
    plot(den_ms, 
         main = paste("Mean Strength Distribution of permuted versions", sep = ""))
    graphics::mtext(paste("Observed Value - ", 
                round(mean(igraph::strength(species_original_network)),5), sep = ""), 
          side = 3 , col = "red", cex = 0.8)
  }
  
  if ("transitivity" %in% network_metrics) {
    
    permutations_transitivity <- unlist(lapply(networks_list, function(i) igraph::transitivity(i)))
    den_tran <- stats::density(permutations_transitivity, na.rm = TRUE)
    plot(den_tran, 
         main = paste("Transitivity Distribution of permuted versions", sep = ""))
    graphics::mtext(paste("Observed Value - ", 
                round(igraph::transitivity(species_original_network),5), sep = ""), 
          side = 3 , col = "red", cex = 0.8)
  }
  
  if ("diameter" %in% network_metrics) {
    
    permutations_diameter <- unlist(lapply(networks_list, function(i) igraph::diameter(i)))
    den_diam <- stats::density(permutations_diameter, na.rm = TRUE)
    plot(den_diam, 
         main = paste("Diameter Distribution of permuted versions", sep = ""))
    graphics::mtext(paste("Observed Value - ", 
                round(igraph::diameter(species_original_network),5), sep = ""), 
          side = 3 , col = "red", cex = 0.8)
  }
  
}