#' Function to obtain permuted networks from raw datastream
#'
#' @param species_raw A dataframe consisting of raw GPS observations
#' @param temporal_thresh Temporal threshold in minutes
#' @param spatial_thresh Spatial threshold
#' @param n_permutations Number of permuted versions to obtain
#' @param n_cores Number of cores for parallel processing with default 1
#'
#' @return An object of class "list_permuted_networks" of size n_permutations where each element is a network of class igraph obtained by permuting raw datastream 
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
#' @param network_metrics_functions_list A list consisting of function definitions of the network metrics that the user wants to evaluate. Each element in the list should have an assigned name.
#'  Default = c("edge_density" = function(x) igraph::edge_density(x), "diameter" = function(x) igraph::diameter(x, weights = NA), "transitivity" = function(x) igraph::transitivity(x))
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
                                        network_metrics_functions_list = c("edge_density" = function(x) igraph::edge_density(x),
                                                                           "diameter" = function(x) igraph::diameter(x, weights = NA),
                                                                           "transitivity" = function(x) igraph::transitivity(x)),
                                        ...){
  networks_list <- x
  
  if(!inherits(networks_list, "list_permuted_networks")){
    stop("List passed is not of class 'list_permuted_networks'")
  }
  
  j <-1
  for(f in network_metrics_functions_list){
    permutations_metric <- unlist(lapply(networks_list, function(i) f(i)))
    graphics::hist(permutations_metric, 
         main = paste0(names(network_metrics_functions_list)[j]," distribution"), xlab = "Metric values (Permuted versions)")
    graphics::mtext(paste("Observed value : ", 
                          round(f(species_original_network),5), sep = ""), 
                    side = 3 , col = "red", cex = 0.8)
    j <- j +1
  }
}
