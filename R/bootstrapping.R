
#' To obtain bootstrapped versions of a network's adjacency matrix
#'
#' @param network An igraph object
#' @param n_nodes Number of nodes to be selected in bootstrapped versions (default : All nodes)
#' @param n_versions Number of bootstrapped versions required
#' @param seed seed number
#'
#' @return A list of class bootstrapped_pvalue_matrix consisting of two elements. The first element contains the adjacency matrix of the original network 
#'         and the second element contains bootstrapped versions of the adjacency matrices.
#' @export
#'
#' @examples
#' data(elk_network_2010)
#' obtain_bootstrapped_samples(elk_network_2010, n_versions = 100)
obtain_bootstrapped_samples <- function(network, n_nodes = igraph::gorder(network), n_versions = 1000, seed = 12345){
  
  #Obtain a sub-network of order "size"
  if(n_nodes < igraph::gorder(network)){
    sample_nodes <- sample(igraph::V(network), size = n_nodes)
    network <- igraph::induced_subgraph(network, sample_nodes, impl = "auto")
  }
  
  # Obtain adjacency matrix from the network
  network_matrix <- igraph::as_adjacency_matrix(network, type = "both", sparse = FALSE, attr = "weight")
  bootstrapped_versions <- net_bootstrap(network_matrix, n_versions = n_versions) 
  
  return(bootstrapped_versions)
}


#' To obtain two non-overlapping bootstrapped versions and obtain p-values for the significance of difference between them
#'
#' @param network An igraph object
#' @param n_versions Number of bootstrapped versions to be used (default = 1000)
#' @param seed seed number 
#' @param n.iter Number of iterations at each level
#' @param network_metrics_functions_list A list consisting of function definitions of the network metrics that the user wants to evaluate. Each element in the list should have an assigned name.
#'  Default = c("edge_density" = function(x) igraph::edge_density(x), "diameter" = function(x) igraph::diameter(x, weights = NA), "transitivity" = function(x) igraph::transitivity(x))
#'
#' @return A matrix of p-values whose rows correspond to the sub-sample size and columns correspond to the chosen network metric. The sub-sample size values (corresponding to rows) 
#' occur in multiples of 5 and range from 5 to a maximum of half the number of nodes in the network
#' @export
#'
#' @examples
#' \donttest{
#' data(elk_network_2010)
#' mean_pvalue_matrix <- bootstrapped_difference_pvalues(elk_network_2010, n_versions = 100)
#' plot(mean_pvalue_matrix)
#' }
bootstrapped_difference_pvalues <- function(network, 
                                            n_versions = 1000, 
                                            seed = 12345, 
                                            n.iter = 10, 
                                            network_metrics_functions_list = c("edge_density" = function(x) igraph::edge_density(x),
                                                                               "diameter" = function(x) igraph::diameter(x, weights = NA),
                                                                               "transitivity" = function(x) igraph::transitivity(x))){
  
  subsample_size <- 5*(1:floor(floor(igraph::gorder(network)/2)/5))
  
  mean_metrics_pvalue <- matrix(NA, nrow = length(subsample_size), ncol = length(network_metrics_functions_list))
  
  for(i in 1:length(subsample_size)){
    metrics_pvalue <- p_value_matrix(network, size_subnet = subsample_size[i], n.iter, network_metrics_functions_list, n_versions)
    mean_value <- apply(metrics_pvalue, 2, mean, na.rm=TRUE)
    mean_metrics_pvalue[i,] <- mean_value
  }
  colnames(mean_metrics_pvalue) <- names(network_metrics_functions_list)
  rownames(mean_metrics_pvalue) <- as.character(5*(1:floor(floor(igraph::gorder(network)/2)/5)))
  
  class(mean_metrics_pvalue) = "bootstrapped_pvalue_matrix"
  
  return(mean_metrics_pvalue)
}

#' To plot the results obtained from bootstrapped_difference_pvalues function
#'
#' @param x A matrix of p-values obtained from bootstrapped_difference_pvalues function
#' @param ... Further arguments are ignored.
#'
#' @return No return value, called for side effects. The plot shows p-values between 0 and 1 corresponding to each sample size. 
#' @method plot bootstrapped_pvalue_matrix
#' @export
#'
#' @examples
#' \donttest{
#' data(elk_network_2010)
#' mean_pvalue_matrix <- bootstrapped_difference_pvalues(elk_network_2010, n_versions = 100)
#' plot(mean_pvalue_matrix)
#' }
plot.bootstrapped_pvalue_matrix <- function(x,...){
  
  bootstrapped_results <- x
  
  if(!inherits(bootstrapped_results,"bootstrapped_pvalue_matrix")){
    stop("Matrix passed is not of class 'bootstrapped_pvalue_matrix'")
  }
  
  subsample_size <- seq(5,5*nrow(bootstrapped_results), 5)
  col_vec <- c("red", "blue", "green", "yellow", "black")
  
  plot(subsample_size, 
         bootstrapped_results[,1], 
         ylim = c(0, 1),
         xlim = c(0, max(subsample_size)),
         type = "b",
         col = col_vec[1],
         pch = 16,
         xlab = "Sample Size",
         ylab = "p - value")
  
  if(ncol(bootstrapped_results) > 1){
    
    for(i in 2:ncol(bootstrapped_results)){
      graphics::points(subsample_size, 
             bootstrapped_results[,i],
             type = "b",
             col = col_vec[i],
             pch = 16)
    }
  }
  
  graphics::abline(h = 0.1, col = "brown")
  graphics::legend("topleft", 
                   legend = colnames(bootstrapped_results), 
                   col = col_vec[1:ncol(bootstrapped_results)],
                   pch = 16)
  
  
}


