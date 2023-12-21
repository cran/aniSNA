#' To perform regression analysis for local network metrics
#'
#' @param network An igraph graph object consisting of observed network
#' @param n_simulations Number of sub-samples to be obtained at each level
#' @param subsampling_proportion A vector depicting proportions of sub-sampled nodes
#' @param network_metrics_functions_list A list consisting of function definitions of the network metrics that the user wants to evaluate. Each element in the list should have an assigned name. Each function 
#' definition should include two parameters, one for the main network and another one for the subnetwork. See default example.
#'  Default = c("degree" = function(net, sub_net) igraph::degree(net, v = igraph::V(sub_net)$name),
#'  "strength" = function(net, sub_net) igraph::strength(net, v = igraph::V(sub_net)$name),
#'  "betweenness" = function(net, sub_net) igraph::betweenness(net, v = igraph::V(sub_net)$name),
#'  "clustering_coefficient" = function(net, sub_net) igraph::transitivity(net, type = "local", vids = igraph::V(sub_net)$name),
#'  "eigenvector_centrality" = function(net, sub_net) igraph::eigen_centrality(net)$vector[igraph::V(sub_net)$name])
#'
#' @return A list of network metrics of class list_regression_matrices. Each element of list is a matrix whose columns 
#'         correspond to subsampling_proportion and rows correspond to n_simulations.
#'         The entries of the matrix provide value of the slope of regression when the 
#'         nodal values in sub-sampled network are regressed upon the values of the same 
#'         nodes in the full network for the corresponding metric.
#' @export
#'
#' @examples
#' \donttest{
#' data(elk_network_2010)
#' elk_regression_analysis <- regression_slope_analyze(elk_network_2010)
#' plot(elk_regression_analysis)
#' }
regression_slope_analyze <- function(network, 
                                     n_simulations = 10,
                                     subsampling_proportion = c(0.1, 0.30, 0.50, 0.70, 0.90),
                                     network_metrics_functions_list = c("degree" = function(net, sub_net) igraph::degree(net, v = igraph::V(sub_net)$name),
                                                                        "strength" = function(net, sub_net) igraph::strength(net, v = igraph::V(sub_net)$name),
                                                                        "betweenness" = function(net, sub_net) igraph::betweenness(net, v = igraph::V(sub_net)$name),
                                                                        "clustering_coefficient" = function(net, sub_net) igraph::transitivity(net, type = "local", vids = igraph::V(sub_net)$name),
                                                                        "eigenvector_centrality" = function(net, sub_net) igraph::eigen_centrality(net)$vector[igraph::V(sub_net)$name])){
  
  regression_slope <- list()
  regression_slope <- lapply(1:length(network_metrics_functions_list), function(i){ 
    regression_slope[[names(network_metrics_functions_list)[i]]] <- matrix(0, 
                                                                           n_simulations, 
                                                                           length(subsampling_proportion), 
                                                                           dimnames = list(as.character(c(1:n_simulations)), as.character(subsampling_proportion*100)))
  })
  names(regression_slope) <- names(network_metrics_functions_list)
  
  for (i in 1:n_simulations) {
    for (j in 1:length(subsampling_proportion)) {
      random_sample_nodes <- as.vector(sample(igraph::V(network), size = subsampling_proportion[j] * igraph::gorder(network)))
      sub_network <- igraph::induced_subgraph(network, random_sample_nodes, impl = "auto")
      
      k <- 1
      for(f in network_metrics_functions_list){
        tryCatch(
          (regression_slope[[names(network_metrics_functions_list)[k]]][i,j] <- stats::coef(stats::lm(f(sub_network, sub_network) ~ f(network, sub_network)))[2]),
          error = function(e){NA}
        )
        k <- k+1
      }
    }
  }
  
  class(regression_slope) <- "list_regression_matrices"
  
  return(regression_slope)
}




#' To plot regression analysis results
#'
#' @param x A list of matrices obtained from regression_slope_analyze function
#' @param ... Further arguments are ignored
#'
#' @return No return value, called for side effects. The plots show regression slope values corresponding to proportion of individuals in the sample.
#' @export
#' @method plot list_regression_matrices
#'
#' @examples
#' \donttest{
#' data(elk_network_2010)
#' elk_regression_analysis <- regression_slope_analyze(elk_network_2010)
#' plot(elk_regression_analysis)
#' }
plot.list_regression_matrices <- function(x,...){
  
  regression_results <- x
  
  if(!inherits(regression_results,"list_regression_matrices")){
    stop("List passed is not of class 'list_regression_matrices'")
  }
  
  
  for(i in 1:length(regression_results)){
    mean_slope = apply(regression_results[[i]], 2, mean, na.rm = TRUE)
    sample_proportions <- as.integer(colnames(regression_results[[1]]))/100
    
    #Remove the index where mean_slope is NA
    ind_remove <- which(is.na(mean_slope))
    if (length(ind_remove) > 0){
      mean_slope = mean_slope[-ind_remove]
      sample_proportions <- sample_proportions[-ind_remove]
      
    } else{
      sample_proportions <- sample_proportions
    }
    
    plot(sample_proportions, rep(1, length(sample_proportions)), 
         type = 'l',
         col = "black",
         lty = "dashed",
         main = paste(" Slope of Linear Regression for \n", stringr::str_to_sentence(names(regression_results), locale = "en")[i],sep = ""),
         ylim = c(min(0, mean_slope), 1),
         xlim = (c(0,1)),
         ylab = "Slope",
         xlab = "Proportion of Individuals in Subsample")
    
    graphics::lines(sample_proportions, sample_proportions, 
          col = 'black',
          lty = "dashed")
    
    graphics::lines(sample_proportions, mean_slope, lwd = 2,
          col = "red")
    
  }
}