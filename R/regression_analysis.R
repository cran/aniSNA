#' To perform regression analysis for local network metrics
#'
#' @param network An igraph graph object consisting of observed network
#' @param n_simulations Number of sub-samples to be obtained at each level
#' @param subsampling_proportion A vector depicting proportions of sub-sampled nodes
#' @param network_metrics A vector depicting names of local global network metrics 
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
#' regression_slope_analyze(elk_network_2010)
#' }
regression_slope_analyze <- function(network, 
                                n_simulations = 10,
                                subsampling_proportion = c(0.1, 0.30, 0.50, 0.70, 0.90),
                                network_metrics = c("degree", "strength", "betweenness", "clustering_coefficient", "eigenvector_centrality")){
  
  regression_slope <- list()
  regression_slope <- lapply(1:length(network_metrics), function(i){ 
    regression_slope[[network_metrics[i]]] <- matrix(0, 
                                                       n_simulations, 
                                                       length(subsampling_proportion), 
                                                       dimnames = list(as.character(c(1:n_simulations)), as.character(subsampling_proportion*100)))
  })
  names(regression_slope) <- network_metrics
  
  for (i in 1:n_simulations) {
    for (j in 1:length(subsampling_proportion)) {
      random_sample_nodes <- as.vector(sample(igraph::V(network), size = subsampling_proportion[j] * igraph::gorder(network)))
      sub_network <- igraph::induced_subgraph(network, random_sample_nodes, impl = "auto")
      
      if("degree" %in% network_metrics){regression_slope$degree[i,j] <- stats::coef(stats::lm(igraph::degree(sub_network) ~ igraph::degree(network, igraph::V(sub_network)$name)))[2]}
      if("strength" %in% network_metrics){regression_slope$strength[i,j] <- stats::coef(stats::lm(igraph::strength(sub_network) ~ igraph::strength(network, igraph::V(sub_network)$name)))[2]}
      if("betweenness" %in% network_metrics){regression_slope$betweenness[i,j] <- stats::coef(stats::lm(igraph::betweenness(sub_network) ~ igraph::betweenness(network, igraph::V(sub_network)$name)))[2]}
      if("clustering_coefficient" %in% network_metrics){
        #Sometimes, clustering_coefficients can be computed for very low levels of sub sampling. The program will return NA if that is the case. 
        tryCatch(
          (regression_slope$clustering_coefficient[i,j] <- stats::coef(stats::lm(igraph::transitivity(sub_network, type = "local") ~ igraph::transitivity(network, type = "local", vids = igraph::V(sub_network)$name)))[2]),
          error = function(e){NA}
        )
      }
      if("eigenvector_centrality" %in% network_metrics){
        #scaling the eigenvector centrality
        eigen_cen_network <- igraph::eigen_centrality(network)$vector
        eigen_cen_network_scaled <- (eigen_cen_network - mean(eigen_cen_network))/(max(eigen_cen_network) - min(eigen_cen_network))
        regression_slope$eigenvector_centrality[i,j] <- stats::coef(stats::lm(as.vector(scale(igraph::eigen_centrality(sub_network)$vector, scale = FALSE)) ~ eigen_cen_network_scaled[igraph::V(sub_network)$name], singular.ok = TRUE))[2]}
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