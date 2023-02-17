#' To perform correlation analysis for local network metrics
#'
#' @param network An igraph graph object consisting of observed network
#' @param n_simulations Number of sub-samples to be obtained at each level
#' @param subsampling_proportion A vector depicting proportions of sub-sampled nodes
#' @param network_metrics A vector depicting names of local global network metrics 
#'
#' @return A list of network metrics of class list_correlation_matrices. Each element of list is a matrix whose columns 
#'         correspond to subsampling_proportion and rows correspond to n_simulations.
#'         The entries of the matrix provide value of correlation between the nodes in 
#'         full network and the sub-sampled network for the corresponding metric. 
#' @export
#'
#' @examples
#' \donttest{
#' data(elk_network_2010)
#' correlation_analyze(elk_network_2010)
#' }
correlation_analyze <- function(network, 
                               n_simulations = 10,
                               subsampling_proportion = c(0.1, 0.30, 0.50, 0.70, 0.90),
                               network_metrics = c("degree", "strength", "betweenness", "clustering_coefficient", "eigenvector_centrality")){
  
  correlation_values <- list()
  correlation_values <- lapply(1:length(network_metrics), function(i){ 
    correlation_values[[network_metrics[i]]] <- matrix(0, 
                                           n_simulations, 
                                           length(subsampling_proportion), 
                                           dimnames = list(as.character(c(1:n_simulations)), as.character(subsampling_proportion*100)))
  })
  names(correlation_values) <- network_metrics
  
  for (i in 1:n_simulations) {
    for (j in 1:length(subsampling_proportion)) {
      random_sample_nodes <- as.vector(sample(igraph::V(network), size = subsampling_proportion[j] * igraph::gorder(network)))
      sub_network <- igraph::induced_subgraph(network, random_sample_nodes, impl = "auto")
      
      if("degree" %in% network_metrics){correlation_values$degree[i,j] <- stats::cor(igraph::degree(network, igraph::V(sub_network)$name), igraph::degree(sub_network))}
      if("strength" %in% network_metrics){correlation_values$strength[i,j] <- stats::cor(igraph::strength(network, igraph::V(sub_network)$name), igraph::strength(sub_network))}
      if("betweenness" %in% network_metrics){correlation_values$betweenness[i,j] <- stats::cor(igraph::betweenness(network, igraph::V(sub_network)$name), igraph::betweenness(sub_network))}
      if("clustering_coefficient" %in% network_metrics){correlation_values$clustering_coefficient[i,j] <- stats::cor(igraph::transitivity(network, type = "local", vids = igraph::V(sub_network)$name), igraph::transitivity(sub_network, type = "local"), use = "pairwise.complete.obs")}
      if("eigenvector_centrality" %in% network_metrics){
        eigen_cen_whole <- igraph::eigen_centrality(network)$vector
        correlation_values$eigenvector_centrality[i,j] <- stats::cor(eigen_cen_whole[igraph::V(sub_network)$name], igraph::eigen_centrality(sub_network)$vector)}
    }
  }
  
  class(correlation_values) = "list_correlation_matrices"
  
  return(correlation_values)
}


#' To plot correlation analysis results
#'
#' @param x A list of matrices obtained from correlation_analyze function
#' @param ... Further arguments are ignored
#'
#' @return No return value, called for side effects. The plots show mean and standard deviation of correlation coefficients obtained over multiple iterations.
#' @export
#' @method plot list_correlation_matrices
#'
#' @examples
#' \donttest{
#' data(elk_network_2010)
#' elk_correlation_analysis <- correlation_analyze(elk_network_2010)
#' plot(elk_correlation_analysis)
#' }
plot.list_correlation_matrices <- function(x,...){
  
  correlation_results <- x
  
  if(!inherits(correlation_results,"list_correlation_matrices")){
    stop("List passed is not of class 'list_correlation_matrices'")
  }
  

  for(i in 1:length(correlation_results)){

    #To include mean and standard deviation in the plots
    mean_correlation = apply(correlation_results[[i]], 2, mean, na.rm = TRUE)
    sd_correlation = apply(correlation_results[[i]], 2, stats::sd, na.rm = TRUE)
    sample_proportions <- as.integer(colnames(correlation_results[[1]]))/100
    
    #Remove the index where mean or sd is NA
    ind_remove <- unique(c(which(is.na(mean_correlation)), which(is.na(sd_correlation))))
    
    if (length(ind_remove) > 0){
      mean_correlation = mean_correlation[-ind_remove]
      sd_correlation = sd_correlation[-ind_remove]
      sample_proportions <- sample_proportions[-ind_remove]
    } else{
      sample_proportions <- sample_proportions
    }

    plot(sample_proportions, mean_correlation,
         type = 'n',
         xlab = "Proportion of Individuals in Subsample",
         ylab = "Correlation",
         main = names(correlation_results)[i],
         ylim = c(min(0, mean_correlation - sd_correlation, na.rm = TRUE), 1))

    upper_limit <- mean_correlation + sd_correlation
    upper_limit[upper_limit > 1] <- 1

    lower_limit <- mean_correlation - sd_correlation
    lower_limit[lower_limit < -1] <- -1

    graphics::polygon(c(rev(sample_proportions), sample_proportions),
            c(rev(upper_limit), lower_limit),
            col = 'grey80', border = NA)
    graphics::lines(sample_proportions, mean_correlation,
          col = 'black')
    graphics::lines(c(rev(sample_proportions), sample_proportions),
          c(rev(upper_limit), lower_limit),
          lty = 'dashed',
          col = 'red')
  }
}


