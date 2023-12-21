#' To obtain confidence intervals for node-level network metrics
#'
#' @param network An igraph graph object consisting of observed network.
#' @param n_versions Number of bootstrapped versions to be used. (default = 100)
#' @param network_metrics_functions_list A list consisting of function definitions of the network metrics that the user wants to evaluate. Each element in the list should have an assigned name. Each function 
#' definition should include two parameters, one for the main network and another one for the subnetwork. See default example.
#' network_metrics_functions_list = c("degree" = igraph::degree, 
#'"strength" = igraph::strength , 
#'"betweenness" = igraph::betweenness, 
#'"clustering_coefficient" = function(x) \{
#'trans <- igraph::transitivity(x, type = "local", isolates = "zero")
#'names(trans) <- igraph::V(x)$name;return(trans)
#'\},
#'"eigenvector_centrality" = function(x) igraph::eigen_centrality(x)$vector
#')
#' @param n_cores Number of cores for parallel processing with default 1.
#'
#' @return A list of dataframes of class list_node_level_CI. Each element of list is a dataframe having five columns and
#'         having number of rows equal to number of nodes in the network. The five columns correspond to node_number,
#'         node_name, metric_value, lower_CI, upper_CI.
#'         correspond to subsampling_proportion and rows correspond to n_simulations.
#'         The entries of the matrix provide value of correlation between the nodes in 
#'         full network and the sub-sampled network for the corresponding metric. 
#' @export
#'
#' @examples
#'  \donttest{
#' data(elk_network_2010)
#' elk_node_level_CI <- obtain_node_level_CI(elk_network_2010)
#' plot(elk_node_level_CI)
#' }
obtain_node_level_CI <- function(network, 
                                 n_versions = 100,
                                 network_metrics_functions_list = c("degree" = igraph::degree, 
                                                                    "strength" = igraph::strength , 
                                                                    "betweenness" = igraph::betweenness, 
                                                                    "clustering_coefficient" = function(x){
                                                                      trans <- igraph::transitivity(x, type = "local", vids = igraph::V(x) ,isolates = "zero");
                                                                      names(trans) <- igraph::V(x)$name;
                                                                      return(trans)
                                                                    },
                                                                    "eigenvector_centrality" = function(x) igraph::eigen_centrality(x)$vector),
                                 n_cores = 1){
  
  bs_samples_matrix <- obtain_bootstrapped_samples(network, n_versions = n_versions)
  bs_samples_network <- sapply(bs_samples_matrix[[2]], function(x) igraph::graph_from_adjacency_matrix(x, mode = "undirected", weighted = TRUE))
  
  node_level_all_metrics <- parallel::mclapply(1:length(network_metrics_functions_list), function(j) {
    metric_CI <- calculate_CI_nodes(network, bs_samples_network, network_metrics_functions_list[[j]])
    #make a dataframe with values of node metric, names and index
    species_metric_CI <- data.frame("node_number" = 1:igraph::vcount(network), 
                                    "node_name" = igraph::V(network)$name, 
                                    "metric_value" = metric_CI$Observed_val,
                                    "lower_CI" = as.numeric(metric_CI$CI_lower),
                                    "upper_CI" = as.numeric(metric_CI$CI_upper))
    species_metric_CI <- species_metric_CI[order(species_metric_CI$metric_value, decreasing = TRUE), ]
    return(species_metric_CI)},
    mc.cores = n_cores)
  
  names(node_level_all_metrics) <- names(network_metrics_functions_list)
  
  class(node_level_all_metrics) =  "list_node_level_CI"
  
  return(node_level_all_metrics)
}

#' To plot the results for node-level confidence intervals
#'
#' @param x A list of dataframes obtained from obtain_node_level_CI function.
#' @param ... Further arguments are ignored.
#'
#' @return No return value, called for side effects. 
#' The plots show 95% confidence intervals along with the observed metric value for each of the nodes in the network.
#' @export
#' @method plot list_node_level_CI
#'
#'
#' @examples
#' \donttest{
#' data(elk_network_2010)
#' elk_node_level_CI <- obtain_node_level_CI(elk_network_2010)
#' plot(elk_node_level_CI)
#' }
plot.list_node_level_CI <- function(x,...){
  node_level_all_metrics <- x
  
  if(!inherits(node_level_all_metrics,"list_node_level_CI")){
    stop("List passed is not of class 'list_node_level_CI'")
  }
  
  for(j in 1:length(node_level_all_metrics)){
    if(max(node_level_all_metrics[[j]]$metric_value, node_level_all_metrics[[j]]$upper_CI) < 1){
      max_limit <- seq(0,max(node_level_all_metrics[[j]]$metric_value, node_level_all_metrics[[j]]$upper_CI), length.out = 5)
    }else max_limit <- 0:max(node_level_all_metrics[[j]]$metric_value, node_level_all_metrics[[j]]$upper_CI)
    
    # Create plotrix plot with confidence intervals
    plot(node_level_all_metrics[[j]]$metric_value, 
         xaxt="n", 
         yaxt ="n",
         type = 'n',
         xlab = "Node Number",
         ylab = names(node_level_all_metrics)[[j]],
         ylim = c(min(node_level_all_metrics[[j]]$metric_value, node_level_all_metrics[[j]]$lower_CI), max(node_level_all_metrics[[j]]$metric_value, node_level_all_metrics[[j]]$upper_CI)),
         xlim = c(1, length(node_level_all_metrics[[j]]$node_number)))
    
    graphics::abline(h =  pretty(max_limit), 
                     v = 1:length(node_level_all_metrics[[j]]$node_number), lty = 1, col = "grey89")
    
    plotrix::plotCI(x = 1:length(node_level_all_metrics[[j]]$node_number),
                    y = node_level_all_metrics[[j]]$metric_value,
                    li =node_level_all_metrics[[j]]$lower_CI,
                    ui = node_level_all_metrics[[j]]$upper_CI, 
                    scol = "black", 
                    add = TRUE)
    graphics::points(node_level_all_metrics[[j]]$metric_value, pch=19, col = "orange", cex = 0.7)
    graphics::axis(side = 2, at=pretty(max_limit), las=2)
    graphics::axis(side = 1, at=1:length(node_level_all_metrics[[j]]$node_number), 
                   labels=node_level_all_metrics[[j]]$node_number, las=2, tick = TRUE, cex.axis = 0.85)
  }
}




