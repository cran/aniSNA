#' To generate subsamples of the permuted networks and obtain network metrics of those subsamples
#'
#' @param networks_list A list of igraph objects obtained by permuting the observed network
#' @param subsampling_proportion A vector depicting the levels (in proportion) at which subsamples to be taken
#' @param network_metrics_functions_list A list consisting of function definitions of the network metrics that the user wants to evaluate. Each element in the list should have an assigned name.
#'  Default = c("edge_density" = function(x) igraph::edge_density(x), "diameter" = function(x) igraph::diameter(x, weights = NA), "transitivity" = function(x) igraph::transitivity(x))
#'
#' @return A list of network metrics of class "Subsampled_Permuted_Network_Metrics". Each element of list is a matrix whose columns 
#'         correspond to subsampling_proportion and rows correspond to the number of networks in networks_list. 
#'         The entries of the matrix provide values of the corresponding metric.
#' @export
#'
#' @examples
#' \donttest{
#' data(elk_2010_permutations)
#' elk_subsamples_permuted_networks <- subsamples_permuted_networks(elk_2010_permutations)
#' plot(elk_subsamples_permuted_networks, elk_network_2010)
#' }
subsamples_permuted_networks <- function(networks_list,
                                         subsampling_proportion = c(0.1, 0.30, 0.50, 0.70, 0.90),
                                         network_metrics_functions_list = c("edge_density" = function(x) igraph::edge_density(x),
                                                                            "diameter" = function(x) igraph::diameter(x, weights = NA),
                                                                            "transitivity" = function(x) igraph::transitivity(x))){
  
  subsampling_result <- list() 
  subsampling_result <- lapply(1:length(network_metrics_functions_list), function(i){ 
    subsampling_result[[names(network_metrics_functions_list)[i]]] <- matrix(0, length(networks_list), 
                                                                             length(subsampling_proportion), 
                                                                             dimnames = list(as.character(c(1:length(networks_list))), as.character(subsampling_proportion*100)))
  })
  
  names(subsampling_result) <- names(network_metrics_functions_list)
  
  for(i in 1:length(networks_list)){
    for(j in 1:length(subsampling_proportion)){
      
      random_sample_nodes <- as.vector(sample(igraph::V(networks_list[[i]]), size = subsampling_proportion[j] * igraph::gorder(networks_list[[i]])))
      sub_network <- igraph::induced_subgraph(networks_list[[i]], random_sample_nodes, impl = "auto")
      
      metric_values <- network_metrics_evaluate(sub_network, network_metrics_functions_list)
      
      for(net_met in 1:length(network_metrics_functions_list)){
        subsampling_result[[net_met]][i,j] <- metric_values[[net_met]]
      }
    }
  }
  class(subsampling_result) <- "Subsampled_Permuted_Network_Metrics"
  return(subsampling_result)
}





#' To plot sub-sampling results of the original network and permuted networks
#'
#' @param x A list of matrices obtained from subsamples_permuted_networks function of class "Subsampled_Permuted_Network_Metrics"
#' @param network An igraph graph object consisting of observed network
#' @param n_simulations For subsampling results of original network, this determines the number of sub-samples to be obtained at each level
#' @param subsampling_proportion A vector depicting the levels (in proportion) at which subsamples to be taken. This parameter should be the same as the
#' subsampling_proportion parameter passed for subsamples_permuted_networks to obtain x.
#' @param network_metrics_functions_list This is the same argument that is passed for obtaining the results from the function subsamples_permuted_networks.
#' A list consisting of function definitions of the network metrics that the user wants to evaluate. Each element in the list should have an assigned name. 
#' Default = c("edge_density" = function(x) igraph::edge_density(x), "diameter" = function(x) igraph::diameter(x, weights = NA), "transitivity" = function(x) igraph::transitivity(x))
#' @param ... Further arguments are ignored
#'
#' @return No return value, called for side effects. The boxplots show side-by-side comparison of network metrics distribution from subsamples of observed network and subsamples from permuted networks. 
#' 
#' 
#' @method plot Subsampled_Permuted_Network_Metrics
#' @export
#'
#' @examples
#' \donttest{
#' data(elk_2010_permutations, elk_network_2010)
#' elk_subsamples_permuted_networks <- subsamples_permuted_networks(elk_2010_permutations)
#' plot(elk_subsamples_permuted_networks, elk_network_2010)
#' }
plot.Subsampled_Permuted_Network_Metrics <- function(x, network,
                                                     n_simulations = 100, 
                                                     subsampling_proportion = c(0.1, 0.30, 0.50, 0.70, 0.90),
                                                     network_metrics_functions_list = c("edge_density" = function(x) igraph::edge_density(x),
                                                                                        "diameter" = function(x) igraph::diameter(x, weights = NA),
                                                                                        "transitivity" = function(x) igraph::transitivity(x)),...){
  
  
  species_permuted_results = x
  
  if(!inherits(species_permuted_results, "Subsampled_Permuted_Network_Metrics")){
    stop("x passed is not of class 'Subsampled_Permuted_Network_Metrics'")
  }
  
  variable <- value <- category <- NULL #To predefine variables to be used for ggplot2
  
  species_sample_results <- subsampled_network_metrics(network, n_simulations, subsampling_proportion, network_metrics_functions_list)
  
  species_permuted_results <- lapply(species_permuted_results, 
                                     function(matrix){cbind(matrix, rep("Permuted", nrow(matrix)))})
  species_sample_results <- lapply(species_sample_results, 
                                   function(matrix){cbind(matrix, rep("Observed", nrow(matrix)))})
  
  
  species_combine <- vector(mode = "list", length = length(species_sample_results))
  species_combine <- lapply(1: length(species_sample_results), 
                            function(i){rbind(species_permuted_results[[i]], species_sample_results[[i]])})
  names(species_combine) <- names(species_sample_results)
  species_combine <- lapply(species_combine, 
                            function(matrix){colnames(matrix) <- c(subsampling_proportion, "category")
                            matrix <- as.data.frame(matrix,stringsAsFactors = FALSE)
                            })
  
  species_long <- lapply(species_combine, function(df){df <- reshape::melt(df, id = "category")
  df$value <- as.numeric(df$value)
  return(df)})
  
  for(i in 1:length(species_permuted_results)){
    
    plot_metrics <- ggplot2::ggplot(species_long[[names(species_permuted_results)[i]]], ggplot2::aes(x=variable, y=value, fill=category)) +
      ggplot2::geom_boxplot()+ 
      ggplot2::geom_hline(yintercept =  network_metrics_functions_list[[names(species_permuted_results)[i]]](network), color = "red", linetype = "dashed") +
      ggplot2::theme(legend.position=c(0.8, 0.9), legend.title = ggplot2::element_blank(),
                     plot.title = ggplot2::element_text(size = 10,hjust = 0.5, face="bold"),
                     axis.title=ggplot2::element_text(size=9)) +
      ggplot2::ggtitle(names(species_permuted_results)[i])+
      ggplot2::ylab("Value") + ggplot2::xlab("Sub-sample size (in %)")+
      ggplot2::scale_fill_manual(values=c("orange", "skyblue"))
    
    print(plot_metrics)
    
  }              
}