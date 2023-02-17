
#' Calculates and prints network summary statistics
#'
#' @param network An undirected network with nodes representing animal IDs and edges representing associations between them.
#'
#' @return No return value, called for side effects. The function prints values of network metrics to the console. 
#' @export
#'
#' @examples
#' data(elk_network_2010)
#' get_network_summary(elk_network_2010)
get_network_summary <- function(network) {
  cat("The number of vertices are ", igraph::gorder(network), "\n") 
  cat("The number of edges are ", igraph::gsize(network), "\n") 
  cat("Vertex Attributes are : ", igraph::list.vertex.attributes(network), "\n")
  cat("Edge Attributes are : ", igraph::list.edge.attributes(network), "\n")
  cat("The edge density of the network is : ", igraph::edge_density(network), "\n")
  cat("The mean degree is ", mean(igraph::degree(network)), "\n")
  cat("The mean strength is ", mean(igraph::strength(network)), "\n")
  cat("The diameter is ", igraph::diameter(network, weights = NA), "\n")
  cat("The transitivity is ", igraph::transitivity(network), "\n")
  cat("The mean geodesic distance is ", igraph::mean_distance(network), "\n")
}


#' Visualize Animal Network
#'
#' @param species_network An igraph graph object consisting of observed network.
#' @param seed Seed to be set for layout.
#'
#' @return No return value, called for side effects. The plots depict a visualisation of network structure. 
#' @export
#'
#' @examples 
#' data(elk_network_2010)
#' plot_network(elk_network_2010)
plot_network <- function(species_network, seed = 1){
  set.seed(seed)
  layout <- igraph::layout.fruchterman.reingold(species_network)
  plot(species_network,
       vertex.size = 7, 
       vertex.frame.color = "black",
       vertex.label = NA, 
       edge.lty=c("solid"),
       layout = layout)
}