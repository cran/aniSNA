calculate_CI_nodes <- function(network, bs_samples_network, network_metric){
  
  confidence_intervals_node_df <- data.frame("node_id" = character(igraph::vcount(network)), 
                                             "CI_lower" = integer(igraph::vcount(network)), 
                                             "CI_upper" = integer(igraph::vcount(network)))
  
  observed_val <- network_metrics_evaluate(network, list(network_metric))[[1]]
  metric_value <- lapply(bs_samples_network, network_metric)
  
  k <-1
  for(i in igraph::V(network)$name){
    node_metric_vector <- unlist(sapply(1:length(bs_samples_network), function(x){
      if(i %in% igraph::V(bs_samples_network[[x]])$name){ value_vector <- rep(metric_value[[x]][i], table(igraph::V(bs_samples_network[[x]])$name)[[i]]) 
      }else value_vector <- NA
      return(value_vector)
    }))
    #Obtain 95% confidence intervals
    quant <- stats::quantile(node_metric_vector, probs=c(0.25,0.75), na.rm = TRUE)
    confidence_intervals_node_df[k,] <- c(i, quant[1], quant[2])
    k <- k+1
  }
  
  confidence_intervals_node_df$Observed_val <- observed_val
  
  return(confidence_intervals_node_df)
}