

#' To obtain width of confidence intervals using bootstrapped versions at each level of sub-sampling
#'
#' @param network An igraph object
#' @param n_versions Number of bootstrapped versions to be used. (default = 100)
#' @param seed seed number
#' @param n.iter Number of iterations at each level. (default = 10)
#' @param network_metrics Network metrics to be evaluated. This should be supplied as a character vector and the values 
#' should be chosen from "mean_degree", "mean_strength", "density", "diameter", "transitivity". (default = c("mean_degree", "mean_strength", "density", "diameter", "transitivity"))
#' @param scaled_metrics Optional. A vector subset of network_metrics with the names of metrics that should be scaled. 
#' Values can be chosen from c("mean_degree", "mean_strength", "diameter").
#'
#'
#' @return A matrix of class Width_CI_matrix containing width of Confidence Intervals where each row corresponds to the sub-sample size and columns correspond to the chosen network metric.
#' @export
#'
#' @examples
#' \donttest{
#' data(elk_network_2010)
#' width_CI(elk_network_2010, n_versions = 100)
#' }
width_CI <- function(network,
                     n_versions = 100,
                     seed = 12345, 
                     n.iter = 10, 
                     network_metrics = c("mean_degree", "mean_strength", "density", "diameter", "transitivity"),
                     scaled_metrics = NULL){
  
  sample_size_values <- seq(10, igraph::gorder(network), 10)
  mean_value_CI_len <- data.frame(temp = numeric(0))
  for(i in 1:length(network_metrics)){
    mean_value_CI_len[[network_metrics[i]]] <- numeric(0)
  }
  mean_value_CI_len <- mean_value_CI_len[,-1]
  
  j <- 1
  for(s in sample_size_values){
    metrics_CI_len <- CI_matrix(network, size_subnet = s, n_versions, n.iter, network_metrics)
    mean_value_CI_len[j,] <- apply(metrics_CI_len, 2, mean, na.rm=TRUE)
    j <- j+1
  }
  mean_value_CI_len <- cbind(sample_size_values, mean_value_CI_len)
  
  if(is.null(scaled_metrics) == FALSE){
   
    if(all(scaled_metrics %in% network_metrics)){
      for(k in 1:length(scaled_metrics)){
        mean_value_CI_len[[paste0(scaled_metrics[k], "_", "scaled")]] <- mean_value_CI_len[[scaled_metrics[k]]]/mean_value_CI_len$sample_size_values
      }
    }else{
      stop("scaled_metrics is not a subset of network_metrics")
    } 
  }
  
  class(mean_value_CI_len) =  "Width_CI_matrix"
  
  return(mean_value_CI_len)
}


#' To plot the results obtained from width_CI function
#'
#' @param x A matrix of width of Confidence Intervals obtained from width_CI function
#' @param ... Further arguments are ignored.
#'
#' @return No return value, called for side effects. Plots show width of confidence intervals corresponding to number of individuals in the sub-sample.
#' @method plot Width_CI_matrix
#' @export
#'
#' @examples
#' \donttest{
#' data(elk_network_2010)
#' width_CI_elk <- width_CI(elk_network_2010, n_versions = 100)
#' plot(width_CI_elk)
#' }
plot.Width_CI_matrix <- function(x,...){
  
  width_CI_results <- x
  
  if(!inherits(width_CI_results,"Width_CI_matrix")){
    stop("Matrix passed is not of class 'Width_CI_matrix'")
  }
  
  width_CI_results <- as.data.frame(do.call(cbind, width_CI_results))
    
  for(i in 2:ncol(width_CI_results)){
      plot(width_CI_results$sample_size_values, 
           width_CI_results[,i], 
           type = "b",
           col = "black",
           main = colnames(width_CI_results)[i],
           pch = 16,
           xlab = "Sample Size",
           ylab = "Width of Confidence Interval")
  }
}
