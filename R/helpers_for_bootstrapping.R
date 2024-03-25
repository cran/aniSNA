net_bootstrap=function(M, n_versions=1000, seed=12345) {
  N=nrow(M)
  b.M=list()
  set.seed(seed)
  for (i in 1:n_versions) {
    s=sample(1:N,replace=T)
    b.M[[i]]=M;diag(b.M[[i]])=sample(b.M[[i]][-seq(1,N^2,N+1)],N,replace=T);b.M[[i]]=b.M[[i]][s,s];diag(b.M[[i]])=diag(M)[s]
  }
  set.seed(NULL)
  return(list(orig.net=M, boot.nets=b.M))
}

compare_boot_nets=function(bootnets1, bootnets2, network_metrics_functions_list) { 
  ans=NULL
  boot1.stats=sapply(1:length(bootnets1$boot.nets), function(i) network_metrics_evaluate_from_adjacency_matrix(bootnets1$boot.nets[[i]], network_metrics_functions_list)) # generates all the bootstrap statistics
  boot2.stats=sapply(1:length(bootnets2$boot.nets), function(i) network_metrics_evaluate_from_adjacency_matrix(bootnets2$boot.nets[[i]], network_metrics_functions_list))
  
  if(length(network_metrics_functions_list) == 1){
    boot1.stats <- matrix(boot1.stats, nrow = 1)
    boot2.stats <- matrix(boot2.stats, nrow = 1)
  }
  
  for (i in 1:length(network_metrics_functions_list)) {
    t.stat=(mean(boot1.stats[i,])-mean(boot2.stats[i,]))/sqrt(stats::var(boot1.stats[i,]) + stats::var(boot2.stats[i,]))
    ans[i]=2*stats::pt(-abs(t.stat),df=2*length(boot1.stats[i,])-2) # two-sided t-test
  }
  names(ans)=names(network_metrics_functions_list) # returns a named vector of p-values
  return(ans)
}


#Function to extract 2 subsamples of equal size and non overlapping
two_samples_network <- function(network, size){
  if(size > igraph::gorder(network)/2){
    stop("Size of sample is greater than 50%, choose a smaller size\n")
    return(NULL)
  }
  sampled_nodes_1 <- sample(igraph::V(network), size = size)
  sampled_nodes_2 <- sample(igraph::V(network)[-sampled_nodes_1], size = size)
  
  sub_network_1 <- igraph::induced_subgraph(network, sampled_nodes_1, impl = "auto")
  sub_network_2 <- igraph::induced_subgraph(network, sampled_nodes_2, impl = "auto")
  return(list(sub_network_1, sub_network_2))
}

#Function to calculate network statistics with input of adjacency matrix
network_metrics_evaluate_from_adjacency_matrix=function(x, funcs) {
  x <- igraph::graph_from_adjacency_matrix(x, mode = "undirected", weighted = TRUE)
  return(unlist(lapply(funcs, function(f) f(x))))
}

#Function to pass a network and obtain p value matrix
p_value_matrix <- function(network, size_subnet, n.iter = 10, network_metrics_functions_list, n_versions){
  
  metrics_pvalue <- data.frame(temp = numeric(0))
  for(i in 1:length(network_metrics_functions_list)){
    metrics_pvalue[[names(network_metrics_functions_list)[i]]] <- numeric(0)
  }
  metrics_pvalue <- metrics_pvalue[,-1] #remove the temporary column
  
  for(t in 1:n.iter){
    set.seed(t)
    g_sub <- two_samples_network(network, size_subnet)
    #compare these two networks
    g1_matrix <- igraph::as_adjacency_matrix(g_sub[[1]], type = "both",  sparse = FALSE, attr = "weight")
    g2_matrix <- igraph::as_adjacency_matrix(g_sub[[2]], type = "both",  sparse = FALSE, attr = "weight")
    g1_boot=net_bootstrap(g1_matrix, n_versions = n_versions) 
    g2_boot=net_bootstrap(g2_matrix, n_versions = n_versions)
    metrics_pvalue[t,] <- as.vector(compare_boot_nets(g1_boot, g2_boot, network_metrics_functions_list))
  }
  return(metrics_pvalue)
}

#Function to obtain a subnetwork of the network passed as first argument of size passed as second argument
sub_sample_network <- function(network, size){
  if(size > igraph::gorder(network)){
    stop("Size of subnetwork is greater than actual network, choose a smaller size\n")
    return(NULL)
  }
  sampled_nodes <- sample(igraph::V(network), size = size)
  sub_network <- igraph::induced_subgraph(network, sampled_nodes, impl = "auto")
  return(sub_network)
}

#A function that takes booststrapped samples and a function to obtain the required network statistics
#and obtains a vector of width of confidence intervals for the estimates
CI_width_nets=function(bootnets, network_metrics_functions_list, CI_size) { 
  ans=NULL
  boot.stats=sapply(1:length(bootnets$boot.nets), function(i) network_metrics_evaluate_from_adjacency_matrix(bootnets$boot.nets[[i]], network_metrics_functions_list)) # generates all the bootstrap statistics
  
  if(length(network_metrics_functions_list) == 1) boot.stats <- matrix(boot.stats, nrow = 1)
  
  for (i in 1:length(boot.stats[,1])) {
    quant <- stats::quantile(boot.stats[i,], probs=c((1-CI_size)/2, 0.5 + CI_size/2), na.rm = TRUE)
    ans[i] <- quant[2] - quant[1] #width of CI
  }
  names(ans)=names(boot.stats[,1])
  return(ans)
}

CI_matrix <- function(network, size_subnet, n_versions, n.iter = 10, network_metrics_functions_list, CI_size){
  
  width_CI <- data.frame(temp = numeric(0))
  for(i in 1:length(network_metrics_functions_list)){
    width_CI[[names(network_metrics_functions_list)[i]]] <- numeric(0)
  }
  width_CI <- width_CI[,-1]
  
  for(t in 1:n.iter){
    set.seed(t)
    #1. Generates subsample of given size n.iter number of times
    g_sub <- sub_sample_network(network, size_subnet)
    
    g_matrix <- igraph::as_adjacency_matrix(g_sub, type = "both", sparse = FALSE, attr = "weight")
    
    #2. Generates booststrapped versions from each of those subnetworks
    g_boot = net_bootstrap(g_matrix, n_versions = n_versions)
    
    #3. Computes CI width for each subnetwork's bootstrapped versions
    width_CI[t,] <- CI_width_nets(g_boot, network_metrics_functions_list, CI_size)
  }
  #4. Takes in all the CI lengths and returns a dataframe with rows equal to number of iterations 
  #and columns equal to number of network statistics
  return(width_CI)
}

