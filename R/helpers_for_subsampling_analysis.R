#Function to evaluate network metrics from the function definitions passed by the user
network_metrics_evaluate=function(x, funcs) {return(lapply(funcs, function(f) f(x)))}