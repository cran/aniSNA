#The following function permutes the order of dates
obtain_permutation_dates=function(species_raw, n_permutations = 100, n_cores = 1) {
  
  #order the species by animal ID
  species_raw <- species_raw[order(species_raw$animal_id),]
  
  #Extract date, hour and minute
  species_raw$date = as.character(lubridate::date(species_raw$datetime)) #it is faster to compare strings
  species_raw$hour = lubridate::hour(species_raw$datetime)
  species_raw$minute = lubridate::minute(species_raw$datetime)
  
  permuted_dates_list <- vector(mode = "list", length = n_permutations)
  
  # first create a subset data frame containing only the unique pairs of ids and dates
  u=unique(subset(species_raw, select=c("animal_id", "date")))
  
  permuted_dates_list <- parallel::mclapply(permuted_dates_list, 
                                  function(permuted_dates_list) 
                                  {u$permuted_date=unlist(sapply(unique(u$animal_id), function(animal_id) sample(u$date[u$animal_id==animal_id])))
                                  permuted_dates = u$permuted_date[match(species_raw$date, u$date)]
                                  #Make a permuted datetime column
                                  permuted_dates_list$permuted_datetime <- as.vector(paste(permuted_dates, " ", 
                                                                                           formatC(species_raw$hour, width = 2, flag = "0"), ":",
                                                                                           formatC(species_raw$minute, width = 2, flag = "0"), sep = ""))
                                  return(permuted_dates_list$permuted_datetime)
                                  }, mc.cores=n_cores)
  return(permuted_dates_list)
}
