
#' Function to obtain a network structure from interactions dataframe
#'
#' @param species_raw A dataframe consisting of raw GPS observations.
#' It should have at least four columns namely "animal_id", "datetime", "latitude|_rad", and "longitude_rad". 
#' "latitude|_rad", and "longitude_rad" are latitude and longitude values in radians respectively. See function "get_coordinates_in_radian"
#' to get these values.
#' @param interactions A dataframe of interactions obtained from raw GPS observations using the function "get_interactions"
#' @param n_cores Number of cores for parallel processing, default is 1
#'
#' @return An object of class igraph 
#' @export
#'
#' @examples
#' data(elk_data_2010, elk_interactions_2010)
#' network_from_interactions(elk_data_2010, elk_interactions_2010)
#' 
#' @importFrom rlang .data
network_from_interactions <- function(species_raw, interactions, n_cores = 1) {
  for (i in which(interactions$Animal_A > interactions$Animal_B)) {
    swap_variable_id <- interactions$Animal_A[i]
    interactions$Animal_A[i] <- interactions$Animal_B[i]
    interactions$Animal_B[i] <- swap_variable_id

    swap_variable_time <- interactions$Timestamp_A[i]
    interactions$Timestamp_A[i] <- interactions$Timestamp_B[i]
    interactions$Timestamp_B[i] <- swap_variable_time
  }
  
  by_animal_AandB <- interactions %>% dplyr::group_by(.data$Animal_A, .data$Animal_B)
  number_interactions <- by_animal_AandB %>% dplyr::summarise(n = dplyr::n())
  number_interactions$weight <- unlist(parallel::mclapply(1:nrow(number_interactions),
    function(i) {
      get_AI(i, number_interactions, species_raw)
    },
    mc.cores = n_cores
  ))
  species_network <- network_obtain(species_raw, number_interactions)
  
  return(species_network)
}

#' @importFrom magrittr %>%
NULL
