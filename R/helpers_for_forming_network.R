
# Function to calculate association index based on modified simple ratio index for GPS telemetry data
get_AI <- function(i, number_interactions, species_raw) {

  # get the ids for two individuals in the selected pair
  animal_A <- number_interactions$Animal_A[i]
  animal_B <- number_interactions$Animal_B[i]

  # get all the dates when individuals A and B are observed
  dates_A <- lubridate::date(species_raw$datetime[species_raw$animal_id == animal_A])
  dates_B <- lubridate::date(species_raw$datetime[species_raw$animal_id == animal_B])

  # Count number of sightings of A(B) only when B(A) is also observed
  total_observations <- (sum(dates_A %in% dates_B) + sum(dates_B %in% dates_A) - number_interactions$n[i]) # removing number_interactions between A and B
  # as it has been accounted for twice while counting

  return(number_interactions$n[i] / total_observations)
}

# Function to obtain network structure from raw GPS observations and interactions captured.
network_obtain <- function(species_raw, weighted_interactions) {
  species_network <- igraph::graph_from_data_frame(weighted_interactions,
    directed = FALSE,
    vertices = unique(species_raw$animal_id)
  )
  return(species_network)
}