

#' To convert latitude and longitude values from degrees to radians
#'
#' @param species_raw A DataFrame consisting of GPS observations. The DataFrame must have a "latitude" column and a "longitude" column
#' whose values are specified in degrees. 
#'
#' @return The same DataFrame that has been passed as the argument with two additional columns namely "latitude_rad" and "longitude_rad"
#' @export
#'
get_coordinates_in_radian <- function(species_raw) {
  species_raw$latitude_rad <-  0.01745329251 * species_raw$latitude
  species_raw$longitude_rad <-  0.01745329251 * species_raw$longitude
  return(species_raw)
}


#' To obtain interactions from raw GPS observations
#'
#' @param species_raw A DataFrame consisting of GPS observations. 
#' It should have at least four columns namely "animal_id", "datetime", "latitude|_rad", and "longitude_rad". 
#' "latitude|_rad", and "longitude_rad" are latitude and longitude values in radians respectively. See function "get_coordinates_in_radian"
#' to get these values.
#' @param temporal_thresh Temporal threshold in minutes with default 7 minutes
#' @param spatial_thresh The maximum distance in meters within which two animals are considered interacting
#' @param n_cores Number of cores for parallel processing with default 1
#'
#' @return A dataframe consisting of five columns. The first two columns contain animal ids, third and fourth column contain timestamp of their observations and the final column contains the distance between the two individuals.
#' @export
#'
#' @examples
#' \donttest{
#' data(elk_data_2010)
#' get_interactions(elk_data_2010, temporal_thresh = 7, spatial_thresh = 15)
#' }
get_interactions <- function(species_raw, temporal_thresh = 7, spatial_thresh, n_cores = 1) {

  # Sorting the observations by datetime values
  species_raw <- species_raw[with(species_raw, order(species_raw$datetime)), ][1:nrow(species_raw), ]

  # For each row i in data, we get corresponding subsequent rows that are within the user provided temporal and spatial threshold
  list_ijs <- parallel::mclapply(1:(nrow(species_raw) - 1), function(i) {
      interacting_pairs(
      i - 1, species_raw$datetime,
      species_raw$latitude_rad,
      species_raw$longitude_rad, temporal_thresh, spatial_thresh
    )
  },
  mc.cores = n_cores
  )

  # Drop the null lists
  keep <- !sapply(list_ijs, is.null)
  keep <- c(keep, FALSE)
  interactions_new <- data.frame(
    Animal_A = rep(species_raw$animal_id[keep], unlist(sapply(list_ijs[keep], function(x) ncol(x)))),
    Animal_B = species_raw$animal_id[unlist(sapply(list_ijs[keep], function(x) x[1, ]))],
    Timestamp_A = rep(species_raw$datetime[keep], unlist(sapply(list_ijs[keep], function(x) ncol(x)))),
    Timestamp_B = species_raw$datetime[unlist(sapply(list_ijs[keep], function(x) x[1, ]))],
    distance = unlist(sapply(list_ijs[keep], function(x) x[2, ])), stringsAsFactors = FALSE
  )
  interactions_new <- interactions_new[which(interactions_new$Animal_A != interactions_new$Animal_B), ]

  return(interactions_new)
}