#' To obtain spatial threshold for calculating interactions from raw GPS observations. The threshold is obtained as the 
#' distance interval that captures maximum number of inter-individual interactions. 
#'
#' @param species_interactions A dataframe consisting of individual interactions within maximum possible distance
#' @param interval_size Minimum interval size within which the number of interactions should be calculated 
#'
#' @return Spatial threshold in meters
#' @export
#'
#' @examples
#' data(elk_all_interactions_2010)
#' get_spatial_threshold(elk_all_interactions_2010, interval_size = 2)
get_spatial_threshold <- function(species_interactions, interval_size) {
  max_distance <- max(species_interactions$distance)
  breaks <- seq(0, max_distance, by = interval_size)
  distance.cut <- cut(species_interactions$distance, breaks, right = FALSE)
  distance.freq <- table(distance.cut)

  plot(breaks[-1], as.numeric(distance.freq),
    type = "b",
    col = "blue",
    main = paste("Number of observations in each interval", sep = ""),
    ylab = "Number of observations",
    xlab = "Distance"
  )

  nobs_wrt_distance <- as.data.frame(cbind(breaks[-1], as.numeric(distance.freq)))
  colnames(nobs_wrt_distance) <- c("distance", "interactions_count")
  
  first_mode <- nobs_wrt_distance$distance[which.max(nobs_wrt_distance$interactions_count)]
  return(first_mode)
}
