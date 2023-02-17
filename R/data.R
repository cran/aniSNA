### Dataset documentation : elk_data_2010 ###
#'
#' Data to showcase functions in our package
#'
#' Contains GPS telemetry observations of the species elk in year 2010
#'
#' @format A dataframe with 123568 rows and 4 variables:
#' \describe{
#' \item{animal_id}{Unique ID of individuals in the observed sample}
#' \item{datetime}{Date and timestamp of the observation}
#' \item{latitude_rad}{Latitude of individual observation in radians}
#' \item{longitude_rad}{Longitude of individual observation in radians}
#' }
#'
#' @examples
#' data(elk_data_2010)
"elk_data_2010"


#' Dataset of interactions from elk_data_2010 using first mode as the spatial threshold
#' 
#' @format A dataframe with 2393 rows and 5 variables
#' \describe{
#' \item{Animal_A}{First animal ID}
#' \item{Animal_B}{Second animal ID}
#' \item{Timestamp_A}{Observation timestamp of first animal}
#' \item{Timestamp_B}{Observation timestamp of second animal}
#' \item{distance}{Distance in metres between the two animals}
#' }
#' 
#' @examples 
#' data(elk_interactions_2010)
"elk_interactions_2010"


#' Dataset of all possible interactions from elk_data_2010
#' 
#' @format A dataframe with 7615 rows and 5 variables
#' \describe{
#' \item{Animal_A}{First animal ID}
#' \item{Animal_B}{Second animal ID}
#' \item{Timestamp_A}{Observation timestamp of first animal}
#' \item{Timestamp_B}{Observation timestamp of second animal}
#' \item{distance}{Distance in metres between the two animals}
#' }
#' 
#' @examples 
#' data(elk_all_interactions_2010)
"elk_all_interactions_2010"


#' An igraph object depicting the network obtained from elk_interactions_2010
#' 
#' @format An igraph object with 57 nodes and 114 edges
#' 
#' @examples
#' igraph::E(elk_network_2010)
"elk_network_2010"


#' A list of 100 igraph objects obtained by permuting the raw elk_data_2010 and obtaining network from those
#' 
#' @format A list of 100 igraph objects
#' 
#' @examples
#' data(elk_2010_permutations)
"elk_2010_permutations"