#' Create Observations in Microbiome FHIR report
#'
#' @param patient_name
#' @param ymd_hms_issued
#' @param snomed_id
#' @param percent_abundance
#' @param species_id
#'
#' @return
#' @export
#' @import tidyverse
#' @import jsonlite
#'
#' @examples
create_obs <- function(patient_name = "Roel",
                       ymd_hms_issued = "2019-08-16T07:03:00Z",
                       snomed_id = NA,
                       percent_abundance,
                       species_id){

  species_name <- species_id %>% str_split(pattern = "__") %>% unlist() %>% tail(1)

  data.frame(
    resourceType = "Observation",
    id = species_id,
    status = "final",
    #type = ,
    subject = I(data.frame(reference = NA,
                           display = patient_name)),
    issued = lubridate::ymd_hms(ymd_hms_issued),
    #collection = ,
    valueCodeableConcept = I(data.frame(coding = I(list(data.frame(system = "http://snomed.info/sct",
                                                                   code = snomed_id,
                                                                   display = species_name))))),
    valueQuantity = I(data.frame(value = percent_abundance,
                                 unit = "percent",
                                 code = "%"))
  ) %>% return()
}
