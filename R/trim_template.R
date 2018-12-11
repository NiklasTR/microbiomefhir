#' Trim a FHIR template to add speciment information
#'
#' @param template
#'
#' @return
#' @export
#'
#' @examples
trim_template <- function(template){
  template$contained <- template$contained[1,]

  return(template)
}
