#' Load example Microbiome data
#'
#' @param name
#' @param all
#'
#' @return
#' @export
#' @import phyloseq
#' @import curatedMetagenomicData
#' @import tidyverse
#'
#' @examples
load_example <- function(name = "LomanNJ_2013.metaphlan_bugs_list.stool", all = FALSE){
  if(all == TRUE){
    print("This will be a lot of data")
  list_datasets <- curatedMetagenomicData()[grep("metaphlan_bugs_list.stool",curatedMetagenomicData())]
  data = curatedMetagenomicData(x=list_datasets,dryrun=FALSE)
  } else{

  list_datasets <- curatedMetagenomicData()[grep(name, curatedMetagenomicData())]
  data = curatedMetagenomicData(x=list_datasets,dryrun=FALSE)
  }
  data %>%
  mergeData(.) %>%
    ExpressionSet2phyloseq(.) %>%
    return()
}
