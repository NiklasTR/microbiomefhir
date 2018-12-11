#' Parse FHIR report from phyloseq object
#'
#' @param template_path
#'
#' @return
#' @export
#' @import here
#' @import jsonlite
#' @import tidyverse
#'
#' @examples
parse_json <- function(phyloseq_obj, file, column = 1, template_path = here("data/fhir_template.RData")){
  #loading template
  load(template_path)

  df <- otu_table(phyloseq_obj)[, column] %>% as.data.frame() %>%
    magrittr::set_colnames("abundance") %>%
    rownames_to_column("id")

  template <- obj_syn %>%
    trim_template()

  #creating results
  res_df <- create_res(df$id)

  # creating obs
  obs_df <- base::mapply(create_obs, percent_abundance = df$abundance, species_id = df$id, SIMPLIFY = FALSE)
  # appending obs
  obs_df <- append(obs_df, template$contained, 0)


  # setting ids - most of the time stool, bacterium 1, 2, 3, ...
  #template$contained$id <- c(template$contained$id, df$id)
  template$contained <- obs_df
  # setting obs
  template$result <- res_df

  # writing JSON
  template %>% toJSON(dataframe = "rows", pretty = TRUE) %>% write(file)

}
