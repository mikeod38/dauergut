#' format_dauer
#'
#' formats dauer data for analysis requires variables 'strains' and 'foods'
#' @param df csv or data frame (can be piped)
#' @param p.dauer include post an/or partial dauers. Arguments "exclude" (default) - not counted,
#' "non" - animals are considered non-dauer in cases where 100% dauers are observed,
#' "dauer" - animals are counted as dauers in cases where only partial dauers form (daf-16)
#' @importFrom magrittr "%>%"
#' @importFrom magrittr "%<>%"
#' @export
#' @examples df %>% format_dauer(p.dauer = "exclude")
#' 
#' 
format_dauer <- function(df, p.dauer) {
  
  df %<>% dplyr::filter(genotype %in% strains, food %in% foods) %>%
    dplyr::mutate(genotype = factor(genotype, levels=strains),
       non.dauer = as.numeric(paste(pd + non)),
       strainDate = interaction(genotype, day),
       plateID = interaction(strainDate, plate),
       food = factor(food, levels = foods))
  
if(p.dauer == "exclude" | missing(p.dauer)) {
  df %<>% dplyr::mutate(
    n = as.numeric(paste(dauer + non)),
    pct = as.numeric(paste(dauer/n)))
} else { 
  if(p.dauer == "non") {
    df %<>% dplyr::mutate(
      n = as.numeric(paste(dauer + non.dauer)),
      pct = as.numeric(paste(dauer/n)))
  } else {
    if(p.dauer == "dauer") {
      df %<>% dplyr::mutate(
        n = as.numeric(paste(dauer + non.dauer)),
        pct = as.numeric(paste((dauer+pd)/n)))
  }
}
}
}


