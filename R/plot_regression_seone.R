#' plot_regression_seone
#'
#' plot dauer formation (x) vs. absolute latitidue (y)
#' @param df input dataset. Requires the following columns: meanX (pheromone), meanY, semX, and abs_lat
#' @export 
#' @examples plot_regression_seone(df, "mean.HT","SEM.HT")
plot_regression_seone<-function(df, meanX, meanY, semX) {
  xSEmin = df[[meanX]]-df[[semX]]
  xSEmax = df[[meanX]]+df[[semX]]
  p <- ggplot(df, aes_string(x=meanX, y=meanY),environment = environment())
  p + geom_point() + 
    geom_errorbarh(aes(xmin=xSEmin, xmax=xSEmax)) +
    coord_cartesian(xlim=c(-.05,0.9)) + 
    theme_my +
    geom_smooth(aes_string(x=meanX, y=meanY, group=1), method="lm", se=TRUE, alpha = 0.2) +
    ggrepel::geom_text_repel(aes(label=genotype),
                    box.padding = unit(0.35, "lines"),
                    point.padding = unit(0.5, "lines"),
                    segment.color = 'grey50'
    ) +
    labs(
      title = 'Dauer formation and isolation latitude are uncorrelated',
      x = "High temp dauer formation
      (proportion dauers)",
      y = "Absolute latitude
      (deg)"
    )
}