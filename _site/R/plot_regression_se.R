#' plot_regression_se
#'
#' plot correlation between pheromone dauer formation and high temp dauer formation using SEM, with a 
#' heat map for isolation latitude. 
#' @param df input dataset. Requires the following columns: meanX (pheromone), meanY (27degree),
#' semX and semY. Also a lat_abs column for latitude of each strain. 
#' @export 
#' @examples plot_regression_se(C327, "mean.C3", "mean.HT","SEM.C3", "SEM.HT")
plot_regression_se<-function(df, meanX, meanY, semX, semY) {
  ySEmin = df[[meanY]]-df[[semY]]
  ySEmax = df[[meanY]]+df[[semY]]
  xSEmin = df[[meanX]]-df[[semX]]
  xSEmax = df[[meanX]]+df[[semX]]
  p <- ggplot(df, aes_string(x=meanX, y=meanY),environment = environment())
  p + geom_point() + 
    theme_my +
    scale_color_viridis(option="plasma", end=0.9, direction = -1, name="absolute 
                        latitude") + 
    geom_errorbar(aes(ymin=ySEmin, ymax=ySEmax)) + 
    geom_errorbarh(aes(xmin=xSEmin, xmax=xSEmax)) +
    coord_cartesian(ylim = c(-.05, 1), xlim=c(-.05,1)) + 
    geom_smooth(aes_string(x=meanX, y=meanY, group=1), method="lm", se=TRUE, alpha = 0.2) +
    ggrepel::geom_text_repel(aes(label=genotype, colour=lat_abs),
                    box.padding = unit(0.35, "lines"),
                    point.padding = unit(0.5, "lines"),
                    segment.color = 'grey50'
    ) +
    labs(
      title = 'Pheromone and Temperature dependent dauer formation are
      genetically separable',
      x = "C3 dauer formation
      (proportion dauers)",
      y = "High temp dauer formation
      (proportion dauers)"
    )
}