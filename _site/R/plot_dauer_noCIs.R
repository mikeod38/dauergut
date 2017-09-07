#' plot_dauer_noCIs
#'
#' Function generates boxplots for dauer data. Optionally fills boxes based on alpha cutoff 
#' level from tukey_contrasts or dunnett_contrasts
#' @param df input dataset. Requires a "genotype" column. Response needs to be "pct" for
#' dauer plots.
#' @param title title for the plot
#' @param fillsig optional argument to include box fill colored by significance level specified in level argument
#' @param contrasts pairwise contrast object output from tukey_contrasts() or dunnett_contrasts()
#' @param level alpha level for filling box plots. White boxes are cases when p > level. 
#' @export 
#' @examples plot_dauer_noCIs(df, title = "dauer plot", plot.contrasts = plot.contrasts, ypos =  1.075)
plot_dauer_noCIs<-function (df, title, fillsig, contrasts, level) {
  if(missing(fillsig)) {
    ggplot(df, aes(x=genotype, y=pct)) +
      geom_boxplot(width=0.5) + geom_point(cex=1,alpha = 0.5) +
      labs(title = title,
           y = "proportion dauer",
           x = "genotype"
      ) + theme_my +
      scale_y_continuous(breaks=c(0,0.25,0.5,0.75, 1.0))
  } else {
    contrasts$sig<-
      if_else(
        (contrasts["p.value"] < level) & (contrasts["odds.ratio"] < 1),
        "#FDE725FF",
        if_else(
          (contrasts["p.value"] < level) & (contrasts["odds.ratio"] >1),
          "#440154FF",
          "white"))
    sig.list<-contrasts$sig
    ggplot(df, aes(x=genotype, y=pct)) +
      geom_boxplot(width=0.5,fill = c("white", sig.list)) +
      labs(title = title,
           y = "proportion dauer",
           x = "genotype"
      ) + theme_my +
      scale_y_continuous(breaks=c(0,0.25,0.5,0.75, 1.0))
  }
}

