require(data.table)
require(ggplot2)

plot_cuminc <- function(cuminc){
  ci <- copy(cuminc)
  ci[,t0:=t0+1]
  ci <- rbind(as.list(rep(0,7)),ci)
  ci_t <- melt(ci, id.vars = c("t0"), measure.vars = list(c("Y_int", "Y_nc")),
                   value.name = c("risk"), variable.name = "treatment")
  ci_t[,treatment:=fcase(treatment=="Y_int","Intervention",
                             treatment=="Y_nc","Natural Course")]
  ci_t[,treatment:=factor(treatment,levels = c("Natural Course", "Intervention"))]
  
  ggplot(ci_t, aes(x=t0, y=risk)) + 
    geom_step(aes(color=treatment), direction = "vh") +
    scale_color_hue(direction = -1, h.start=90) +
    ylim(c(0,NA)) +
    xlab("Time") + 
    ylab("Adjusted cumulative incidence of outcome") + 
    theme_bw() + 
    labs(color="Treatment")
}
