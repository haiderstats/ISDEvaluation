#File Name: plotSurvivalCurves.R
#Date Created: June 20, 2018
#Author: Humza Haider
#Email: hshaider@ualberta.ca

#Purpose and General Comments:
############################################################################################################################################
#Dependencies
library(ggplot2)
library(reshape2)
plotSurvivalCurves = function(survivalCurves, indexToPlot = 1, color = c(), xlim = c()){
  colorOK = T
  if(length(color) == 0)
    colorOK = F
  else if(length(color) != length(indexToPlot)){
    warning("If you would like to select custom colors please make sure the number of colors
            matches the number of curves.")
    colorOK =F 
  }
  time = survivalCurves$time
  curves = survivalCurves[,indexToPlot +1,drop=F]
  plotTimes = seq(min(time),max(time), length.out = length(time)*100)
  plotProbs = as.data.frame(sapply(curves,
                                   function(curve){
                                     curve = ifelse(curve < 1e-20,0,curve)
                                     survivialSpline = splinefun(time, curve, method = "hyman")
                                     return(pmax(survivialSpline(plotTimes),0))
                                   }
  ))
  data = cbind.data.frame(plotTimes,plotProbs)
  longFormData = melt(data,measure.vars = names(data)[-1], variable.name = "Index")
  plot = ggplot(data = longFormData, aes(x = plotTimes,y = value, colour = Index))+
    geom_line(size = 1.5)
  if(colorOK)
    plot = plot + scale_color_manual(values = color) 
  if(length(xlim)==2){
    plot = plot+ xlim(c(xlim[1],xlim[2]))
  }
  plot = plot +scale_y_continuous( limits= c(0,1),breaks = c(0,.1,.2,.3,.4,.5,.6,.7,.8,.9,1))+
    theme_bw() +
    theme(text = element_text(size=18, face=  "bold"),
          axis.title = element_text(size = 20),
          axis.title.x = element_text(margin = margin(t = 15)),
          axis.title.y = element_text(margin = margin(r = 15))) + 
    labs(y = "Survival Probaility",x = "Time" )

  return(plot)
}




