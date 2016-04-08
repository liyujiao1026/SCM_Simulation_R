# SCM_plot.R script 
#
# the purpose of this script is to generate the plot in shiny app
#
# first version: 160404
# this version:  160404
# last change by: Yujiao Li

#==================================================================================================#

library(tidyr)
library(ggplot2)

# raw.data is generated from DataIntegrate.plot()

Plot.SCM <- function(raw.data, invTime, randomOption = FALSE){   
            
            data.plot <- gather(raw.data, Unit, Value,-Time)
            
            plot <- ggplot(data.plot, aes(x = Time, y = Value, colour = Unit, group = Unit)) +
                        geom_line(lwd = 1, alpha = 0.3,linetype = 1) +
                        
                        
                        
                        geom_line(
                                    data = data.plot[data.plot$Unit == "Treat.Obs",],
                                    aes(x = Time, y = Value), colour = "black",
                                    lwd = 1.5, linetype = "solid" ) +
                        
                        geom_line(
                                    data = data.plot[data.plot$Unit == "Treat.Hat",],
                                    aes(x = Time, y = Value), colour = "royalblue3",
                                    lwd = 1.5,linetype = "longdash" ) +
                        
                        geom_line(
                                    data = data.plot[data.plot$Unit == "Treat.Obs.Rand",],
                                    aes(x = Time, y = Value), colour = "sienna1",
                                    lwd = 2.2,linetype = "dashed" ) +
         
                        # add points             
                        geom_point(
                                    data = data.plot[data.plot$Unit == "Treat.Hat",],
                                    aes(x = Time, y = Value), colour = "blue",size = 2 ) +
                        geom_point(
                                    data = data.plot[data.plot$Unit == "Treat.Obs",],
                                    aes(x = Time, y = Value), colour = "black",size = 2 ) +
                        geom_point(
                                    data = data.plot[data.plot$Unit == "Treat.Obs.random",],
                                    aes(x = Time, y = Value), colour = "darkorange",size = 2 ) +
                        
                        # add line and text
                        geom_vline(xintercept = invTime,linetype = "longdash") +
                        
                        ggtitle("Synthetic Control Method\nSimulation") +
                        
                        annotate(
                                    "text", label = paste0("Intervention \n  =  ",invTime),
                                    x = invTime, y = (0.05 + min(data.plot$Value)),
                                    size = 4, colour = "red") +
                        theme(legend.position = "none")
            
            
                        if (randomOption == TRUE) {
                                    plot <- plot + 
                                                geom_line(
                                                data = data.plot[data.plot$Unit == "Treat.Hat_SCM",],
                                                aes(x = Time, y = Value), colour = "darkred",
                                                lwd = 1.8, linetype = "dashed") +
                                    
                                                geom_line(
                                                data = data.plot[data.plot$Unit == "Treat.Hat_SCM.Rand",],
                                                aes(x = Time, y = Value), colour = "deeppink",
                                                lwd = 1.5, linetype = "solid")
                        
                        }
            
            return(plot)
            
}


