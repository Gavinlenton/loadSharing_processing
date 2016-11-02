## Plots the input data using a column chart design
##   data: a data frame.
##   ylimits: a row vector specifying the min and max limits of the y-axis
##   xData: data to plot on x-axis
##   yData: data to plot on y-axis
##   fillData: labels to use in the legend
##   scale_cont_breaks: vector specifying the spacing of the y-axis (e.g., 0:9*10 starts at zero and increments in 10 until 90)
plot_bars_ROM = function(data = NULL, ylimits, xData = NULL, yData = NULL, fillData = NULL, 
                         scale_cont_breaks) {
  
  # Use ggplot to create the plot
  ggplot(data, aes(x=xData, y=yData, fill=fillData)) + 
    geom_bar(position=position_dodge(0.9), colour = "black", stat="identity") +
    geom_errorbar(aes(ymin=yData-se, ymax=yData+se),
                  width=.2,                    # Width of the error bars
                  position=position_dodge(0.9)) +
    xlab("Armour type") + # Labels
    ylab("Range of motion (deg)") +
    coord_cartesian(ylim=ylimits) +
    scale_fill_manual(name="Mass", # Legend label, use darker colors
                      breaks=c("15", "30"),
                      labels=c("15 kg", "30 kg"),
                      values = c("#D5D5D5","#545354")) +
    scale_y_continuous(breaks=scale_cont_breaks) + theme_light(base_size = 12, base_family = "Calibri")
}