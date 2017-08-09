# Basic data analysis for ROM data
 # Clear workspace
rm(list = ls())
# ctrl + l to clear the screen

# set working directory
setwd("~/Google Drive/LS_main_data_collection/ElaboratedData/ROM/")

# Load files
hip_flexion_ROM <- read.csv("hip_flexion_ROM.csv", header = TRUE, sep = ",", na.strings = "NaN")
shoulder_FF_ROM <- read.csv("shoulder_FF_ROM.csv", header = TRUE, sep = ",", na.strings = "NaN")
UUA_ROM <- read.csv("UUA_ROM.csv", header = TRUE, sep = ",", na.strings = "NaN")
trunk_flexion_ROM <- read.csv("trunk_flexion_ROM.csv", header = TRUE, sep = ",", na.strings = "NaN")

# Load libraries for multiple imputation by chained equations
library(Rcpp)
library(mice)
library(VIM)

# Multiple plots function -------------------------------------------------
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

missingValueImputation <- function(data2impute){
  
  # Check the pattern of missing values
  patternMissing = md.pattern(data2impute)
  
  # Plot missing values
  #missing_plot <- aggr(data2impute, col=c('navyblue','yellow'),
                       #numbers=TRUE, sortVars=TRUE,
                       #labels=names(data2impute), cex.axis=.7,
                       #gap=3, ylab=c("Missing data","Pattern"))
  
  # maxit – Refers to no. of iterations taken to impute missing values
  # method – Refers to method used in imputation. we used predictive mean matching.
  # m  – Refers to 5 imputed data sets
  imputed_Data <- mice(data2impute, m=2, maxit = 100, method = 'pmm', seed = 500);
  
  #check imputed values
  #imputed_Data$imp$peak_moment
  
  completeData <- complete(imputed_Data,1)  
  
  return(completeData);
}

# Function to detect outliers and remove - uses Tukey method (1.5 * IQR)
outlierKD <- function(dt, var) {
  var_name <- eval(substitute(var),eval(dt))
  na1 <- sum(is.na(var_name))
  m1 <- mean(var_name, na.rm = T)
  par(mfrow=c(2, 2), oma=c(0,0,3,0))
  boxplot(var_name, main="With outliers")
  hist(var_name, main="With outliers", xlab=NA, ylab=NA)
  outlier <- boxplot.stats(var_name)$out
  mo <- mean(outlier)
  var_name <- ifelse(var_name %in% outlier, NA, var_name)
  boxplot(var_name, main="Without outliers")
  hist(var_name, main="Without outliers", xlab=NA, ylab=NA)
  title(paste(i), outer=TRUE)
  na2 <- sum(is.na(var_name))
  cat("Outliers identified:", na2 - na1, "n")
  cat("Proportion (%) of outliers:", round((na2 - na1) / sum(!is.na(var_name))*100, 1), "n")
  cat("Mean of the outliers:", round(mo, 2), "n")
  m2 <- mean(var_name, na.rm = T)
  cat("Mean without removing outliers:", round(m1, 2), "n")
  cat("Mean if we remove outliers:", round(m2, 2), "n")
  response <- readline(prompt="Do you want to remove outliers and to replace with NA? [yes/no]: ")
  if(response == "y" | response == "yes"){
    dt[as.character(substitute(var))] <- invisible(var_name)
    name <- paste(i,"noOut", sep = "_")
    names(dt)[ncol(dt)] <- name
    #assign(as.character(as.list(match.call())$dt), dt, envir = .GlobalEnv)
    cat("Outliers successfully removed")
    return(dt[[name]])
  } else{
    cat("Nothing changed", "n")
    return(invisible(var_name))
  }
}

# Reshape the data for repeated measures ANOVA
hip_flexion_ROM.l = melt(hip_flexion_ROM, 
                         "Participants", 
                         measure.vars = c("TBAS15", "CRYE15", "TYR15", "USMC15", "CORE15", "SORD15", 
                                          "TBAS30", "CRYE30", "TYR30", "USMC30", "CORE30", "SORD30"), 
                         variable.name = "Armour", 
                         value.name = "hip_flex_ROM")

UUA_ROM.l = melt(UUA_ROM, 
                 "Participants", 
                 measure.vars = c("TBAS15", "CRYE15", "TYR15", "USMC15", "CORE15", "SORD15",
                                  "TBAS30", "CRYE30", "TYR30", "USMC30", "CORE30", "SORD30"), 
                 variable.name = "Armour", 
                 value.name = "UUA_ROM")

shoulder_FF_ROM.l = melt(shoulder_FF_ROM, 
                         "Participants", 
                         measure.vars = c("TBAS15", "CRYE15", "TYR15", "USMC15", "CORE15", "SORD15",
                                          "TBAS30", "CRYE30", "TYR30", "USMC30", "CORE30", "SORD30"), 
                         variable.name = "Armour", 
                         value.name = "shoulder_FF_ROM")

trunk_flexion_ROM.l = melt(trunk_flexion_ROM, 
                           "Participants", 
                           measure.vars = c("TBAS15", "CRYE15", "TYR15", "USMC15", "CORE15", "SORD15",
                                            "TBAS30", "CRYE30", "TYR30", "USMC30", "CORE30", "SORD30"), 
                           variable.name = "Armour", 
                           value.name = "trunk_flex_ROM")

# Add mass and armour type columns
hip_flexion_ROM.l$Mass = gl(2, 120, labels = c("15", "30"))
hip_flexion_ROM.l$Armour_type = gl(6, 20, 240, labels = c("TBAS", "cARM1", "cARM2", "pARM1", "pARM2", "pARM3"))

UUA_ROM.l$Mass = gl(2, 120, labels = c("15", "30"))
UUA_ROM.l$Armour_type = gl(6, 20, 240, labels = c("TBAS", "cARM1", "cARM2", "pARM1", "pARM2", "pARM3"))

shoulder_FF_ROM.l$Mass = gl(2, 120, labels = c("15", "30"))
shoulder_FF_ROM.l$Armour_type = gl(6, 20, 240, labels = c("TBAS", "cARM1", "cARM2", "pARM1", "pARM2", "pARM3"))

trunk_flexion_ROM.l$Mass = gl(2, 120, labels = c("15", "30"))
trunk_flexion_ROM.l$Armour_type = gl(6, 20, 240, labels = c("TBAS", "cARM1", "cARM2", "pARM1", "pARM2", "pARM3"))


# Initial data checks -----------------------------------------------------

# Check all variables for outliers
summary(hip_flexion_ROM)
summary(trunk_flexion_ROM)
summary(shoulder_FF_ROM)
summary(UUA_ROM)

i <- "UUA";

# Loop through vars and remove outliers
# Remove outliers
hip_flexion_ROM.l[[6]]= outlierKD(hip_flexion_ROM.l, hip_flexion_ROM.l[[3]]);
trunk_flexion_ROM.l[[6]]= outlierKD(trunk_flexion_ROM.l, trunk_flexion_ROM.l[[3]]);
UUA_ROM.l[[6]]= outlierKD(UUA_ROM.l, UUA_ROM.l[[3]]);
shoulder_FF_ROM.l[[6]]= outlierKD(shoulder_FF_ROM.l, shoulder_FF_ROM.l[[3]]);

#Name the column
names(UUA_ROM.l)[[6]] = paste(i, "noOut", sep = "_");

# Remove previous variable with missing values
hip_flexion_ROM.l$hip_flex_ROM = NULL;
trunk_flexion_ROM.l$trunk_flex_ROM = NULL;
UUA_ROM.l$UUA_ROM = NULL;
shoulder_FF_ROM.l$shoulder_FF_ROM = NULL;

# Missing value imputation
hip_flexion_ROM.final <- missingValueImputation(hip_flexion_ROM.l);
trunk_flexion_ROM.final <- missingValueImputation(trunk_flexion_ROM.l);
UUA_ROM.final <- UUA_ROM.l;
shoulder_FF_ROM.final <- missingValueImputation(shoulder_FF_ROM.l);



# Load libraries for ANOVA
library(dplyr)
library(mosaic)
library(pastecs)
library(ez)
library(reshape2)
library(effsize)

# Repeated measure ANOVA --------------------------------------------------
# Run the ANOVA
# Hip flexion
model_hip_flex = ezANOVA(data = hip_flexion_ROM.final, 
                         dv = hip_flex_noOut, 
                         wid = Participants, 
                         within = .(Armour_type, Mass), 
                         detailed = TRUE, 
                         type = 3, 
                         return_aov = TRUE)

# UUA
model_UUA = ezANOVA(data = UUA_ROM.final, 
                         dv = UUA_noOut, 
                         wid = Participants, 
                         within = .(Armour_type, Mass), 
                         detailed = TRUE, 
                         type = 3, 
                         return_aov = TRUE)

# Shoulder forward flexion
model_shoulder_FF= ezANOVA(data = shoulder_FF_ROM.final, 
                    dv = shoulder_FF_noOut, 
                    wid = Participants, 
                    within = .(Armour_type, Mass), 
                    detailed = TRUE, 
                    type = 3, 
                    return_aov = TRUE)

# Trunk flexion
model_trunk_flex = ezANOVA(data = trunk_flexion_ROM.final, 
                    dv = trunk_flex_noOut, 
                    wid = Participants, 
                    within = .(Armour_type, Mass), 
                    detailed = TRUE, 
                    type = 3, 
                    return_aov = TRUE)

# Post hoc tests



# Post hoc hip flexion ----------------------------------------------------


# Bonferroni for armour type
pairwise.t.test(hip_flexion_ROM.final$hip_flex_noOut,
                hip_flexion_ROM.final$Armour_type,
                paired = TRUE,
                p.adjust.method = "bonferroni")

# Bonferroni for mass
pairwise.t.test(hip_flexion_ROM.final$hip_flex_noOut,
                hip_flexion_ROM.final$Mass,
                paired = TRUE,
                p.adjust.method = "bonferroni")

summary_hip_flex_mass <- summarySE(data = hip_flexion_ROM.final, measurevar="hip_flex_noOut", groupvars=c("Mass"))
summary_hip_flex_armour_mass <- summarySE(data = hip_flexion_ROM.final, measurevar="hip_flex_noOut", groupvars=c("Armour_type","Mass"))

# Interaction simple effects analysis
# Divide the data into mass levels
hip_flexion_15 = subset(hip_flexion_ROM.final, Mass == "15")
hip_flexion_30 = subset(hip_flexion_ROM.final, Mass == "30")

# Bonferroni for armour type
pairwise.t.test(hip_flexion_15$hip_flex_noOut,
                hip_flexion_15$Armour_type,
                paired = TRUE,
                p.adjust.method = "bonferroni")

pairwise.t.test(hip_flexion_30$hip_flex_noOut,
                hip_flexion_30$Armour_type,
                paired = TRUE,
                p.adjust.method = "bonferroni")



# Post hoc trunk flexion --------------------------------------------------

# Bonferroni for mass
pairwise.t.test(trunk_flexion_ROM.final$trunk_flex_noOut,
                trunk_flexion_ROM.final$Mass,
                paired = TRUE,
                p.adjust.method = "bonferroni")

summary_trunk_flex_mass <- summarySE(data = trunk_flexion_ROM.final, measurevar="trunk_flex_noOut", groupvars=c("Mass"))


# Post hoc shoulder FF ----------------------------------------------------

## ANOVA only showed significant value for mass
# Bonferroni for mass
pairwise.t.test(shoulder_FF_ROM.final$shoulder_FF_noOut,
                shoulder_FF_ROM.final$Mass,
                paired = TRUE,
                p.adjust.method = "bonferroni")

summary_shoulder_FF_mass <- summarySE(data = shoulder_FF_ROM.final, measurevar="shoulder_FF_noOut", groupvars=c("Mass"))



# Post hoc shoulder abduction ---------------------------------------------

## ANOVA only showed significant value for armour type
# Bonferroni for armour type
pairwise.t.test(UUA_ROM.final$UUA_noOut,
                UUA_ROM.final$Armour_type,
                paired = TRUE,
                p.adjust.method = "bonferroni")

summary_UUA_armour<- summarySE(data = UUA_ROM.final, measurevar="UUA_noOut", groupvars=c("Armour_type"))


# Effect sizes ------------------------------------------------------------

# Hip flexion
hip_flex <- subset(joint_moments, select = c(Speed, Armour_type, Mass, hip_flex_noOut))

hip_flex_cARM1 <- subset(hip_flexion_ROM.final, Armour_type == "cARM1"); hip_flex_pARM1 <- subset(hip_flexion_ROM.final, Armour_type == "pARM1")
hip_flex_cARM2 <- subset(hip_flexion_ROM.final, Armour_type == "cARM2"); hip_flex_pARM2 <- subset(hip_flexion_ROM.final, Armour_type == "pARM2")
hip_flex_pARM3 <- subset(hip_flexion_ROM.final, Armour_type == "pARM3"); hip_flex_TBAS <- subset(hip_flexion_ROM.final, Armour_type == "TBAS")

summary_hip_flex_armour <- summarySE(data = hip_flexion_ROM.final, measurevar="hip_flex_noOut", groupvars=c("Armour_type"))

# Effect sizes
cohen.d(hip_flex_cARM1$hip_flex_noOut,hip_flex_TBAS$hip_flex_noOut, hedges.correction = TRUE)
cohen.d(hip_flex_cARM2$hip_flex_noOut,hip_flex_TBAS$hip_flex_noOut, hedges.correction = TRUE)
cohen.d(hip_flex_pARM1$hip_flex_noOut,hip_flex_TBAS$hip_flex_noOut, hedges.correction = TRUE)
cohen.d(hip_flex_pARM2$hip_flex_noOut,hip_flex_TBAS$hip_flex_noOut, hedges.correction = TRUE)
cohen.d(hip_flex_pARM3$hip_flex_noOut,hip_flex_TBAS$hip_flex_noOut, hedges.correction = TRUE)

# Trunk flexion
trunk_flex <- subset(joint_moments, select = c(Speed, Armour_type, Mass, trunk_flex_noOut))

trunk_flex_cARM1 <- subset(trunk_flexion_ROM.final, Armour_type == "cARM1"); trunk_flex_pARM1 <- subset(trunk_flexion_ROM.final, Armour_type == "pARM1")
trunk_flex_cARM2 <- subset(trunk_flexion_ROM.final, Armour_type == "cARM2"); trunk_flex_pARM2 <- subset(trunk_flexion_ROM.final, Armour_type == "pARM2")
trunk_flex_pARM3 <- subset(trunk_flexion_ROM.final, Armour_type == "pARM3"); trunk_flex_TBAS <- subset(trunk_flexion_ROM.final, Armour_type == "TBAS")

summary_trunk_flex_armour <- summarySE(data = trunk_flexion_ROM.final, measurevar="trunk_flex_noOut", groupvars=c("Armour_type"))

# Effect sizes
cohen.d(trunk_flex_cARM1$trunk_flex_noOut,trunk_flex_TBAS$trunk_flex_noOut, hedges.correction = TRUE)
cohen.d(trunk_flex_cARM2$trunk_flex_noOut,trunk_flex_TBAS$trunk_flex_noOut, hedges.correction = TRUE)
cohen.d(trunk_flex_pARM1$trunk_flex_noOut,trunk_flex_TBAS$trunk_flex_noOut, hedges.correction = TRUE)
cohen.d(trunk_flex_pARM2$trunk_flex_noOut,trunk_flex_TBAS$trunk_flex_noOut, hedges.correction = TRUE)
cohen.d(trunk_flex_pARM3$trunk_flex_noOut,trunk_flex_TBAS$trunk_flex_noOut, hedges.correction = TRUE)

# Shoulder forward flexion
shoulder_FF <- subset(joint_moments, select = c(Speed, Armour_type, Mass, shoulder_FF_noOut))

shoulder_FF_cARM1 <- subset(shoulder_FF_ROM.final, Armour_type == "cARM1"); shoulder_FF_pARM1 <- subset(shoulder_FF_ROM.final, Armour_type == "pARM1")
shoulder_FF_cARM2 <- subset(shoulder_FF_ROM.final, Armour_type == "cARM2"); shoulder_FF_pARM2 <- subset(shoulder_FF_ROM.final, Armour_type == "pARM2")
shoulder_FF_pARM3 <- subset(shoulder_FF_ROM.final, Armour_type == "pARM3"); shoulder_FF_TBAS <- subset(shoulder_FF_ROM.final, Armour_type == "TBAS")

summary_shoulder_FF_armour <- summarySE(data = shoulder_FF_ROM.final, measurevar="shoulder_FF_noOut", groupvars=c("Armour_type"))

# Effect sizes
cohen.d(shoulder_FF_cARM1$shoulder_FF_noOut,shoulder_FF_TBAS$shoulder_FF_noOut, hedges.correction = TRUE)
cohen.d(shoulder_FF_cARM2$shoulder_FF_noOut,shoulder_FF_TBAS$shoulder_FF_noOut, hedges.correction = TRUE)
cohen.d(shoulder_FF_pARM1$shoulder_FF_noOut,shoulder_FF_TBAS$shoulder_FF_noOut, hedges.correction = TRUE)
cohen.d(shoulder_FF_pARM2$shoulder_FF_noOut,shoulder_FF_TBAS$shoulder_FF_noOut, hedges.correction = TRUE)
cohen.d(shoulder_FF_pARM3$shoulder_FF_noOut,shoulder_FF_TBAS$shoulder_FF_noOut, hedges.correction = TRUE)

# Shoulder abduction
UUA <- subset(joint_moments, select = c(Speed, Armour_type, Mass, UUA_noOut))

UUA_cARM1 <- subset(UUA_ROM.final, Armour_type == "cARM1"); UUA_pARM1 <- subset(UUA_ROM.final, Armour_type == "pARM1")
UUA_cARM2 <- subset(UUA_ROM.final, Armour_type == "cARM2"); UUA_pARM2 <- subset(UUA_ROM.final, Armour_type == "pARM2")
UUA_pARM3 <- subset(UUA_ROM.final, Armour_type == "pARM3"); UUA_TBAS <- subset(UUA_ROM.final, Armour_type == "TBAS")

summary_UUA_armour <- summarySE(data = UUA_ROM.final, measurevar="UUA_noOut", groupvars=c("Armour_type"))

# Effect sizes
cohen.d(UUA_cARM1$UUA_noOut,UUA_TBAS$UUA_noOut, hedges.correction = TRUE)
cohen.d(UUA_cARM2$UUA_noOut,UUA_TBAS$UUA_noOut, hedges.correction = TRUE)
cohen.d(UUA_pARM1$UUA_noOut,UUA_TBAS$UUA_noOut, hedges.correction = TRUE)
cohen.d(UUA_pARM2$UUA_noOut,UUA_TBAS$UUA_noOut, hedges.correction = TRUE)
cohen.d(UUA_pARM3$UUA_noOut,UUA_TBAS$UUA_noOut, hedges.correction = TRUE)


# Plot means ------------------------------------------------------

summary_hip_flex <- summarySE(data = hip_flexion_ROM.final, measurevar="hip_flex_noOut", groupvars=c("Mass","Armour_type"))
summary_trunk_flex <- summarySE(data = trunk_flexion_ROM.final, measurevar="trunk_flex_noOut", groupvars=c("Mass","Armour_type"))
summary_shoulder_FF <- summarySE(data = shoulder_FF_ROM.final, measurevar="shoulder_FF_noOut", groupvars=c("Mass","Armour_type"))
summary_UUA <- summarySE(data = UUA_ROM.final, measurevar="UUA_noOut", groupvars=c("Mass","Armour_type"))

# Plot with error bars represent standard error of the mean
# Hip flexion
plot_HF <- ggplot(summary_hip_flex, aes(x=Armour_type, y=hip_flex_noOut, fill=Mass)) + 
  geom_bar(position=position_dodge(0.9), colour = "black", stat="identity") +
  geom_errorbar(aes(ymin=hip_flex_noOut-se, ymax=hip_flex_noOut+se),
                width=.2,                    # Width of the error bars
                position=position_dodge(0.9)) +
  xlab("Armour type") +
  ylab("Range of motion (deg)") +
  ggtitle("A") +
  coord_cartesian(ylim=c(60,90)) +
  scale_fill_manual(name="Mass", # Legend label, use darker colors
                    breaks=c("15", "30"),
                    labels=c("15 kg", "30 kg"),
                    values = c("#D5D5D5","#545354")) +
  scale_y_continuous(breaks=0:9*10) + 
  theme_classic(base_size = 12, base_family = "Calibri") +
  theme(legend.background = element_rect(), legend.position=c(0.7, 0.90), legend.direction = "horizontal",
        panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(), axis.ticks.x = element_blank(), axis.line.y = element_line(colour = "Black"),
        axis.line.x = element_line(colour = "Black"), axis.title.x  = element_blank(), axis.title.y = element_text(size=14), 
        plot.title = element_text(hjust = -0.01, vjust=1))


# shoulder ff plot --------------------------------------------------------


# Shoulder forward flexion
plot_SFF <- ggplot(summary_shoulder_FF, aes(x=Armour_type, y=shoulder_FF_noOut, fill=Mass)) + 
  geom_bar(position=position_dodge(0.9), colour = "black", stat="identity") +
  geom_errorbar(aes(ymin=shoulder_FF_noOut-se, ymax=shoulder_FF_noOut+se),
                width=.2,                    # Width of the error bars
                position=position_dodge(0.9)) +
  xlab("Armour type") +
  ggtitle("C") +
  ylab("Range of motion (deg)") +
  coord_cartesian(ylim=c(110,150)) +
  scale_fill_manual(name="Mass", # Legend label, use darker colors
                    breaks=c("15", "30"),
                    labels=c("15 kg", "30 kg"),
                    values = c("#D5D5D5","#545354")) +
  scale_y_continuous(breaks=0:19*10) + 
  theme_classic(base_size = 12, base_family = "Calibri") +
  theme(legend.background = element_rect(),
        panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(), axis.ticks.x = element_blank(), axis.line.y = element_line(colour = "Black"),
        axis.line.x = element_line(colour = "Black"), axis.title.y  = element_text(size=14), legend.position = "None", 
        axis.title.x  = element_blank(), plot.title = element_text(hjust = -0.01, vjust=2.12))


# trunk flexion plot ------------------------------------------------------

# Trunk flexion
plot_TF <- ggplot(summary_trunk_flex, aes(x=Armour_type, y=trunk_flex_noOut, fill=Mass)) + 
  geom_bar(position=position_dodge(0.9), colour = "black", stat="identity") +
  geom_errorbar(aes(ymin=trunk_flex_noOut-se, ymax=trunk_flex_noOut+se),
                width=.2,                    # Width of the error bars
                position=position_dodge(0.9)) +
  xlab("Armour type") +
  ggtitle("B") +
  ylab("Range of motion (deg)") +
  coord_cartesian(ylim=c(60,100)) +
  scale_fill_manual(name="Mass", # Legend label, use darker colors
                    breaks=c("15", "30"),
                    labels=c("15 kg", "30 kg"),
                    values = c("#D5D5D5","#545354")) +
  scale_y_continuous(breaks=0:10*10) + 
  theme_classic(base_size = 12, base_family = "Calibri") +
  theme(legend.background = element_rect(),
        panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(), axis.ticks.x = element_blank(), axis.line.y = element_line(colour = "Black"),
        axis.line.x = element_line(colour = "Black"), axis.title.y = element_blank(), legend.position = "None",
        axis.title.x  = element_blank(), plot.title = element_text(hjust = -0.01, vjust=2.12))


# shoulder abduction plot -------------------------------------------------

# shoulder abduction
plot_UUA <- ggplot(summary_UUA, aes(x=Armour_type, y=UUA_noOut, fill=Mass)) + 
  geom_bar(position=position_dodge(0.9), colour = "black", stat="identity") +
  geom_errorbar(aes(ymin=UUA_noOut-se, ymax=UUA_noOut+se),
                width=.2,                    # Width of the error bars
                position=position_dodge(0.9)) +
  xlab("Armour type") +
  ggtitle("D") +
  ylab("Range of motion (deg)") +
  coord_cartesian(ylim=c(120,170)) +
  scale_fill_manual(name="Mass", # Legend label, use darker colors
                    breaks=c("15", "30"),
                    labels=c("15 kg", "30 kg"),
                    values = c("#D5D5D5","#545354")) +
  scale_y_continuous(breaks=0:19*10) + 
  theme_classic(base_size = 12, base_family = "Calibri") +
  theme(legend.background = element_rect(),
        panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(), axis.ticks.x = element_blank(), axis.line.y = element_line(colour = "Black"),
        axis.line.x = element_line(colour = "Black"), axis.title.y = element_blank(), legend.position = "None", 
        axis.title.x  = element_blank(), plot.title = element_text(hjust = -0.01, vjust=2.12))

# Plot all on same figure
png(file="ROM_summary_new.png",width = 8, height = 5, units = 'in', res = 300)
multiplot(plot_HF, plot_SFF, plot_TF, plot_UUA, cols=2)

dev.off()

