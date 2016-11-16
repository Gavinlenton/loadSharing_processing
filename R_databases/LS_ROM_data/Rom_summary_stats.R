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


# Initial data checks -----------------------------------------------------


# Check all variables for outliers
summary(hip_flexion_ROM)
summary(trunk_flexion_ROM)
summary(shoulder_FF_ROM)
summary(UUA_ROM)

# Convert super low values to NA
hip_flexion_ROM[,2:14][hip_flexion_ROM[,2:14] < 55] <- NA
shoulder_FF_ROM[,2:14][shoulder_FF_ROM[,2:14] < 90] <- NA
UUA_ROM[,2:14][UUA_ROM[,2:14] < 90] <- NA
trunk_flexion_ROM[,2:14][trunk_flexion_ROM[,2:14] < 70] <- NA

# Omit NAs listwise
hip_flexion_ROM <- na.omit(hip_flexion_ROM)
trunk_flexion_ROM <- na.omit(trunk_flexion_ROM)
shoulder_FF_ROM <- na.omit(shoulder_FF_ROM)
UUA_ROM <- na.omit(UUA_ROM)

# Load libraries for ANOVA
library(mosaic)
library(pastecs)
library(ez)
library(reshape2)

# Repeated measure ANOVA --------------------------------------------------

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
hip_flexion_ROM.l$Mass = gl(2, 102, labels = c("15", "30"))
hip_flexion_ROM.l$Armour_type = gl(6, 17, 204, labels = c("TBAS", "cARM1", "cARM2", "pARM1", "pARM2", "pARM3"))

UUA_ROM.l$Mass = gl(2, 120, labels = c("15", "30"))
UUA_ROM.l$Armour_type = gl(6, 20, 240, labels = c("TBAS", "cARM1", "cARM2", "pARM1", "pARM2", "pARM3"))

shoulder_FF_ROM.l$Mass = gl(2, 108, labels = c("15", "30"))
shoulder_FF_ROM.l$Armour_type = gl(6, 18, 216, labels = c("TBAS", "cARM1", "cARM2", "pARM1", "pARM2", "pARM3"))

trunk_flexion_ROM.l$Mass = gl(2, 114, labels = c("15", "30"))
trunk_flexion_ROM.l$Armour_type = gl(6, 19, 228, labels = c("TBAS", "cARM1", "cARM2", "pARM1", "pARM2", "pARM3"))

# Run the ANOVA
# Hip flexion
model_hip_flex = ezANOVA(data = hip_flexion_ROM.l, 
                         dv = hip_flex_ROM, 
                         wid = Participants, 
                         within = .(Armour_type, Mass), 
                         detailed = TRUE, 
                         type = 3, 
                         return_aov = TRUE)

# UUA
model_UUA = ezANOVA(data = UUA_ROM.l, 
                         dv = UUA_ROM, 
                         wid = Participants, 
                         within = .(Armour_type, Mass), 
                         detailed = TRUE, 
                         type = 3, 
                         return_aov = TRUE)

# Shoulder forward flexion
model_shoulder_FF= ezANOVA(data = shoulder_FF_ROM.l, 
                    dv = shoulder_FF_ROM, 
                    wid = Participants, 
                    within = .(Armour_type, Mass), 
                    detailed = TRUE, 
                    type = 3, 
                    return_aov = TRUE)

# Trunk flexion
model_trunk_flex = ezANOVA(data = trunk_flexion_ROM.l, 
                    dv = trunk_flex_ROM, 
                    wid = Participants, 
                    within = .(Armour_type, Mass), 
                    detailed = TRUE, 
                    type = 3, 
                    return_aov = TRUE)

# Post hoc tests


# Post hoc hip flexion ----------------------------------------------------


# Bonferroni for armour type
pairwise.t.test(hip_flexion_ROM.l$hip_flex_ROM,
                hip_flexion_ROM.l$Armour_type,
                paired = TRUE,
                p.adjust.method = "bonferroni")

# Bonferroni for mass
pairwise.t.test(hip_flexion_ROM.l$hip_flex_ROM,
                hip_flexion_ROM.l$Mass,
                paired = TRUE,
                p.adjust.method = "bonferroni")

# Interaction simple effects analysis
# Divide the data into mass levels
hip_flexion_15 = subset(hip_flexion_ROM.l, Mass == "15")
hip_flexion_30 = subset(hip_flexion_ROM.l, Mass == "30")

# Bonferroni for armour type
pairwise.t.test(hip_flexion_15$hip_flex_ROM,
                hip_flexion_15$Armour_type,
                paired = TRUE,
                p.adjust.method = "bonferroni")

pairwise.t.test(hip_flexion_30$hip_flex_ROM,
                hip_flexion_30$Armour_type,
                paired = TRUE,
                p.adjust.method = "bonferroni")


# Post hoc trunk flexion --------------------------------------------------


# Bonferroni for armour type
pairwise.t.test(trunk_flexion_ROM.l$trunk_flex_ROM,
                trunk_flexion_ROM.l$Armour_type,
                paired = TRUE,
                p.adjust.method = "bonferroni")

# Bonferroni for mass
pairwise.t.test(trunk_flexion_ROM.l$trunk_flex_ROM,
                trunk_flexion_ROM.l$Mass,
                paired = TRUE,
                p.adjust.method = "bonferroni")

# Interaction simple effects analysis
# Divide the data into mass levels
trunk_flexion_15 = subset(trunk_flexion_ROM.l, Mass == "15")
trunk_flexion_30 = subset(trunk_flexion_ROM.l, Mass == "30")

# Bonferroni for armour type
pairwise.t.test(trunk_flexion_15$trunk_flex_ROM,
                trunk_flexion_15$Armour_type,
                paired = TRUE,
                p.adjust.method = "bonferroni")

pairwise.t.test(trunk_flexion_30$trunk_flex_ROM,
                trunk_flexion_30$Armour_type,
                paired = TRUE,
                p.adjust.method = "bonferroni")


# Post hoc shoulder FF ----------------------------------------------------

## ANOVA only showed significant value for mass
# Bonferroni for mass
pairwise.t.test(shoulder_FF_ROM.l$shoulder_FF_ROM,
                shoulder_FF_ROM.l$Mass,
                paired = TRUE,
                p.adjust.method = "bonferroni")


# Post hoc shoulder abduction ---------------------------------------------
## ANOVA only showed significant value for armour type
# Bonferroni for mass
pairwise.t.test(UUA_ROM.l$UUA_ROM,
                UUA_ROM.l$Armour_type,
                paired = TRUE,
                p.adjust.method = "bonferroni")

UUA_pARM2 = subset(UUA_ROM.l, Armour_type == "pARM2")
UUA_TBAS = subset(UUA_ROM.l, Armour_type == "TBAS")
UUA_pARM3 = subset(UUA_ROM.l, Armour_type == "pARM3")

mean_pARM2 = mean(UUA_pARM2$UUA_ROM); sd_pARM2 = sd(UUA_pARM2$UUA_ROM)
mean_TBAS = mean(UUA_TBAS$UUA_ROM); sd_TBAS = sd(UUA_TBAS$UUA_ROM)
mean_pARM3 = mean(UUA_pARM3$UUA_ROM); sd_pARM3 = sd(UUA_pARM3$UUA_ROM)

# Summary_SE function -----------------------------------------------------


# Generate descriptive statistics
library(ggplot2)

## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
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

# Plot means ------------------------------------------------------

summary_hip_flex <- summarySE(data = hip_flexion_ROM.l, measurevar="hip_flex_ROM", groupvars=c("Mass","Armour_type"))
summary_trunk_flex <- summarySE(data = trunk_flexion_ROM.l, measurevar="trunk_flex_ROM", groupvars=c("Mass","Armour_type"))
summary_shoulder_FF <- summarySE(data = shoulder_FF_ROM.l, measurevar="shoulder_FF_ROM", groupvars=c("Mass","Armour_type"))
summary_UUA <- summarySE(data = UUA_ROM.l, measurevar="UUA_ROM", groupvars=c("Mass","Armour_type"))

# Use armour_type as a factor rather than numeric
summary_hip_flex2 <- summary_hip_flex
summary_hip_flex2$Armour_type <- factor(summary_hip_flex2$Armour_type)

summary_trunk_flex2 <- summary_trunk_flex
summary_trunk_flex2$Armour_type <- factor(summary_trunk_flex2$Armour_type)

summary_shoulder_FF2 <- summary_shoulder_FF
summary_shoulder_FF2$Armour_type <- factor(summary_shoulder_FF2$Armour_type)

summary_UUA2 <- summary_UUA
summary_UUA2$Armour_type <- factor(summary_UUA2$Armour_type)

# Plot with error bars represent standard error of the mean
# Hip flexion
plot_HF <- ggplot(summary_hip_flex2, aes(x=Armour_type, y=hip_flex_ROM, fill=Mass)) + 
  geom_bar(position=position_dodge(0.9), colour = "black", stat="identity") +
  geom_errorbar(aes(ymin=hip_flex_ROM-se, ymax=hip_flex_ROM+se),
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
plot_SFF <- ggplot(summary_shoulder_FF2, aes(x=Armour_type, y=shoulder_FF_ROM, fill=Mass)) + 
  geom_bar(position=position_dodge(0.9), colour = "black", stat="identity") +
  geom_errorbar(aes(ymin=shoulder_FF_ROM-se, ymax=shoulder_FF_ROM+se),
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
plot_TF <- ggplot(summary_trunk_flex2, aes(x=Armour_type, y=trunk_flex_ROM, fill=Mass)) + 
  geom_bar(position=position_dodge(0.9), colour = "black", stat="identity") +
  geom_errorbar(aes(ymin=trunk_flex_ROM-se, ymax=trunk_flex_ROM+se),
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
plot_UUA <- ggplot(summary_UUA2, aes(x=Armour_type, y=UUA_ROM, fill=Mass)) + 
  geom_bar(position=position_dodge(0.9), colour = "black", stat="identity") +
  geom_errorbar(aes(ymin=UUA_ROM-se, ymax=UUA_ROM+se),
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
png(file="ROM_summary.png",width = 8, height = 5, units = 'in', res = 300)
multiplot(plot_HF, plot_SFF, plot_TF, plot_UUA, cols=2)

dev.off()

