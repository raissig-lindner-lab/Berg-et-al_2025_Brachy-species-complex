##### R Script Berg et al. 2025 - Stomatal anatomy and gas exchange dynamics in the Brachypodium distachyon complex

#### Physiology plots -----------------------------------------------------------------------------------

### load libraries
library(licornetics)
library(tidyverse)
library(ggpubr)
library(readxl)
library(ggtext)


setwd("./Li-COR files")


### 1. absolute stomatal conductance gsw
p1<- licorplots(c("bd", "bs", "bh"), timestamps = c(20,40,60), legend_title = "**Species**", 
                legend_labels = c("*B. distachyon*", "*B. stacei*", "*B. hybridum*"),
                remove_outliers = "yes", timeframe = 16:75,
                colours = met.brewer("Hokusai3", n=3),
                y_axis_limits = c(0, 0.35))

## plot normalised by stomatal density (i. e. per stoma values)
p1.2<- licorplots(c("bd", "bs", "bh"), timestamps = c(20,40,60), legend_title = "**Species**", 
                  legend_labels = c("*B. distachyon*", "*B. stacei*", "*B. hybridum*"),
                  stomden = c(100.3367, 102.693603, 59.5959596),
                  remove_outliers = "yes", timeframe = 16:75,
                  colours = met.brewer("Hokusai3", n=3),
                  y_axis_limits = c(0, 0.000006))



### 3. carbon assimilation A
p3<- licorplots(c("bd", "bs", "bh"), timestamps = c(20,40,60), legend_title = "**Species**", 
                legend_labels = c("B. distachyon", "B. stacei", "B. hybridum"),
                remove_outliers = "yes", timeframe = 16:75, 
                colours = met.brewer("Hokusai3", n=3), type = "A",
                y_axis_limits = c(-5, 24))
## plot normalised by stomatal density (i. e. per stoma values)
p3.2<- licorplots(c("bd", "bs", "bh"), timestamps = c(20,40,60), legend_title = "**Species**", 
                  legend_labels = c("*B. distachyon*", "*B. stacei*", "*B. hybridum*"),
                  stomden = c(100.3367, 102.693603, 59.5959596),
                  remove_outliers = "yes", timeframe = 16:75, 
                  colours = met.brewer("Hokusai3", n=3), type = "A",
                  y_axis_limits = c(-0.0001, 0.0004))




### 4. intrinsic water-use efficiency iWUE
p4<- licorplots(c("bd", "bs", "bh"), timestamps = c(20,40), legend_title = "**Species**", 
                legend_labels = c("*B. distachyon*", "*B. stacei*", "*B. hybridum*"),
                remove_outliers = "yes", timeframe = 16:55, 
                colours = met.brewer("Hokusai3", n=3), type = "WUE",
                y_axis_limits = c(0, 150))
## plot normalised by stomatal density (i. e. per stoma values)
p4.2<- licorplots(c("bd", "bs", "bh"), timestamps = c(20,40), legend_title = "**Species**", 
                  legend_labels = c("*B. distachyon*", "*B. stacei*", "*B. hybridum*"),
                  stomden = c(100.3367, 102.693603, 59.5959596),
                  #remove_outliers = "yes", 
                  timeframe = 16:55, 
                  colours = met.brewer("Hokusai3", n=3), type = "WUE",
                  y_axis_limits = c(0, 0.0018))









#### Steady state physiological values and kinetics --------------------------------------------------------------

setwd("./Li-COR files")

licordata <- licorvalues(c("bd", "bs", "bh"), transition=list(c(5:18), c(21:38), c(41:58), c(61:75)), 
             label = c("*B. distachyon*", "*B. stacei*", "*B. hybridum*"), 
             colours = met.brewer("Hokusai3", n = 3), 
             remove_outliers = "yes")

licordata$transition_zone <- as.character(licordata$transition_zone)

licordata <- licordata %>% arrange(factor(transition_zone, levels = c('5:18', '21:38', '41:58', '61:75')))




## code for plot
licorvalues(c("bd", "bs", "bh"), transition=list(c(21:38), c(41:58), c(61:75)), 
            label = c("*B. distachyon*", "*B. stacei*", "*B. hybridum*"), 
            colours = met.brewer("Hokusai3", n = 3), 
            remove_outliers = "yes")







##### steady state values and kinetics normalised by stomatal density ----------------------------------------

setwd("./Li-COR files")

licordata2 <- licorvalues(c("bd", "bs", "bh"), transition=list(c(5:18), c(21:38), c(41:58), c(61:75)),
                          stomden = c(100.3367, 102.693603, 59.5959596),
                          label = c("*B. distachyon*", "*B. stacei*", "*B. hybridum*"), 
                          remove_outliers = "yes",
                          colours = met.brewer("Hokusai3", n=3))

licordata2$transition_zone <- as.character(licordata2$transition_zone)

licordata2 <- licordata2 %>% arrange(factor(transition_zone, levels = c('5:18', '21:38', '41:58', '61:75')))



### code for plot
licorvalues(c("bd", "bs", "bh"), transition=list(c(21:38), c(41:58), c(61:75)), 
            stomden = c(100.3367, 102.693603, 59.5959596),
            label = c("*B. distachyon*", "*B. stacei*", "*B. hybridum*"), 
            remove_outliers = "yes",
            colours = met.brewer("Hokusai3", n=3))











#### Stomatal density of kinetics/Li-COR leaves -----------------------------------------------------------------------

##load packages
library(tidyverse)
library(readxl)
library(agricolae)
library(MetBrewer)

### if the working directory was changed to the Li-COR files folder:
setwd("..")

### set colour palette
colourpalette <- met.brewer("Hokusai3", n=3)
#colourpalette <- c("#000000", "#666666", "#999999", "#CCCCCC", "#EEEEEE")


### load .xlsx file with data
stomden <- read_excel("density_licor-leaves-2023.xlsx", sheet=1, col_names = T)


### calculate mean stomatal density for each individual
stomden1 <- stomden %>% group_by(genotype, identifier, individual) %>% summarise(mean_dens = mean(stomatal_density),
                                                                                 nr_stomata = sum(stomata))


### test for significant different between means of individuals
##ANOVA (is there a significant difference between some of the groups?)
anovaden2 <- aov(stomden1$mean_dens ~stomden1$genotype)
summary(anovaden2)

##Tukey's test to find out which groups differ from each other
tukeyden2 <- HSD.test(anovaden2, trt="stomden1$genotype")
tukeyden2$groups


##create data frame with significance levels
sigden2 <- data.frame(genotype=c("B. distachyon", "B. stacei", "B. hybridum"),
                      significance=c("a", "a", "b"))



##reorder data
stomden1$genotype <- ordered(stomden1$genotype , levels=c("B. distachyon", "B. stacei", "B. hybridum"))



##plot stomatal density data for each genotype as box and dot plot
pden <- ggplot(stomden1, mapping=aes(x=genotype, y=mean_dens))+
  geom_boxplot()+
  geom_text(sigden2, mapping=aes(y=130, label=paste(significance)), show.legend = F)+
  geom_jitter(aes(colour=genotype), width=0.15, height = 0, alpha = 0.8, show.legend = F)+
  scale_colour_manual(values = colourpalette)+
  scale_y_continuous(limits = c(0, 130))+
  theme_classic()+
  scale_x_discrete(labels = c("*B. distachyon*", "*B. stacei*", "*B. hybridum*"))+
  theme(axis.text.x = element_markdown())+
  labs(x=NULL, y=expression(paste("Stomatal density [stomata mm"^-2, "]")))


























#### Stomatal anatomy --------------------------------------------------------------------------------------------

##load packages
library(tidyverse)
library(readxl)
library(agricolae)
library(MetBrewer)


### if the working directory was changed to the Li-COR files folder:
setwd("..")

###load excel file (.xls or .xlsx) with measurements
stomata<- read_excel("anatomy-measurements-2023.xlsx", col_names = T)
stomata <- na.omit(stomata)
stomata$length <- as.numeric(stomata$length)
stomata$width <- as.numeric(stomata$width)

##calculate length/width ratio
stomata$ratio <- stomata$length/stomata$width




### calculate means per individual
stomeans <- stomata %>% group_by(genotype, identifier, individual) %>% summarise(mean_length = mean(length), sd_length = sd(length), 
                                                                                 mean_width = mean(width), sd_width = sd(width),
                                                                                 mean_ratio = mean(ratio), sd_ratio = sd(ratio),
                                                                                 stom_nrs = length(number))






### Stomatal length -----------------------------------------

## test for significant differences between genotypes based on means of individuals
## ANOVA (is there a significant difference between some of the groups?)
anoval2 <- aov(stomeans$mean_length ~stomeans$genotype)
summary(anoval2)

##Tukey's test to find out which groups differ from each other
tukeyl2 <- HSD.test(anoval2, trt="stomeans$genotype")
##show results
tukeyl2$groups

##create data frame with significance levels for each 'genotype'
sigl2 <- data.frame(genotype=c("B. distachyon", "B. stacei", "B. hybridum"),
                    significance=c("a", "a", "b"))


##pick colours
colourpalette <- met.brewer("Hokusai3", n=3)
#colourpalette <- c("#000000", "#666666", "#999999", "#CCCCCC", "#EEEEEE")



##plot stomatal length for each genotype with means of individuals
stomeans<- stomeans %>% mutate(genotype = fct_relevel(genotype, "B. distachyon", "B. stacei", "B. hybridum"))

plen2 <- ggplot(stomeans, mapping=aes(x=genotype, y=mean_length))+
  geom_boxplot(show.legend = F)+
  geom_jitter(mapping=aes(colour=genotype), width = 0.15, alpha = 0.8, height = 0, show.legend = F)+
  geom_text(sigl2, mapping=aes(y=40, label=paste(significance)))+
  scale_colour_manual(values = colourpalette)+
  scale_y_continuous(limits = c(0, 40))+
  theme_classic()+
  scale_x_discrete(labels = c("*B. distachyon*", "*B. stacei*", "*B. hybridum*"))+
  theme(axis.text.x = element_markdown())+
  labs(x=NULL, y="Stomatal length [µm]")







### Stomatal width ---------------------------------------------

## test for significant differences between the genotypes based on the means per individual
## ANOVA (is there a significant difference between some of the groups?)
anovaw2 <- aov(stomeans$mean_width ~stomeans$genotype)
summary(anovaw2)

## Tukey's test to find out which groups differ from each other
tukeyw2 <- HSD.test(anovaw2, trt="stomeans$genotype")
## show results
tukeyw2$groups

## create data frame with significance levels for each 'genotype'
sigw2 <- data.frame(genotype=c("B. distachyon", "B. stacei", "B. hybridum"),
                    significance=c("a", "a", "b"))


## pick colours
colourpalette <- met.brewer("Hokusai3", n=3)


## plot width for all genotypes based on means per individual
stomeans<- stomeans %>% mutate(genotype = fct_relevel(genotype, "B. distachyon", "B. stacei", "B. hybridum"))

pwid2 <- ggplot(stomeans, mapping=aes(x=genotype, y=mean_width))+
  geom_boxplot()+
  geom_jitter(aes(colour=genotype), width = 0.15, height = 0, alpha = 0.8, show.legend = F)+
  geom_text(sigw2, mapping=aes(y=8, label=paste(significance)))+
  scale_colour_manual(values = colourpalette)+
  scale_y_continuous(limits = c(0, 8))+
  theme_classic()+
  scale_x_discrete(labels = c("*B. distachyon*", "*B. stacei*", "*B. hybridum*"))+
  theme(axis.text.x = element_markdown())+
  labs(x=NULL, y="Stomatal width [µm]")








### Stomatal length/width ratio ------------------------------------------

## test for significant differences based on means per individual
## ANOVA (is there a significant difference between some of the groups?)
anovar2 <- aov(stomeans$mean_ratio ~stomeans$genotype)
summary(anovar2)

## Tukey's test to find out which groups differ from each other
tukeyr2 <- HSD.test(anovar2, trt="stomeans$genotype")
## show results
tukeyr2$groups

## create data frame with significance levels for each 'genotype'
sigr2 <- data.frame(genotype=c("B. distachyon", "B. stacei", "B. hybridum"),
                    significance=c("a", "a", "b"))



## plot length/width ratio for each genotype based on means per individual
stomeans<- stomeans %>% mutate(genotype = fct_relevel(genotype, "B. distachyon", "B. stacei", "B. hybridum"))

prat2 <- ggplot(stomeans, mapping=aes(x=factor(genotype), y=mean_ratio))+
  geom_boxplot()+
  geom_jitter(aes(colour=genotype), width = 0.15, height = 0, alpha = 0.8, show.legend = F)+
  geom_text(sigr2, mapping=aes(y=5, label=paste(significance)))+
  scale_colour_manual(values = colourpalette)+
  theme_classic()+
  scale_x_discrete(labels = c("*B. distachyon*", "*B. stacei*", "*B. hybridum*"))+
  scale_y_continuous(limits = c(0,5))+
  theme(axis.text.x = element_markdown())+
  labs(x=NULL, y="Stomatal length to width ratio")








#### test for correlation between stomatal density and stomatal anatomy ---------------------------------------------
library(ggpubr)
#devtools::install_github('smin95/smplot2')
library(smplot2)
library(MetBrewer)


### merge data objects for stomatal anatomy and stomatal density data (see above)
stomana <- merge(stomeans, stomden1)
stomana<- stomana %>% mutate(genotype = fct_relevel(genotype, "B. distachyon", "B. stacei", "B. hybridum"))


pdenlen <- ggplot(stomana, aes(x = mean_length, y = mean_dens))+
  geom_point(aes(colour = genotype), size = 3, alpha = 0.7)+
  sm_statCorr(corr_method = "pearson", colour = "darkgray",
              separate_by = '\n',
              label_x = 21, label_y = 55)+
  theme_classic()+
  scale_colour_manual(values = met.brewer("Hokusai3", n=3), 
                      labels = c("*B. distachyon*", "*B. stacei*", "*B. hybridum*"))+
  guides(colour = guide_legend("Species"))+
  theme(legend.position = c(.8, .8),
        legend.background = element_rect(fill = "#e5e4e2"),
        legend.text = element_markdown())+
  labs(x= "Stomatal length [µm]", y=expression(paste("Stomatal density [stomata mm"^-2, "]")))


pdenwid <- ggplot(stomana, aes(x = mean_width, y = mean_dens))+
  geom_point(aes(colour = genotype), size = 3, alpha = 0.7)+
  sm_statCorr(corr_method = "pearson", colour = "darkgray",
              separate_by = '\n',
              label_x = 6.2, label_y = 55)+
  theme_classic()+
  scale_colour_manual(values = met.brewer("Hokusai3", n=3), 
                      labels = c("*B. distachyon*", "*B. stacei*", "*B. hybridum*"))+
  guides(colour = guide_legend("Species"))+
  theme(legend.position = c(.8, .8),
        legend.background = element_rect(fill = "#e5e4e2"),
        legend.text = element_markdown())+
  labs(x= "Stomatal width [µm]", y=expression(paste("Stomatal density [stomata mm"^-2, "]")))
















#### calculate anatomical gsmax based on measurements from leaves used for kinetics -------------------------------

### load packages
library(tidyverse)
library(readxl)
library(agricolae)
library(MetBrewer)
library(ggtext)


### if the working directory was changed to the Li-COR files folder:
setwd("..")

### load excel file (.xls or .xlsx) with measurements
stomata<- read_excel("anatomy-measurements-2023.xlsx", col_names = T)
stomata$length <- as.numeric(stomata$length)
stomata$width <- as.numeric(stomata$width)
stomata <- na.omit(stomata)


### calculate anatomical gsmax
## pore area
stomata$amax <- 0.9*stomata$PL*stomata$PW
## pore depth
stomata$l <- 0.5*stomata$width
## constant for diffusivity of water in air
d <- 0.0000249
## constant for molar volume of air
v <- 0.024464
stomata$anato_gsmax <- (d*stomata$`SD (from \"density_licor-leaves-2023\")`*stomata$amax)/(v*(stomata$l+pi/2*sqrt(stomata$amax/pi)))



### calculate means per individual
stomeans <- stomata %>% group_by(genotype, identifier, individual) %>% summarise(mean_length = mean(length),
                                                                                 mean_width = mean(width),
                                                                                 mean_ratio = mean(ratio),
                                                                                 mean_anatomical = mean(anato_gsmax),
                                                                                 sd_anatomical = sd(anato_gsmax))


### test for significant differences between genotypes based on means of individuals
### ANOVA (is there a significant difference between some of the groups?)
anovaa2 <- aov(stomeans$mean_anatomical ~stomeans$genotype)
summary(anovaa2)

### Tukey's test to find out which groups differ from each other
tukeya2 <- HSD.test(anovaa2, trt="stomeans$genotype")
### show results
tukeya2$groups

### create data frame with significance levels for each 'genotype'
siga2 <- data.frame(genotype=c("B. distachyon", "B. stacei", "B. hybridum"),
                    significance=c("ab", "a", "b"))



### plot anatomical gsmax for each genotype (means of individuals)
colourpalette <- met.brewer("Hokusai3", n=3)

stomeans <- stomeans %>% mutate(genotype = fct_relevel(genotype, "B. distachyon", "B. stacei", "B. hybridum"))

panakin2 <- ggplot(stomeans, mapping=aes(x=genotype, y=mean_anatomical))+
  geom_boxplot(show.legend = F)+
  geom_jitter(mapping=aes(colour=genotype), width = 0.15, alpha = 0.8, height = 0, show.legend = F)+
  geom_text(siga2, mapping=aes(y=0.6, label=paste(significance)))+
  scale_colour_manual(values = colourpalette)+
  scale_y_continuous(limits = c(0, 0.6))+
  theme_classic()+
  scale_x_discrete(labels = c("*B. distachyon*", "*B. stacei*", "*B. hybridum*"))+
  theme(axis.text.x = element_markdown(),
        axis.title.y = element_markdown())+
  labs(x=NULL, y="Anatomical *g*<sub>S</sub>max [mol m<sup>-2</sup> s<sup>-1</sup>]")


### values
data2 <- stomeans %>% group_by(genotype) %>% summarise(mean_anato_gsmax=mean(mean_anatomical),
                                                       sd_anato_gsmax = sd(mean_anatomical))

data2








#### calculate anatomical gsmax per stoma (divided by stomatal density*1000) ---------------------------------------

stomata$anato_gsmax_norm <- stomata$anato_gsmax/(stomata$`SD (from \"density_licor-leaves-2023\")`*1000)

stomeans <- stomata %>% group_by(genotype, identifier, individual) %>% summarise(mean_length = mean(length),
                                                                                 mean_width = mean(width),
                                                                                 mean_ratio = mean(ratio),
                                                                                 mean_anatomical = mean(anato_gsmax),
                                                                                 sd_anatomical = sd(anato_gsmax),
                                                                                 mean_anato_norm = mean(anato_gsmax_norm),
                                                                                 sd_anato_norm = sd(anato_gsmax_norm))


### test for significant differences between genotypes based on means per individual
### ANOVA (is there a significant difference between some of the groups?)
anovaan2 <- aov(stomeans$mean_anato_norm ~stomeans$genotype)
summary(anovaan2)

### Tukey's test to find out which groups differ from each other
tukeyan2 <- HSD.test(anovaan2, trt="stomeans$genotype")
### show results
tukeyan2$groups

### create data frame with significance levels for each 'genotype'
sigan2 <- data.frame(genotype=c("B. distachyon", "B. stacei", "B. hybridum"),
                     significance=c("a", "a", "b"))



### plot anatomical gsmax for each genotype based on means per individual
colourpalette <- met.brewer("Hokusai3", n=3)

stomeans <- stomeans %>% mutate(genotype = fct_relevel(genotype, "B. distachyon", "B. stacei", "B. hybridum"))

panankin2 <- ggplot(stomeans, mapping=aes(x=genotype, y=mean_anato_norm))+
  geom_boxplot(show.legend = F)+
  geom_jitter(mapping=aes(colour=genotype), width = 0.15, alpha = 0.8, height = 0, show.legend = F)+
  geom_text(sigan2, mapping=aes(y=0.000008, label=paste(significance)))+
  scale_colour_manual(values = colourpalette)+
  scale_y_continuous(limits = c(0, 0.000008), breaks=c(0, 0.000002, 0.000004, 0.000006, 0.000008), labels = c(0, 2, 4, 6, 8))+
  theme_classic()+
  scale_x_discrete(labels = c("*B. distachyon*", "*B. stacei*", "*B. hybridum*"))+
  theme(axis.text.x = element_markdown(),
        axis.title.y = element_markdown())+
  labs(x=NULL, y="Anatomical *g*<sub>S</sub>max [µmol stoma<sup>-1</sup> s<sup>-1</sup>]")



### values
data2 <- stomeans %>% group_by(genotype) %>% summarise(mean_anato_gsmax_norm=mean(mean_anato_norm),
                                                       sd_anato_gsmax_norm=sd(mean_anato_norm))
data2