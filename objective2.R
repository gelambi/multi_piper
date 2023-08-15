##########################################################################################################
### Multiscale variability in nutrients and secondary metabolites in a bat-dispersed neotropical fruit ###
##########################################################################################################
rm(list=ls()) # removes objects in workspace
library(dplyr)
library(tidyverse)
library(devtools)
library(glmmTMB)
library(RColorBrewer)
library(viridis)
library(viridisLite)
library(corrplot)
library(MetBrewer)
library(effects)
library(ggpubr)
library(parameters)
library(ggeffects)
library(gridExtra)
library(psych)
library(vegan)
library(emmeans)
library(ggfortify)
library(MuMIn)
library(sjPlot)
library(AICcmodavg)
library(ggridges)
library(performance)
library(ggcorrplot)

##################################################################################################################
### QUESTION 2 (Q2): Do the patterns of association between nutrients and defensive metabolites at the inter- ###
### and intraindividual scales support the removal-rate or nutrient-toxin titration/relative risk hypotheses? ###
##################################################################################################################

#################
### READ DATA ###
#################
data <-  read.csv("data_FINAL.csv")
head(data) # the letters are different alkenylphenols 
data$tree <- as.factor(data$tree)
levels(data$tree) # There is no tree 4, so I will rename treeID for consistency 
data$tree <- factor(data$tree, levels = c("1","2","3","5","6","7","8","9","10","11"), 
                    labels = c("1","2","3","4","5","6","7","8","9","10"))
data$fruitID <- as.factor(data$fruitID)
head(data)

########################
### Correlation plot ###
########################
correlation_data_names <- data %>% 
  rename(
    "Total phenolics" = totalphenolics,
    "Glucose" = dglucose, 
    "Fructose" = fructose, 
    "Total proteins" = proteins,
    "Total alkenylphenols" = alkenylphenols
  )

just_traits <- correlation_data_names[, 4:18]
cor <- cor(just_traits, method = "spearman")
p.mat <- cor.mtest(cor, method = "spearman")

# Create heatmap
p <- ggplot(data = as.data.frame(cor_df), aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() +
  theme_classic() +
  scale_fill_gradient2(low = "#240154FF", mid = "white", high = "#C2DA37FF",
                       midpoint = 0, limits = c(-1, 1), name = "Correlation\nCoefficient") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, color = "black"),
        axis.text.y = element_text(color = "black"),
        legend.position = "right") +
  labs(title = " ",
       fill = "Correlation\nCoefficient",
       x = "", y = "")
p
# Add correlation significance levels to heatmap
figure_5 <- p + annotate("text", x = expand.grid(seq(1:ncol(cor)), seq(1:ncol(cor)))[as.vector(p.mat$p < 0.001), 1],
             y = expand.grid(seq(1:ncol(cor)), seq(1:ncol(cor)))[as.vector(p.mat$p < 0.001), 2],
             label = "***", size = 3) +
  annotate("text", x = expand.grid(seq(1:ncol(cor)), seq(1:ncol(cor)))[as.vector(p.mat$p >= 0.001 & p.mat$p < 0.01), 1],
           y = expand.grid(seq(1:ncol(cor)), seq(1:ncol(cor)))[as.vector(p.mat$p >= 0.001 & p.mat$p < 0.01), 2],
           label = "**", size = 3) +
  annotate("text", x = expand.grid(seq(1:ncol(cor)), seq(1:ncol(cor)))[as.vector(p.mat$p >= 0.01 & p.mat$p < 0.05), 1],
           y = expand.grid(seq(1:ncol(cor)), seq(1:ncol(cor)))[as.vector(p.mat$p >= 0.01 & p.mat$p < 0.05), 2],
           label = "*", size = 3)
ggsave(file="Figure_5.jpg", 
       plot= figure_5,
       width=6,height=5,units="in",dpi=300)

#############################
### INTRAINDIVIDUAL LEVEL ###
#############################
head(data)
data$total_sugars <- data$dglucose + data$fructose
### GLMMs with family = beta

### transform percentage data in proportions: 
data$sugars_beta <- (data$total_sugars)/100
data$proteins_beta <- (data$proteins)/100
data$totalphenolics_beta <- (data$totalphenolics)/100
data$alkenylphenols_beta <- (data$alkenylphenols)/100

# Fit a GLMM with total defensive metabolites as the response variables, and total sugars and proteins as predictor variables.

### Phenolics 
glmm_phenolics_beta <- glmmTMB(totalphenolics_beta ~ sugars_beta + proteins_beta + (1|tree), data = data,  family = beta_family(link="logit"))
summary(glmm_phenolics_beta)
diagnose(glmm_phenolics_beta)
shapiro.test(resid(glmm_phenolics_beta)) # Check for normality of residuals using the Shapiro-Wilk test: normal 
par_phenolics <- parameters(glmm_phenolics_beta) # Extract the estimated model parameters into a data frame and save it as a CSV file.
write.csv(par_phenolics, file = "parameters_phenolics_intra.csv")
plot(allEffects(glmm_phenolics_beta)) # Create a plot of the estimated marginal effects of each predictor variable on the response variable using the allEffects function from the effects package.
r2(glmm_phenolics_beta) 

### Alkenylphenols 
glmm_alkenylphenols_beta <- glmmTMB(alkenylphenols_beta ~ sugars_beta + proteins_beta + (1|tree), data = data, family = beta_family(link="logit"))
summary(glmm_alkenylphenols_beta)
diagnose(glmm_alkenylphenols_beta)
shapiro.test(resid(glmm_alkenylphenols_beta)) # normal
par_alkenylphenols <- parameters(glmm_alkenylphenols_beta)
write.csv(par_alkenylphenols, file = "parameter_alkenylphenols_intra.csv")
plot(allEffects(glmm_alkenylphenols_beta))
r2(glmm_alkenylphenols_beta)

#### Plot PHENOLICS 

### sugars (x) versus phenolics (y)
new_data_sugars <- data.frame(sugars_beta = data$sugars_beta,
                               proteins_beta = mean(data$proteins_beta),  # Set proteins_beta to its mean value
                               tree = data$tree) # Set tree as factor with same levels as in the original data
new_data_sugars$predicted <- predict(glmm_phenolics_beta, newdata = new_data_sugars, type = "response")
#new_data_sugars$lower_ci <- predict(glmm_phenolics_beta, re.form = NA, newdata = new_data_sugars, type = "response", se.fit = TRUE)$fit - 1.96 * predict(glmm_phenolics_beta,re.form = NA, newdata = new_data_sugars, type = "response", se.fit = TRUE)$se.fit
#new_data_sugars$upper_ci <- predict(glmm_phenolics_beta, re.form = NA, newdata = new_data_sugars,  type = "response", se.fit = TRUE)$fit + 1.96 * predict(glmm_phenolics_beta, re.form = NA, newdata = new_data_sugars, type = "response", se.fit = TRUE)$se.fit

plot_sugars_phenolics <- ggplot() + 
  theme_classic(base_size = 13) +
  #geom_ribbon(data = new_data_sugars, aes(x = sugars_beta, ymin = lower_ci, ymax = upper_ci), fill = "grey", alpha = 0.3) +
  geom_line(data = new_data_sugars, aes(x = sugars_beta, y = predicted, group = tree, color = tree), linetype = 2) + 
  geom_point(data = data, aes(x = sugars_beta, y = totalphenolics_beta, color = tree), size = 3) +
  scale_color_viridis(option = "D", discrete=TRUE, name = "Plant ID", labels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10")) + 
  ylab ("Total phenolics\n(mg/mg DW)") +
  xlab (" ") +
  theme(legend.position = "none") 
plot_sugars_phenolics

new_data_proteins_beta <- data.frame(proteins_beta = data$proteins_beta,
                                     sugars_beta = mean(data$sugars_beta),
                                     tree = data$tree) 
new_data_proteins_beta$predicted <- predict(glmm_phenolics_beta, newdata = new_data_proteins_beta, type = "response")
#new_data_proteins_beta$lower_ci <- predict(glmm_phenolics_beta, newdata = new_data_proteins_beta, type = "response", se.fit = TRUE)$fit - 1.96 * predict(glmm_phenolics_beta,re.form = NA, newdata = new_data_proteins_beta, type = "response", se.fit = TRUE)$se.fit
#new_data_proteins_beta$upper_ci <- predict(glmm_phenolics_beta, newdata = new_data_proteins_beta,  type = "response", se.fit = TRUE)$fit + 1.96 * predict(glmm_phenolics_beta, re.form = NA, newdata = new_data_proteins_beta, type = "response", se.fit = TRUE)$se.fit

plot_proteins_beta_phenolics <- ggplot() + 
  theme_classic(base_size = 13) +
  #geom_ribbon(data = new_data_proteins_beta, aes(x = proteins_beta, ymin = lower_ci, ymax = upper_ci), fill = "grey", alpha = 0.3) +
  geom_line(data = new_data_proteins_beta, aes(x = proteins_beta, y = predicted, group = tree, color = tree)) + 
  geom_point(data = data, aes(x = proteins_beta, y = totalphenolics_beta, color = tree), size = 3) +
  scale_color_viridis(option = "D", discrete=TRUE, name = "Plant ID", labels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10")) + 
  ylab (" ") +
  xlab (" ") +
  theme(legend.position = "none")
plot_proteins_beta_phenolics

### effect size 
# Proteins 0.004, Phenolics 0.006
# Proteins 0.008, Phenolics 0.008

0.006 - 0.008 # 0.002
(0.002*100)/0.006 # 33.3

### plot ALKENYLPHENOLS

### sugars (x) versus alkenylphenols (y)
new_data_sugars <- data.frame(sugars_beta = data$sugars_beta,
                              proteins_beta = mean(data$proteins_beta),  # Set proteins_beta to its mean value
                              tree = data$tree) # Set tree as factor with same levels as in the original data
new_data_sugars$predicted <- predict(glmm_alkenylphenols_beta, newdata = new_data_sugars, type = "response")
#new_data_sugars$lower_ci <- predict(glmm_alkenylphenols_beta, re.form = NA, newdata = new_data_sugars, type = "response", se.fit = TRUE)$fit - 1.96 * predict(glmm_alkenylphenols_beta,re.form = NA, newdata = new_data_sugars, type = "response", se.fit = TRUE)$se.fit
#new_data_sugars$upper_ci <- predict(glmm_alkenylphenols_beta, re.form = NA, newdata = new_data_sugars,  type = "response", se.fit = TRUE)$fit + 1.96 * predict(glmm_alkenylphenols_beta, re.form = NA, newdata = new_data_sugars, type = "response", se.fit = TRUE)$se.fit

plot_sugars_alkenylphenols <- ggplot() + 
  theme_classic(base_size = 13) +
  #geom_ribbon(data = new_data_sugars, aes(x = sugars_beta, ymin = lower_ci, ymax = upper_ci), fill = "grey", alpha = 0.3) +
  geom_line(data = new_data_sugars, aes(x = sugars_beta, y = predicted, group = tree, color = tree), linetype = 2) + 
  geom_point(data = data, aes(x = sugars_beta, y = alkenylphenols_beta, color = tree), size = 3) +
  scale_color_viridis(option = "D", discrete=TRUE, name = "Plant ID", labels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10")) + 
  ylab ("Total alkenylphenols\n(mg/mg DW)") +
  xlab ("Total sugars\n(mg/mg DW)") +
  theme(legend.position = "none") 
plot_sugars_alkenylphenols

### proteins_beta (x) versus alkenylphenols (y)
new_data_proteins_beta <- data.frame(proteins_beta = data$proteins_beta,
                                     sugars_beta = mean(data$sugars_beta),
                                     tree = data$tree) # Set tree as factor with same levels as in the original data
new_data_proteins_beta$predicted <- predict(glmm_alkenylphenols_beta, newdata = new_data_proteins_beta, type = "response")
#new_data_proteins_beta$lower_ci <- predict(glmm_alkenylphenols_beta, re.form = NA, newdata = new_data_proteins_beta, type = "response", se.fit = TRUE)$fit - 1.96 * predict(glmm_alkenylphenols_beta,re.form = NA, newdata = new_data_proteins_beta, type = "response", se.fit = TRUE)$se.fit
#new_data_proteins_beta$upper_ci <- predict(glmm_alkenylphenols_beta, re.form = NA, newdata = new_data_proteins_beta,  type = "response", se.fit = TRUE)$fit + 1.96 * predict(glmm_alkenylphenols_beta, re.form = NA, newdata = new_data_proteins_beta, type = "response", se.fit = TRUE)$se.fit

plot_proteins_beta_alkenylphenols <- ggplot() + 
  theme_classic(base_size = 13) +
  #geom_ribbon(data = new_data_proteins_beta, aes(x = proteins_beta, ymin = lower_ci, ymax = upper_ci), fill = "grey", alpha = 0.3) +
  geom_line(data = new_data_proteins_beta, aes(x = proteins_beta, y = predicted, group = tree, color = tree), linetype = 2) + 
  geom_point(data = data, aes(x = proteins_beta, y = alkenylphenols_beta, color = tree), size = 3) +
  scale_color_viridis(option = "D", discrete=TRUE, name = "Plant ID", labels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10")) + 
  ylab (" ") +
  xlab ("Total proteins\n(mg/mg DW)") +
  theme(legend.position = "none")
plot_proteins_beta_alkenylphenols

intradividual_totalsugars <- ggarrange(plot_sugars_phenolics,
                                       plot_proteins_beta_phenolics,
                                       plot_sugars_alkenylphenols,
                                       plot_proteins_beta_alkenylphenols,
                                       ncol = 2,
                                       nrow = 2,
                                       common.legend = TRUE,
                                       legend="bottom",
                                       align = "hv")

intradividual_totalsugars
ggsave(file="Figure_6_randomeffect.jpg", 
       plot= intradividual_totalsugars,
       width=8,height=8,units="in",dpi=300)

#######################################################
### INTERINDIVIDUAL LEVEL: JUST SUGARS AND PROTEINS ###
#######################################################
data_average <- read.csv("data_average.csv")
data_average$sugars <- data_average$dglucose + data_average$fructose
data_average$tree <- as.factor(data_average$tree)
head(data_average)

data_average$sugars_beta <- (data_average$dglucose)/100
data_average$proteins_beta <- (data_average$proteins)/100
data_average$totalphenolics_beta <- (data_average$totalphenolics)/100
data_average$alkenylphenols_beta <- (data_average$alkenylphenols)/100

### PHENOLICS
glmm_phenolics_beta_intra <- glmmTMB(totalphenolics_beta ~ sugars_beta + proteins_beta, data = data_average, beta_family(link="logit"))
summary(glmm_phenolics_beta_intra)
diagnose(glmm_phenolics_beta_intra)
shapiro.test(resid(glmm_phenolics_beta_intra)) # normal
par_phenolics <- parameters(glmm_phenolics_beta_intra)
write.csv(par_phenolics, file = "parameters_phenolics_inter.csv")
plot(allEffects(glmm_phenolics_beta_intra))
r2(glmm_phenolics_beta_intra) 

### ALKENYLPHENOLS
glmm_alkenylphenols_beta_inter <- glmmTMB(alkenylphenols_beta ~ sugars_beta + proteins_beta, data = data_average, beta_family(link="logit"))
summary(glmm_alkenylphenols_beta_inter)
diagnose(glmm_alkenylphenols_beta_inter)
shapiro.test(resid(glmm_alkenylphenols_beta_inter))
diagnose(glmm_alkenylphenols_beta_inter) 
par_alkenylphenols_inter <- parameters(glmm_alkenylphenols_beta_inter)
write.csv(par_alkenylphenols_inter, file = "par_alkenylphenols_inter.csv")
plot(allEffects(glmm_alkenylphenols_beta_inter))
r2(glmm_alkenylphenols_beta_inter)

### plots, PHENOLICS

### glucose (x) versus phenolics (y)
new_data_average_glucose <- data.frame(sugars_beta = data_average$sugars_beta,
                                       proteins_beta = mean(data_average$proteins_beta))  # Set proteins_beta to its mean value
new_data_average_glucose$predicted <- predict(glmm_phenolics_beta_intra, re.form = NA, newdata = new_data_average_glucose, type = "response")
new_data_average_glucose$lower_ci <- predict(glmm_phenolics_beta_intra, re.form = NA, newdata= new_data_average_glucose, type = "response", se.fit = TRUE)$fit - 1.96 * predict(glmm_phenolics_beta_intra,re.form = NA, newdata = new_data_average_glucose, type = "response", se.fit = TRUE)$se.fit
new_data_average_glucose$upper_ci <- predict(glmm_phenolics_beta_intra, re.form = NA, newdata = new_data_average_glucose,  type = "response", se.fit = TRUE)$fit + 1.96 * predict(glmm_phenolics_beta_intra, re.form = NA, newdata = new_data_average_glucose, type = "response", se.fit = TRUE)$se.fit

plot_sugars_phenolics <- ggplot() + 
  theme_classic(base_size = 13) +
  geom_ribbon(data = new_data_average_glucose, aes(x = sugars_beta, ymin = lower_ci, ymax = upper_ci), fill = "grey", alpha = 0.3) +
  geom_line(data = new_data_average_glucose, aes(x = sugars_beta, y = predicted), color = "black",  linetype = 2) + 
  geom_point(data = data_average, aes(x = sugars_beta, y = totalphenolics_beta, color = tree), size = 3) +
  scale_color_viridis(option = "D", discrete=TRUE, name = "Plant ID", labels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10")) + 
  ylab ("Total phenolics\n(mg/mg DW)") +
  xlab (" ") +
  theme(legend.position = "none")
plot_sugars_phenolics

### proteins_beta (x) versus phenolics (y)
new_data_average_proteins_beta <- data.frame(proteins_beta = data_average$proteins_beta,
                                             sugars_beta = mean(data_average$sugars_beta))
new_data_average_proteins_beta$predicted <- predict(glmm_phenolics_beta_intra, re.form = NA, newdata = new_data_average_proteins_beta, type = "response")
new_data_average_proteins_beta$lower_ci <- predict(glmm_phenolics_beta_intra, re.form = NA, newdata = new_data_average_proteins_beta, type = "response", se.fit = TRUE)$fit - 1.96 * predict(glmm_phenolics_beta_intra,re.form = NA, newdata = new_data_average_proteins_beta, type = "response", se.fit = TRUE)$se.fit
new_data_average_proteins_beta$upper_ci <- predict(glmm_phenolics_beta_intra, re.form = NA, newdata = new_data_average_proteins_beta,  type = "response", se.fit = TRUE)$fit + 1.96 * predict(glmm_phenolics_beta_intra, re.form = NA, newdata = new_data_average_proteins_beta, type = "response", se.fit = TRUE)$se.fit

plot_proteins_beta_phenolics <- ggplot() + 
  theme_classic(base_size = 13) +
  geom_ribbon(data = new_data_average_proteins_beta, aes(x = proteins_beta, ymin = lower_ci, ymax = upper_ci), fill = "grey", alpha = 0.3) +
  geom_line(data = new_data_average_proteins_beta, aes(x = proteins_beta, y = predicted), color = "black", linetype = 2) + 
  geom_point(data = data_average, aes(x = proteins_beta, y = totalphenolics_beta, color = tree), size = 3) +
  scale_color_viridis(option = "D", discrete=TRUE, name = "Plant ID", labels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10")) + 
  ylab (" ") +
  xlab (" ") +
  theme(legend.position = "none")
plot_proteins_beta_phenolics

### effect size 
head(new_data_proteins_beta)
new_data_proteins_beta$predicted_intraspecific <- predict(glmm_phenolics_beta_intra,newdata = new_data_proteins_beta, type = "response")

# Proteins 0.004, Phenolics 0.005
# Proteins 0.006, Phenolics 0.006

0.006 - 0.005 # 0.001
(0.001*100)/0.006 # 16.66667

### plots, ALKENYLPHENOLS

### glucose (x) versus alkenylphenols (y)
new_data_sugars <- data.frame(sugars_beta = data_average$sugars_beta,
                               proteins_beta = mean(data_average$proteins_beta))  # Set proteins_beta to its mean value
new_data_sugars$predicted <- predict(glmm_alkenylphenols_beta_inter, re.form = NA, newdata = new_data_sugars, type = "response")
new_data_sugars$lower_ci <- predict(glmm_alkenylphenols_beta_inter, re.form = NA, newdata = new_data_sugars, type = "response", se.fit = TRUE)$fit - 1.96 * predict(glmm_alkenylphenols_beta_inter,re.form = NA, newdata = new_data_sugars, type = "response", se.fit = TRUE)$se.fit
new_data_sugars$upper_ci <- predict(glmm_alkenylphenols_beta_inter, re.form = NA, newdata = new_data_sugars,  type = "response", se.fit = TRUE)$fit + 1.96 * predict(glmm_alkenylphenols_beta_inter, re.form = NA, newdata = new_data_sugars, type = "response", se.fit = TRUE)$se.fit

plot_sugars_alkenylphenols <- ggplot() + 
  theme_classic(base_size = 13) +
  geom_ribbon(data = new_data_sugars, aes(x = sugars_beta, ymin = lower_ci, ymax = upper_ci), fill = "grey", alpha = 0.3) +
  geom_line(data = new_data_sugars, aes(x = sugars_beta, y = predicted), color = "black",  linetype = 2) + 
  geom_point(data = data_average, aes(x = sugars_beta, y = alkenylphenols_beta, color = tree), size = 3) +
  scale_color_viridis(option = "D", discrete=TRUE, name = "Plant ID", labels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10")) + 
  ylab ("Total alkenylphenols\n(mg/mg DW)") +
  xlab ("Total sugars\n(mg/mg DW)") +
  theme(legend.position = "none")
plot_sugars_alkenylphenols

### proteins_beta (x) versus alkenylphenols (y)
new_data_proteins_beta <- data.frame(proteins_beta = data_average$proteins_beta,
                                     sugars_beta = mean(data_average$sugars_beta))  # Set proteins_beta to its mean value
new_data_proteins_beta$predicted <- predict(glmm_alkenylphenols_beta_inter, re.form = NA, newdata = new_data_proteins_beta, type = "response")
new_data_proteins_beta$lower_ci <- predict(glmm_alkenylphenols_beta_inter, re.form = NA, newdata = new_data_proteins_beta, type = "response", se.fit = TRUE)$fit - 1.96 * predict(glmm_alkenylphenols_beta_inter,re.form = NA, newdata = new_data_proteins_beta, type = "response", se.fit = TRUE)$se.fit
new_data_proteins_beta$upper_ci <- predict(glmm_alkenylphenols_beta_inter, re.form = NA, newdata = new_data_proteins_beta,  type = "response", se.fit = TRUE)$fit + 1.96 * predict(glmm_alkenylphenols_beta_inter, re.form = NA, newdata = new_data_proteins_beta, type = "response", se.fit = TRUE)$se.fit

plot_proteins_beta_alkenylphenols <- ggplot() + 
  theme_classic(base_size = 13) +
  geom_ribbon(data = new_data_proteins_beta, aes(x = proteins_beta, ymin = lower_ci, ymax = upper_ci), fill = "grey", alpha = 0.3) +
  geom_line(data = new_data_proteins_beta, aes(x = proteins_beta, y = predicted), color = "black", linetype = 2) + 
  geom_point(data = data_average, aes(x = proteins_beta, y = alkenylphenols_beta, color = tree), size = 3) +
  scale_color_viridis(option = "D", discrete=TRUE, name = "Plant ID", labels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10")) + 
  ylab (" ") +
  xlab ("Total proteins\n(mg/mg DW)") +
  theme(legend.position = "none")
plot_proteins_beta_alkenylphenols

new_data_proteins_beta$predicted <- predict(glmm_alkenylphenols_beta_inter, re.form = NA, newdata = new_data_proteins_beta, type = "response")

# Proteins 0.004, Alk 0.2
# Proteins 0.008, Alk 0.05

0.2 - 0.05 # 0.15
(0.15*100)/0.05 # 300% 

interindivual_totalsugars <- ggarrange(plot_sugars_phenolics,
                                       plot_proteins_beta_phenolics,
                                       plot_sugars_alkenylphenols,
                                       plot_proteins_beta_alkenylphenols,
                                       ncol = 2,
                                       nrow = 2,
                                       common.legend = TRUE,
                                       legend="bottom",
                                       align = "hv")

interindivual_totalsugars
ggsave(file="Figure_7.jpg", 
       plot= interindivual_totalsugars,
       width=8,height=8,units="in",dpi=300)

#### k-fold model validation ####

# Set the number of folds for cross-validation
k <- 10

# Create an empty list to store the predictions and validation results
fold_predictions <- list()
validation_results <- numeric(k)

# Set the seed for reproducibility
set.seed(123)

# Create the fold indices using createDataPartition
fold_indices <- createDataPartition(data$totalphenolics_beta, p = 0.8, times = k)

# Define the formula for the model
formula <- totalphenolics_beta ~ sugars_beta + proteins_beta + (1|tree)

# Perform k-fold cross-validation
for (i in 1:k) {
  # Create training and testing datasets for the current fold
  fold_train_data <- data[fold_indices[[i]], ]
  fold_test_data <- data[-fold_indices[[i]], ]
  
  # Fit the model using the training data for the current fold
  model <- glmmTMB(formula, data = fold_train_data, family = beta_family(link = "logit"))
  
  # Make predictions on the testing set
  predictions <- predict(model, newdata = fold_test_data, type = "response")
  
  # Calculate RMSE
  rmse <- sqrt(mean((fold_test_data$totalphenolics_beta - predictions)^2, na.rm = TRUE))
  
  # Store the predictions and validation result for this fold
  fold_predictions[[i]] <- predictions
  validation_results[i] <- rmse
}

summary(model)
print(model)
# Combine the fold predictions into a single vector
all_predictions <- unlist(fold_predictions)

# Create a vector of the true response values
true_values <- rep(data$totalphenolics_beta[unlist(fold_indices)], each = nrow(fold_test_data))

# Calculate the overall RMSE
overall_rmse <- sqrt(mean((true_values - all_predictions)^2, na.rm = TRUE))

# Print the validation results and overall RMSE
print(validation_results)
print(overall_rmse)

data$totalphenolics
########################################
########################################
########################################
########################################

