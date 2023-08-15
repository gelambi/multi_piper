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
library(moments)
library(ggforce)
library(gridExtra) 

##################################################################################################
### (1) Describe the multi-scale patterns of variation for nutrients and defensive metabolites ### 
##################################################################################################

#################
### READ DATA ###
#################
data <-  read.csv("data_FINAL.csv")
head(data) # the letters are different alkenylphenols 
data$tree <- as.factor(data$tree)
levels(data$tree)
data$tree <- factor(data$tree, levels = c("1","2","3","5","6","7","8","9","10","11"), labels = c("1","2","3","4","5","6","7","8","9","10"))
data$fruitID <- as.factor(data$fruitID)
str(data)

### distribution of the chemical variables ###
shapiro.test(data$totalphenolics)
shapiro.test(data$dglucose)
shapiro.test(data$fructose) # non-normal
shapiro.test(data$proteins)
shapiro.test(data$alkenylphenols) # non-normal
##########################
### SUMMARY STATISTICS ###
##########################

### Calculate different metrics per tree
summary <- data %>%
  group_by(tree) %>%
  filter %>%
  summarise(
    mean_dglucose = mean(dglucose),
    mean_fructose = mean(fructose),
    mean_proteins = mean(proteins),
    mean_alkenylphenols = mean(alkenylphenols),
    mean_totalphenolics = mean(totalphenolics),
    sd_dglucose = sd(dglucose),
    sd_fructose = sd(fructose),
    sd_proteins = sd(proteins),
    sd_alkenylphenols = sd(alkenylphenols),
    sd_totalphenolics = sd(totalphenolics),
    cv_dglucose = sd_dglucose / mean_dglucose * 100,
    cv_fructose = sd_fructose / mean_fructose * 100,
    cv_proteins = sd_proteins / mean_proteins * 100,
    cv_alkenylphenols = sd_alkenylphenols / mean_alkenylphenols * 100,
    cv_totalphenolics = sd_totalphenolics / mean_totalphenolics * 100,
    skew_dglucose = skewness(dglucose),
    skew_fructose = skewness(fructose),
    skew_proteins = skewness(proteins),
    skew_alkenylphenols = skewness(alkenylphenols),
    skew_totalphenolics = skewness(totalphenolics)
  ) %>%
  gather(key = "variable", value = "value", -tree, factor_key = TRUE)

### create a set plots for each variable by tree
head(summary)
summary$variable <- as.character(summary$variable)

summary <- summary %>%
  mutate(metric = case_when(
    startsWith(variable, "mean_") ~ "mean",
    startsWith(variable, "sd_") ~ "sd",
    startsWith(variable, "cv_") ~ "cv",
    startsWith(variable, "skew_") ~"skewness",
    TRUE ~ NA_character_
  ))

summary <- summary %>%
  mutate(molecule = case_when(
    endsWith(variable, "glucose") ~ "glucose",
    endsWith(variable, "proteins") ~ "proteins",
    endsWith(variable, "fructose") ~ "fructose",
    endsWith(variable, "alkenylphenols") ~ "alkenylphenols",
    endsWith(variable, "totalphenolics") ~"totalphenolics",
    TRUE ~ NA_character_
  ))

head(summary)
summary$molecule <- factor(summary$molecule, levels = c("glucose", "fructose", "proteins", "totalphenolics", "alkenylphenols"))

supp.labs <- c("Total\nalkenylphenols", "Glucose", "Fructose", "Total\nproteins", "Total\nphenolics")
names(supp.labs) <- c("alkenylphenols", "glucose","fructose","proteins", "totalphenolics")
labeller_func <- as_labeller(supp.labs)
summary$metric <- as.factor(summary$metric)

mean <- summary %>%
  filter(metric %in% c("mean"))
mean_1 <- mean %>%
  filter(molecule %in% c("alkenylphenols", "glucose", "fructose"))
meang_1 <- ggplot(mean_1, aes(x=molecule, y=value)) +
  geom_jitter(aes(x=molecule, y=value, color = tree), size = 4, width = 0.2) + 
  geom_boxplot(outlier.alpha = 0, fill = "lightgray",alpha = 0.3) +
  scale_color_viridis_d(option = "D", name = "Plant ID", labels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10")) + 
  scale_x_discrete(labels = c("Glucose", "Fructose", "Total\nalkenylphenols")) + 
  labs(x = " ", y = "Mean\nPercentages\n(mg/mg DW)") +
  theme_classic(base_size = 25) +
  theme(legend.position = "none",
        plot.margin = unit(c(5.5, 5.5, 0, 5.5), "pt"))  
meang_1

mean_2 <- mean %>%
  filter(molecule %in% c("totalphenolics", "proteins"))
meang_2 <- ggplot(mean_2, aes(x=molecule, y=value)) +
  geom_jitter(aes(x=molecule, y=value, color = tree), size = 4, width = 0.2) + 
  geom_boxplot(outlier.alpha = 0, fill = "lightgray",alpha = 0.3) +
  scale_color_viridis_d(option = "D", name = "Plant ID", labels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10")) + 
  scale_x_discrete(labels = c("Total\nproteins", "Total\nphenolics")) + 
  labs(x = " ", y = " ") +
  theme_classic(base_size = 25) +
  theme(legend.position = "none",
        plot.margin = unit(c(5.5, 5.5, 0, 5.5), "pt"))  
meang_2

mean_graphs <- ggarrange(meang_1,
                         meang_2,
                         ncol = 2,
                         nrow = 1)
mean_graphs +   theme(legend.position = "none",
                      plot.margin = unit(c(0, 5.5, 0, 5.5), "pt")) 

###### 

CV <- summary %>%
  filter(metric %in% c("cv"))

CVg <- ggplot(CV, aes(x=molecule, y=value)) +
  geom_jitter(aes(x=molecule, y=value, color = tree), size = 4, width = 0.2) + 
  geom_boxplot(outlier.alpha = 0, fill = "lightgray",alpha = 0.3) +
  scale_color_viridis_d(option = "D", name = "Plant ID", labels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10")) + 
  scale_x_discrete(labels = c("Glucose", "Fructose", "Total\nproteins", "Total\nphenolics", "Total\nalkenylphenols")) + 
  labs(x = " ", y = "CV") +
  theme_classic(base_size = 25) +
  theme(legend.position = "none",
        plot.margin = unit(c(0, 5.5, 0, 5.5), "pt"))  
CVg

###
SD <- summary %>%
  filter(metric %in% c("sd"))

SD_1 <- SD %>%
  filter(molecule %in% c("alkenylphenols", "glucose", "fructose"))
SDg_1 <- ggplot(SD_1, aes(x=molecule, y=value)) +
  geom_jitter(aes(x=molecule, y=value, color = tree), size = 4, width = 0.2) + 
  geom_boxplot(outlier.alpha = 0, fill = "lightgray",alpha = 0.3) +
  scale_color_viridis_d(option = "D", name = "Plant ID", labels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10")) + 
  scale_x_discrete(labels = c("Glucose", "Fructose", "Total\nalkenylphenols")) + 
  labs(x = " ", y = "SD\nPercentages\n(mg/mg DW)") +
  theme_classic(base_size = 25) +
  theme(legend.position = "none",
        plot.margin = unit(c(0, 5.5, 0, 5.5), "pt"))  
SDg_1

SD_2 <- SD %>%
  filter(molecule %in% c("totalphenolics", "proteins"))
SDg_2 <- ggplot(SD_2, aes(x=molecule, y=value)) +
  geom_jitter(aes(x=molecule, y=value, color = tree), size = 4, width = 0.2) + 
  geom_boxplot(outlier.alpha = 0, fill = "lightgray",alpha = 0.3) +
  scale_color_viridis_d(option = "D", name = "Plant ID", labels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10")) + 
  scale_x_discrete(labels = c("Total\nproteins", "Total\nphenolics")) + 
  labs(x = " ", y = " ") +
  theme_classic(base_size = 25) +
  theme(legend.position = "none",
        plot.margin = unit(c(0, 5.5, 0, 5.5), "pt"))  

SDg_2

SD_graphs <- ggarrange(SDg_1,
                       SDg_2,
                       ncol = 2,
                       nrow = 1)
SD_graphs +   theme(legend.position = "none",
                    plot.margin = unit(c(5.5, 5.5, 0, 5.5), "pt")) 

###

Skewness <- summary %>%
  filter(metric %in% c("skewness"))

skewnessg <- ggplot(Skewness, aes(x=molecule, y=value)) +
  geom_jitter(aes(x=molecule, y=value, color = tree), size = 4, width = 0.2) + 
  geom_boxplot(outlier.alpha = 0, fill = "lightgray",alpha = 0.3) +
  scale_color_viridis_d(option = "D", name = "Plant ID", labels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10")) + 
  scale_x_discrete(labels = c("Glucose", "Fructose", "Total\nproteins", "Total\nphenolics", "Total\nalkenylphenols")) + 
  labs(x = " ", y = "Skewness") +
  theme_classic(base_size = 25) +
  theme(legend.position = "bottom",
        plot.margin = unit(c(0, 0, 0, 0), "pt"))  

skewnessg

figure1 <- ggarrange(mean_graphs,
                        SD_graphs,
                        CVg,
                        skewnessg,
                        ncol = 1,
                        nrow = 4,
                        common.legend = TRUE,
                        legend="bottom",
                        align = "h")
ggsave(file="Figure_1.jpg", 
       plot= figure1,
       width=11,height = 18,units="in",dpi=400)

######################
### RIDGELINE PLOT ###
######################

supp.labs <- c("Total\nalkenylphenols", "Glucose", "Fructose", "Total\nproteins", "Total\nphenolics")
names(supp.labs) <- c("alkenylphenols", "dglucose","fructose","proteins", "totalphenolics")

data_pivoted <- data %>% pivot_longer(
  cols = -c(sampleID, fruitID, tree),
  names_to = "molecules",
  values_to="values") %>%
  filter(molecules %in% c("totalphenolics", "dglucose", "fructose", "alkenylphenols", "proteins", "alkenylphenols"))%>%
  filter(!tree %in% c("2", "3"))

data_pivoted$molecules <- factor(data_pivoted$molecules, levels = c("dglucose", "fructose", "proteins", "totalphenolics", "alkenylphenols"))

figure_2 <- ggplot(data_pivoted, aes(x = values, y = tree, fill = factor(tree))) +
  geom_density_ridges(alpha = 0.6, jittered_points = TRUE, point_alpha=1,point_shape=21, point_size = 1) +
  facet_grid(~molecules, labeller = labeller(molecules= supp.labs), scales = "free_x") +
  scale_fill_viridis_d() +
  theme_blank(base_size = 14) +
  xlab("Percentage (%) of nutrients and defensive metabolites\n(mg/mg DW)") +
  ylab("Plant ID") +
  labs(fill = "Tree") + 
  theme(legend.position = "none",
        panel.spacing = unit(0.7, "cm"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
ggsave(file="Figure_2.jpg", 
       plot= figure_2,
       width=7,height=5,units="in",dpi=400)

#############
### NMDS ###
#############
set.seed(3)
data <- data %>% filter(!tree %in% c("2", "3"))
data <- na.exclude(data)
peak_data <- data[, 4:16] # select just the peaks
zero_counts <- colSums(peak_data == 0)
z_data <- as.data.frame(lapply(peak_data, function(x) (x - mean(x)) / sd(x)))
matrix <- as.matrix(z_data) # turn data frame into matrix
nmds_results <- metaMDS(matrix, 
                        k = 2,
                        distance = "euclidean") # Number of iterations
plot(nmds_results, type = "t")
nmds_results$stress
data.scores <- as.data.frame(scores(nmds_results$points))
data.scores$tree <- data$tree
data.scores

centroids <- data.scores %>% 
  group_by(tree) %>% 
  summarize(MDS1 = mean(MDS1), MDS2 = mean(MDS2))

figure_3 <- ggplot(data.scores, aes(x = MDS1, y = MDS2, colour = tree)) +
  geom_point(size = 5, alpha = 0.8) +
  geom_point(data = centroids, aes(x = MDS1, y = MDS2, fill = tree), size = 5, shape = 22, color = "black") +
  theme_classic(base_size = 13) +
  scale_color_viridis(option = "D", discrete=TRUE, name = "Plant ID", labels = c("1", "4", "5", "6", "7", "8", "9", "10")) +
  scale_fill_viridis(option = "D", discrete=TRUE, name = "Plant ID", labels = c("1", "4", "5", "6", "7", "8", "9", "10")) +
  ylab ("NMDS2") +
  xlab ("NMDS1") +
  stat_ellipse() + 
  theme(legend.position = "bottom")
figure_3
ggsave(file="Figure_3.jpg", 
       plot=figure_3,
       width=8,height=8,units="in",dpi=300)

# PERMANOVA
adonis <- adonis2(matrix ~ data.scores$tree, distance = "euclidean", perm=999)
adonis
write.csv(adonis, file = "permanova_results.csv")
# DISPERSION - betadisper
dist_matrix <- vegdist(matrix, method = "bray")
groups <- data$tree
dispersal <- betadisper(dist_matrix, groups)
anova(dispersal)
tukey <- TukeyHSD(dispersal)
write.csv(tukey$group, file = "betadisper_results.csv")

####################################
### Variance Components Analysis ###
####################################

### Gaussian distribution

glmm_glucose_gaussian <- glmmTMB(dglucose ~ (1|tree/fruitID), data = data) ### do not understand the residuals
summary(glmm_glucose_gaussian)
shapiro.test(resid(glmm_glucose_gaussian)) # normal 
diagnose(glmm_glucose_gaussian) # Unusually large Z-statistics (|x|>5)

data$fructose_gaussian <- data$fructose + 0.001 ### the dataset has a 0 in fructose, so a GLMM with a gamma distribution would not run
glmm_fructose_gaussian <- glmmTMB(fructose ~ (1|tree/fruitID), data = data)
summary(glmm_fructose_gaussian) 
shapiro.test(resid(glmm_fructose_gaussian)) # not normal
diagnose(glmm_fructose_gaussian) # Unusually large Z-statistics (|x|>5)

glmm_proteins_gaussian <- glmmTMB(proteins ~ (1|tree/fruitID), data = data)
summary(glmm_proteins_gaussian)
shapiro.test(resid(glmm_proteins_gaussian)) # normal
diagnose(glmm_proteins_gaussian) # Unusually large Z-statistics (|x|>5)

glmm_alkenylphenols_gaussian <- glmmTMB(alkenylphenols ~ (1|tree/fruitID), data = data)
summary(glmm_alkenylphenols_gaussian)
shapiro.test(resid(glmm_alkenylphenols_gaussian)) # not normal
diagnose(glmm_alkenylphenols_gaussian) # Unusually large Z-statistics (|x|>5)

glmm_totalphenolics_gaussian <- glmmTMB(totalphenolics ~ (1|tree/fruitID), data = data)
summary(glmm_totalphenolics_gaussian)
shapiro.test(resid(glmm_totalphenolics_gaussian)) # normal
diagnose(glmm_totalphenolics_gaussian) # Unusually large Z-statistics (|x|>5)

### Gamma distribution 

glmm_glucose_gamma <- glmmTMB(dglucose ~ (1|tree/fruitID), data = data, family=Gamma(link = log))
summary(glmm_glucose_gamma)
shapiro.test(resid(glmm_glucose_gamma)) # normal
diagnose(glmm_glucose_gamma) # Unusually large coefficients (|x|>10)

data$fructose_gamma <- data$fructose + 0.001 ### the dataset has a 0 in fructose, so a GLMM with a gamma distribution would not run
glmm_fructose_gamma <- glmmTMB(fructose_gamma ~ (1|tree/fruitID), data = data, family=Gamma(link = log))
summary(glmm_fructose_gamma)
shapiro.test(resid(glmm_fructose_gamma)) # not normal
diagnose(glmm_fructose_gamma) # Unusually large coefficients (|x|>10)

glmm_proteins_gamma <- glmmTMB(proteins ~ (1|tree/fruitID), data = data, family=Gamma(link = log))
summary(glmm_proteins_gamma)
shapiro.test(resid(glmm_proteins_gamma)) # not normal 
diagnose(glmm_proteins_gamma) # Unusually large coefficients (|x|>10)

glmm_alkenylphenols_gamma <- glmmTMB(alkenylphenols ~ (1|tree/fruitID), data = data, family=Gamma(link = log))
summary(glmm_alkenylphenols_gamma)
shapiro.test(resid(glmm_alkenylphenols_gamma)) # not normal 
diagnose(glmm_alkenylphenols_gamma) # Unusually large Z-statistics (|x|>5)

glmm_totalphenolics_gamma <- glmmTMB(totalphenolics ~ (1|tree/fruitID), data = data, family=Gamma(link = log))
summary(glmm_totalphenolics_gamma)
shapiro.test(resid(glmm_totalphenolics_gamma)) # normal 
diagnose(glmm_totalphenolics_gamma) # Unusually large coefficients (|x|>10)

### Beta distribution 
data$glucose_beta <- (data$dglucose)/100
data$fructose_beta <- (data$fructose)/100 + 0.0001
data$proteins_beta <- (data$proteins)/100
data$totalphenolics_beta <- (data$totalphenolics)/100
data$alkenylphenols_beta <- (data$alkenylphenols)/100

glmm_glucose_beta <- glmmTMB(glucose_beta ~ (1|tree/fruitID), data = data, family = beta_family(link="logit"))
summary(glmm_glucose_beta)
shapiro.test(resid(glmm_glucose_beta)) # normal 
diagnose(glmm_glucose_beta) # Unusually large Z-statistics (|x|>5)

glmm_fructose_beta <- glmmTMB(fructose_beta ~ (1|tree/fruitID), data = data, family = beta_family(link="logit"))
summary(glmm_fructose_beta) 
shapiro.test(resid(glmm_fructose_beta)) # not normal 
diagnose(glmm_fructose_beta) # Unusually large Z-statistics (|x|>5)

glmm_proteins_beta <- glmmTMB(proteins_beta ~ (1|tree/fruitID), data = data, family = beta_family(link="logit"))
summary(glmm_proteins_beta)
shapiro.test(resid(glmm_proteins_beta)) # normal 
diagnose(glmm_proteins_beta) # Unusually large coefficients (|x|>10)

glmm_alkenylphenols_beta <- glmmTMB(alkenylphenols_beta ~ (1|tree/fruitID), data = data, family = beta_family(link="logit"))
summary(glmm_alkenylphenols_beta)
shapiro.test(resid(glmm_alkenylphenols_beta)) # normal 
diagnose(glmm_alkenylphenols_beta) # Unusually large coefficients (|x|>10)

glmm_totalphenolics_beta <- glmmTMB(totalphenolics_beta ~ (1|tree/fruitID), data = data, family = beta_family(link="logit"))
summary(glmm_totalphenolics_beta)
shapiro.test(resid(glmm_totalphenolics_beta)) # normal 
diagnose(glmm_totalphenolics_beta) # Unusually large Z-statistics (|x|>5)

### Model comparison between Gaussian, gamma, and beta models using AICc

AIC(glmm_glucose_gaussian, glmm_glucose_gamma, glmm_glucose_beta)
AIC(glmm_fructose_gaussian, glmm_fructose_gamma, glmm_fructose_beta)
AIC(glmm_proteins_gaussian, glmm_proteins_gamma, glmm_proteins_beta)
AIC(glmm_alkenylphenols_gaussian, glmm_alkenylphenols_gamma, glmm_alkenylphenols_beta)
AIC(glmm_totalphenolics_gaussian, glmm_totalphenolics_gamma, glmm_totalphenolics_beta)

### beta models have the smallest AICc
vpa <- data.frame(
  molecule = c("Glucose", "Fructose", "Total proteins", "Total phenolics", "Total alkenylphenols",
               "Glucose", "Fructose", "Total proteins", "Total phenolics", "Total alkenylphenols"),
  type = c("Nutrient", "Nutrient", "Nutrient", "Defense", "Defense",
           "Nutrient", "Nutrient", "Nutrient", "Defense", "Defense"),
  values = c(51.946, 82.325, 50.025, 65.894, 50.541, 48.054, 17.675, 49.975, 34.106, 49.459),
  variance = c("Intraindividual", "Intraindividual", "Intraindividual", "Intraindividual", "Intraindividual",
               "Intraspecific", "Intraspecific", "Intraspecific", "Intraspecific", "Intraspecific")
)

head(vpa)
vpa$molecule <- factor(vpa$molecule, c("Fructose",
                                       "Total phenolics",
                                       "Glucose",
                                       "Total alkenylphenols",
                                       "Total proteins"))

figure_4 <- ggplot(vpa, aes(x = values, fill =  variance, y = molecule)) +
  theme_classic(base_size = 15) +
  geom_bar(stat = "identity", color = "NA", alpha = 0.9) + 
  scale_fill_viridis_d(name = "Variation") + 
  labs(x = "Variance %", y = "") + 
  theme(legend.position = "top")
ggsave(file="Figure_4.jpg", 
       plot= figure_4,
       width=6,height=8,units="in",dpi=400)

########################################
########################################
########################################
########################################

