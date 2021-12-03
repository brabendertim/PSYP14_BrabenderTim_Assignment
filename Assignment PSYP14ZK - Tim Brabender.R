###############################################################################
#
#
#             Assignment PSYP14 - Zoltan Kekecs - Tim Brabender
# 
#
###############################################################################
#
#
#
#
#                             Part 1 of 3
#
#
#
###############################################################################



# producing the whole library


library(boot) 
library(cAIC4)
library(car)
library(dplyr)
library(gridExtra)	
library(influence.ME) 	
library(imputeTS)
library(lattice) 
library(lmboot)
library(lme4)	
library(lmerTest) 
library(lmtest)
library(lm.beta) 
library(MuMIn)
library(psych)
library(r2glmm) 
library(sandwich)
library(tidyverse) 


# coefficient table function
coef_table = function(model) {
  require(lm.beta)
  mod_sum = summary(model)
  mod_sum_p_values = as.character(round(mod_sum$coefficients[,
                                                             4], 3))
  mod_sum_p_values[mod_sum_p_values != "0" & mod_sum_p_values !=
                     "1"] = substr(mod_sum_p_values[mod_sum_p_values != "0" &
                                                      mod_sum_p_values != "1"],
                                   2, nchar(mod_sum_p_values[mod_sum_p_values !=
                                              "0" & mod_sum_p_values != "1"]))
  mod_sum_p_values[mod_sum_p_values == "0"] = "<.001"
  mod_sum_table = cbind(as.data.frame(round(cbind(coef(model),
                    confint(model), c(0,lm.beta(model)$standardized.coefficients
                      [c(2:length(model$coefficients))])),2)), mod_sum_p_values)
  names(mod_sum_table) = c("b", "95%CI lb", "95%CI ub", "Std.Beta",
                           "p-value")
  mod_sum_table["(Intercept)", "Std.Beta"] = "0"
  return(mod_sum_table)
}



# standardized beta coefficients function
stdCoef.merMod <- function(object) {	
  sdy <- sd(getME(object,"y"))	
  sdx <- apply(getME(object,"X"), 2, sd)	
  sc <- fixef(object)*sdx/sdy	
  se.fixef <- coef(summary(object))[,"Std. Error"]	
  se <- se.fixef*sdx/sdy	
  return(data.frame(stdcoef=sc, stdse=se))	
}	

# loading the dataset

data_sample_1 <- read.csv("https://tinyurl.com/ha-dataset1")

# Check the varibles for coding errors

summary(data_sample_1)

# pain Max is 55 with a scale from 1 to 10
# STAI Min is 4.2 with a  scale from 20 to 80
# IQ levels are unusually ranged, maybe run a histogram


# explore pain scale

hist_pain = data_sample_1 %>% 	
  ggplot() +	
  aes(x = pain) +	
  geom_histogram()

hist_pain

# we should delete the outlier
# find outlier

data_sample_1 %>% 
  filter(pain > 10)

# ID 88 should be deleted

# STAI trait

hist_STAI = data_sample_1 %>% 	
  ggplot() +	
  aes(x = STAI_trait) +	
  geom_histogram()

hist_STAI

# find the outlier

data_sample_1 %>% 
  filter(STAI_trait < 20)

# ID 34 will be dropped

# IQ levels 

hist_IQ = data_sample_1 %>% 	
  ggplot() +	
  aes(x = IQ) +	
  geom_histogram()

hist_IQ

# IQ seems to be about normally distributed with a small unusual peak at ~ 140
# no adjustments necessary here

# check other variables for unusual data


# sex and age
density_sex_age = data_sample_1 %>% 	
  ggplot() +	
  aes(x = age, fill = sex) +	
  geom_density(alpha = 0.3)

density_sex_age



# pain catastrophizing
hist_pain_cat = data_sample_1 %>% 	
  ggplot() +	
  aes(x = pain_cat) +	
  geom_histogram()

hist_pain_cat


# cortisol_serum
hist_serum = data_sample_1 %>% 	
  ggplot() +	
  aes(x = cortisol_serum) +	
  geom_histogram()

hist_serum


# cortisol_saliva
hist_saliva = data_sample_1 %>% 	
  ggplot() +	
  aes(x = cortisol_saliva) +	
  geom_histogram()

hist_saliva


# mindfulness
hist_mindfulness = data_sample_1 %>% 	
  ggplot() +	
  aes(x = mindfulness) +	
  geom_histogram()

hist_mindfulness

# weight
hist_weight = data_sample_1 %>% 	
  ggplot() +	
  aes(x = weight) +	
  geom_histogram()

hist_weight

# income
hist_income = data_sample_1 %>% 
  ggplot() +
    aes(x = household_income) +
      geom_histogram()
  
  
hist_income


# run descriptives
describe(data_sample_1)

# kurtosis on pain variable is excessively high
# STAI trait has a high kurtosis as well


# delete impossible scores ID 88 and ID 34
data_sample_1_cleaned <- data_sample_1[-c(34,88),]
  
# mutate sex to a factor
data_sample_1_cleaned = data_sample_1_cleaned %>% 	
  mutate(sex = factor(sex))	

levels(data_sample_1_cleaned$sex)

# Create first model using age and sex as predictors of pain
mod_1 <- lm(pain ~ age + sex, data = data_sample_1_cleaned)

mod_1

# visualize mod 1 relationships
scatter_pain_age = data_sample_1_cleaned %>% 	
  ggplot() +	
  aes(x = age, y = pain) +	
  geom_point()+	
  geom_smooth(method = "lm")	

scatter_pain_sex = data_sample_1_cleaned %>% 	
  ggplot() +	
  aes(x = sex, y = pain) +	
  geom_point()+	
  geom_smooth(method = "lm")	

grid.arrange(scatter_pain_age, scatter_pain_sex, nrow = 1)	



# Create second model using age, sex, STAI, pain catastrophizing, mindfulness
# and cortisol measures to predict pain
mod_2 <- lm(pain ~ age + sex + STAI_trait + pain_cat + mindfulness + 
              cortisol_serum + cortisol_saliva, data = data_sample_1_cleaned)

mod_2

# Check the variables included in model 2 for coding errors as well as the model
# itself for influential outliers

# cook's distance test
mod_2 %>% 
  plot(which = 4)

mod_2 %>% 
  plot(which = 5)


# threshold for Cook's distance 4/N (4/158)
4/158

# threshold is at about 0.025
# multiple cases are above the threshold, with 47, 74 and 86 sticking out
data_sample_1_cleaned %>% 
  slice(c(47,74,86))

# scores look normal


# Check final model to see if assumptions of linear regression hold true

# Normality of the residuals
mod_2 %>% 
  plot(which = 2)

# 104 and 85 stick out


# histogram of residuals
residuals_mod_2 = enframe(residuals(mod_2))
residuals_mod_2 %>% 
  ggplot()+
    aes(x=value)+
      geom_histogram()

# skew and kurtosis
describe(residuals(mod_2))

# skew and kurtosis is between -1 and 1 so no violation there


# Linearity of the relationship
# check residual plot
mod_2 %>% 
  residualPlots()

# tests all show no significance


# Homogeneity of variance / homoscedasticity
mod_2 %>% 
  plot(which = 3)

mod_2 %>% 
  ncvTest()

# test is not significant
mod_2 %>% 
  bptest()

# bptest is also not significant

# Multicollinearity
mod_2 %>% 
  vif()

# vif tests for cortisol_serum and cortisol_saliva return scores
# above 3. Adjustments will be made

# make appropriate changes and rerun the assumption checks:
# Multicollinearity needs to be accounted for. theres is no structural
# multicollinearity as there are no higher order terms nor interaction terms
# Data multicollinearity seems more likely

# create correlation matrix
data_sample_1_cleaned %>% 
  select(pain, cortisol_serum, cortisol_saliva) %>% 
    pairs.panels(col= "red", lm=T)


# summarize model 2
summary(mod_2)


# one could linearily combine both saliva levels by taking the means
# of both and putting it in a new variable

# create new variable "cortisol"
data_sample_1_cleaned_cortisol_variable <- data_sample_1_cleaned %>% 
  mutate(cortisol = (cortisol_saliva + cortisol_serum)/2)



# build new model with new predictor cortisol
mod_3 <- lm(pain ~ age + sex + STAI_trait + pain_cat + mindfulness + cortisol,
            data = data_sample_1_cleaned_cortisol_variable)

mod_3


# recheck assumptions of the new model


# cook's distance test
mod_3 %>% 
  plot(which = 4)

mod_3 %>% 
  plot(which = 5)


# threshold for Cook's distance 4/N (4/158)
4/158

# threshold is at about 0.025
# multiple cases are above the threshold, with 47, 74 and 86 sticking out
data_sample_1_cleaned_cortisol_variable %>% 
  slice(c(47,74,86))

# scores look normal



# Check final model to see if assumptions of linear regression hold true

# Normality of the residuals
mod_3 %>% 
  plot(which = 2)

# 104 and 85 stick out


# histogram of residuals
residuals_mod_3 = enframe(residuals(mod_3))
residuals_mod_3 %>% 
  ggplot()+
  aes(x=value)+
  geom_histogram()

# skew and kurtosis
describe(residuals(mod_3))

# skew and kurtosis is between -1 and 1 so no violation there


# Linearity of the relationship
# check residual plot
mod_3 %>% 
  residualPlots()

# tests all show no significance


# Homogeneity of variance / homoscedasticity
mod_3 %>% 
  plot(which = 3)

mod_3 %>% 
  ncvTest()

# test is not significant
mod_3 %>% 
  bptest()

# bptest is also not significant

# Multicollinearity
mod_3 %>% 
  vif()


# vif tests all show values below 3 so multicollinearity is not present

# report model test statistics (R^2, F, df, p value)
summary(mod_3)
summary(mod_1)

# Adjusted r-squared: 0.52, F statistic: 29.25, DF: 6/151, p-value < 0.001
coef_table(mod_3)

# report statistics describing the coefficients of the predictors in a table
# format (unstandardized regression coefficients and 95% CIs, standardized
# regression coefficients (B and Beta values), and p values)
summary(mod_3)
confint(mod_3)
lm.beta(mod_3)

# Write up the regression equation of model 2
# pain = 0.40 + (-0.03)*age + 0.16*sexmale + (-0.02)*STAI_trait + 0.13*pain_cat+ (-0.26)*mindfulness + 0.65*cortisol
# x1 .. xn stands for the different predictors


# Compare the two models in terms of how much variance they can explain of 
# pain variability (probably adjusted r^2)
summary(mod_1)$adj.r.squared
summary(mod_3)$adj.r.squared

# model 1 can explain 7.36% of the variance while model 3 can explain 51.92%



# report AIC for both models and the F test statistic and p value of the model
# comparison returned by the anova function
AIC(mod_1)
AIC(mod_3)

# mod 1 scores 574.13 in the AIC while mod 3 scores 474.37 and is therefore 
# showing better fit
anova(mod_1, mod_3)

# F test returns 36.91 and a p value of <0.001 with 4 DF




###############################################################################
#
#
#                   Assignment PSYP14 - Tim Brabender
# 
#
###############################################################################
#
#
#
#
#                             Part 2 of 3
#
#
#
###############################################################################
#
#
#
# create a new model using age, sex, STAI, pain catastrophing, mindfulness,
# serum sortisol, weight, IQ and household income as predictors for pain


initial_model <- lm(pain ~ age + sex + STAI_trait + pain_cat + mindfulness +
                      cortisol_serum + weight +  IQ + household_income, 
                    data = data_sample_1_cleaned_cortisol_variable)


# run model diagnostics for initial model


# cook's distance test
initial_model %>% 
  plot(which = 4)

initial_model %>% 
  plot(which = 5)


# threshold for Cook's distance 4/N (4/158)
4/158

# threshold is at about 0.025
# multiple cases are above the threshold, with 47, 85 and 86 sticking out
data_sample_1_cleaned_cortisol_variable %>% 
  slice(c(47,85,86))

# scores look normal

# Check final model to see if assumptions of linear regression hold true

# Normality of the residuals
initial_model %>% 
  plot(which = 2)

# 104 and 85 stick out


# histogram of residuals
residuals_initial_model = enframe(residuals(initial_model))
residuals_initial_model %>% 
  ggplot()+
  aes(x=value)+
  geom_histogram()

# skew and kurtosis
describe(residuals(initial_model))

# skew and kurtosis is between -1 and 1 so no violation there


# Linearity of the relationship
# check residual plot
initial_model %>% 
  residualPlots()


# tests all show no significance


# Homogeneity of variance / homoscedasticity
initial_model %>% 
  plot(which = 3)

initial_model %>% 
  ncvTest()

# test is not significant
initial_model %>% 
  bptest()

# bptest is also not significant

# Multicollinearity
initial_model %>% 
  vif()


# All vif tests return scores below 3.No adjustments necessary.


# run a backwards regression using age, sex, STAI, pain catastrophing, mindfulness,
# serum sortisol, weight, IQ and household income
backward_model <- step(initial_model, direction = "backward") 


# run full model diagnostics
backward_model %>% 
  plot(which = 4)

backward_model %>% 
  plot(which = 5)

backward_model %>% 
  plot(which = 2)

residuals_backward_model = enframe(residuals(backward_model))
residuals_backward_model %>% 
  ggplot()+
  aes(x=value)+
  geom_histogram()

describe(residuals(backward_model))

backward_model %>% 
  residualPlots()

backward_model %>% 
  plot(which = 3)

backward_model %>% 
  ncvTest()

backward_model %>% 
  bptest()

backward_model %>% 
  vif()


# everything looks good except a big outlier (No. 47)

# run full regression model of the last exercise again
theory_based_model <- lm(pain ~ age + sex + STAI_trait + pain_cat + mindfulness
                         + cortisol, data =
                           data_sample_1_cleaned_cortisol_variable)



# compare both model's AIC scores
AIC(initial_model)
AIC(theory_based_model)
AIC(backward_model)


# run anova
anova(theory_based_model, backward_model)


# both tests show that the new backward model is not significantly better

# test both models on a new dataset
data_sample_2 <- read.csv("https://tinyurl.com/87v6emky")


# create the cortisol variable in dataset 2
data_sample_2_cortisol_variable <- data_sample_2 %>% 
  mutate(cortisol = (cortisol_saliva + cortisol_serum)/2)



# calculate predicted values 	
pred_test_theory <- predict(theory_based_model, data_sample_2_cortisol_variable)	
pred_test_backward <- predict(backward_model, data_sample_2_cortisol_variable)	

# now we calculate the sum of squared residuals 	
RSS_test_theory = sum((data_sample_2_cortisol_variable[,"pain"] - pred_test_theory)^2)	
RSS_test_back = sum((data_sample_2_cortisol_variable[,"pain"] - pred_test_backward)^2)	
RSS_test_theory	
RSS_test_back	

# theory model performed better with less Residual sum of squares


# write up regression equation of the backward model
summary(backward_model)

coef_table(backward_model)

# pain = 1.28 + (-0.04)*age + 0.11*pain_cat + (-0.27)*mindfulness + 
# 0.53*cortisol_serum





###############################################################################
#
#
#                   Assignment PSYP14 - Tim Brabender
# 
#
###############################################################################
#
#
#
#
#                             Part 3 of 3
#
#
#
###############################################################################
#
#
#
# import dataset 3

data_file_3 <- read.csv("https://tinyurl.com/b385chpu")


# check data file for coding errors
summary(data_file_3)
describe(data_file_3)

# there is a negative household income data point which needs to be changed
# there is a coding error for sex 


# check the income variable for the negative number
data_file_3 %>% 
  filter(household_income < 0)


# the negative household income data point is recoded to NA
data_file_3_na <- data_file_3 %>% 
  mutate(household_income = replace(household_income, which(household_income<0),
                                    NA))

head(data_file_3_na)

# NAs replaced by the variable mean
data_file_3_cleaned <- na_mean(data_file_3_na)

head(data_file_3_cleaned)


# recode the coding error concerning sex
data_file_3_cleaned$sex[25]="female"


# mutate sex to a factor
data_file_3_cleaned = data_file_3_cleaned %>% 	
  mutate(sex = factor(sex))


levels(data_file_3_cleaned$sex)

# mutate the hospital ID to a factor
data_file_3_cleaned = data_file_3_cleaned %>% 
  mutate(hospital = factor(hospital))

levels(data_file_3_cleaned$hospital)

# create cortisol variable
data_file_3_cleaned_cortisol_variable <- data_file_3_cleaned %>% 
  mutate(cortisol = (cortisol_saliva + cortisol_serum)/2)


# explorative plotting
hospital_plot_1 = data_file_3_cleaned_cortisol_variable %>%
  ggplot() + aes(y = pain, x = cortisol, color = hospital) +
  geom_point(size = 4) + geom_smooth(method = "lm", se = F,
                                     fullrange = TRUE)+ 
  facet_wrap(~hospital, ncol = 5)
hospital_plot_1


# build linear mixed model on data file 3 accounting for the clustering of the
# data at different hospitals - only fit a random intercept model (random 
# intercept is hospital-ID) and same fixed effect predictors as in part 1
mod_mixed_1 <- lmer(pain ~ age + sex + STAI_trait + pain_cat + mindfulness + 
                      cortisol + (1|hospital), 
                    data = data_file_3_cleaned_cortisol_variable)

################################################################################
# model diagnostics

# cook's distance test

mod_mixed_1 %>% 
  plot(which = 4)

mod_mixed_1 %>% 
  plot(which = 5)

# no funneling

# run influence test 
influence_group <- influence(mod_mixed_1, group="hospital")$alt.fixed

influence_data <- as_tibble(influence_group) %>% 
  gather(colnames(influence_group), value=coefficient, key=predictor)

plot_influence_mixed_model <- influence_data %>% 
  ggplot() + aes(x = 1, y = coefficient, group = predictor) +
  geom_violin() + geom_jitter(width = 0.2) + facet_wrap(~predictor,
                                                        scales = "free")
plot_influence_mixed_model

# there is an influential outlier in the sexmale variable

# Normality of the residuals
qqmath(mod_mixed_1, id = 0.05)

data_file_3_cleaned_cortisol_variable = data_file_3_cleaned_cortisol_variable %>%
  mutate(resid = residuals(mod_mixed_1))


# normality of random effects
qqmath(ranef(mod_mixed_1))

#looks good
random_effects <- as.data.frame(ranef(mod_mixed_1)[[1]])
names(random_effects) = c("intercept")

random_effects %>%
  ggplot() + aes(sample = intercept) + stat_qq() + stat_qq_line()

random_effects %>%
  ggplot() + aes(x = intercept) + geom_histogram()

# concerning histogram
describe(random_effects$intercept)$skew


describe(random_effects$intercept)$kurtosis

# there is strong reason to believe that something is wrong here
# the normality of random effects is not given


# Linearity of the relationship

plot(mod_mixed_1, arg = "pearson")

summary(mod_mixed_1)

# looking at fixed effect predictors and residuals
data_file_3_cleaned_cortisol_variable %>% 
  ggplot()+
  aes(x=cortisol, y=resid)+
  geom_point()

data_file_3_cleaned_cortisol_variable %>% 
  ggplot()+
  aes(x=sex, y=resid)+
  geom_point()

data_file_3_cleaned_cortisol_variable %>% 
  ggplot()+
  aes(x=STAI_trait, y=resid)+
  geom_point()

data_file_3_cleaned_cortisol_variable %>% 
  ggplot()+
  aes(x=pain_cat, y=resid)+
  geom_point()

data_file_3_cleaned_cortisol_variable %>% 
  ggplot()+
  aes(x=mindfulness, y=resid)+
  geom_point()

#nothing wrong here

# Homoscedasticity
homoscedasticity_model = lm(resid^2 ~ hospital, 
                            data = data_file_3_cleaned_cortisol_variable)
summary(homoscedasticity_model)

# p value looks alright (>.05)

# cyclone plot
IQR_hospitals = sapply(split(data_file_3_cleaned_cortisol_variable,
                             f = data_file_3_cleaned_cortisol_variable$hospital),
                       function(x) IQR(x$resid))

rank = rank(IQR_hospitals)

data_file_3_cleaned_cortisol_variable$rank = rep(rank, each =
                                                   length(c("hospital_1", "hospital_10","hospital_2", "hospital_3", 
                                                            "hospital_4", "hospital_5", "hospital_6","hospital_7", 
                                                            "hospital_8", "hospital_9")))



hospital_plot = unique(data_file_3_cleaned_cortisol_variable$hospital[order
                                                                      (data_file_3_cleaned_cortisol_variable$rank)])


ggplot(data_file_3_cleaned_cortisol_variable, aes(y=resid, x=factor(rank), 
                                                  labels=hospital))+
  geom_boxplot()+ scale_x_discrete(labels = hospital_plot) + coord_flip()

# heteroscedasticity is to be suspected

# Multicollinearity

pairs.panels(data_file_3_cleaned_cortisol_variable[, c("age", "sex", "STAI_trait",
                                                       "pain_cat", "mindfulness", "cortisol")], col = "red", lm = T)

# there appears to be some correlation between STAI and pain_cat, pain_cat and 
# mindfulness, and STAI and cortisol. 


# conclusion:
# influence test shows influential outlier on sexmale variable
# Normality looks good
# Normality of random effects however the histogram and kurtosis give reason 
# for concern
# Linearity looks good
# Heteroscedasticity is to be suspected


################################################################################


# note model coefficients and confidence intervals for the coefficients for all
# fixed effect predictors
summary(mod_mixed_1)
confint(mod_mixed_1)
stdCoef.merMod(mod_mixed_1)

# compute variance explained by the fixed effect predictors using marginal R^2
# and varinace explained by the fixed and random effect terms combined using
# conditional R^2
r.squaredGLMM(mod_mixed_1)	


# import data file 4 and create cortisol variable
data_file_4 <- read.csv("https://tinyurl.com/4f8thztv")

data_file_4_cortisol_variable <- data_file_4%>% 
  mutate(cortisol = (cortisol_saliva + cortisol_serum)/2)

# look at data file 4 for coding errors
summary(data_file_4_cortisol_variable)
describe(data_file_4_cortisol_variable)

# make hospital ID and sex a factor
data_file_4_cortisol_variable <- data_file_4_cortisol_variable %>% 
  mutate(sex = factor(sex))

levels(data_file_4_cortisol_variable$sex)

data_file_4_cortisol_variable <- data_file_4_cortisol_variable %>% 
  mutate(hospital = factor(hospital))

levels(data_file_4_cortisol_variable$hospital)

# nothing wrong here


# Include prediction of mixed model in data file 4
data_file_4_cortisol_variable = data_file_4_cortisol_variable %>% 		
  mutate(prediction_mixed_model = predict(mod_mixed_1))

# include mean variable
mod_mean <- lm(pain ~ 1, data = data_file_4_cortisol_variable)

data_file_4_cortisol_variable = data_file_4_cortisol_variable %>% 		
  mutate(prediction_mean = predict(mod_mean))


# check predictions for errors
summary(predict(mod_mixed_1))
describe(predict(mod_mixed_1))



# compute RSS for the mixed model
RSS = sum((data_file_4_cortisol_variable[, "pain"] - 
             data_file_4_cortisol_variable$prediction_mixed_model)^2)	
RSS




TSS = sum((data_file_4_cortisol_variable[, "pain"] - 
             data_file_4_cortisol_variable$prediction_mean)^2)
TSS

1-(RSS/TSS)

# Build a new linear mixed effects model on dataset 3 predicting pain but only
# use the most influential predictor from the previous model. Allow for both
# random intercept and random slope

summary(mod_mixed_1)

backward_mixed <- step(mod_mixed_1, direction = "backward")

summary(backward_mixed)

# model summary shows tht cortisol variable is the most influential predictor

mod_mixed_slope <- lmer(pain~(cortisol | hospital),
                        data_file_3_cleaned_cortisol_variable)


data_file_3_cleaned_cortisol_variable = data_file_3_cleaned_cortisol_variable %>% 		
  mutate(prediction_slope = predict(mod_mixed_slope))


# visualize fitted regression lines for each hospital
data_file_3_cleaned_cortisol_variable %>%
  ggplot() + aes(y = pain, x = cortisol, group = hospital) +
  geom_point(aes(color = hospital), size = 4) + 
  geom_line(color = "red",aes(y = prediction_slope, x = cortisol)) + 
  facet_wrap(~hospital, ncol = 5)




################################################################################
#
#
#
#                           Tim Brabender - 03.12.2021
#
#
#
################################################################################













