#setwd("...")

library(tidyverse)
library(ggplot2)
library(ggpubr)
library(rstatix)
library(magrittr)
library(lsr)


AR_ICR <- read.csv(file = 'AR_weightByICR.csv',header=TRUE)
AR_Stim <- read.csv(file = 'AR_weightBypattern.csv',header=TRUE)

head(AR_ICR)
head(AR_Stim)

AR_ICR %>%
  reorder_levels("Stimulus", order = c("two", "four", "mono")) %>%
  group_by(Stimulus)
  
AR_Stim %>%
  reorder_levels("Stimulus", order = c("grating", "noise", "simple","complex")) %>%
  group_by(Stimulus)

aov_Stim <- anova_test(data = AR_Stim, dv = Weight, wid = Subj, within = Stimulus)
get_anova_table(aov_Stim)

aov_ICR <- anova_test(data = AR_ICR, dv = Weight, wid = Subj, within = Stimulus)
get_anova_table(aov_ICR)

# do permutation ANOVA to cross validate
model<-aovp(Weight ~ Stimulus,data=AR_Stim)
model<-aovp(Weight ~ Stimulus + Error(1/Subj),data=AR_Stim)

print(summary(model))

## do pairwise t-test
t_test_stim <- AR_Stim %>%
  pairwise_t_test(
    Weight ~ Stimulus, paired = TRUE,
    p.adjust.method = "bonferroni"
  )
t_test_stim

t_test_ICR <- AR_ICR %>%
  pairwise_t_test(
    Weight ~ Stimulus, paired = TRUE,
    p.adjust.method = "bonferroni"
  )
t_test_ICR

#effect size calculation
#for stimulus type
cohensD(AR_Stim[AR_Stim$Stimulus=='grating','Weight'], AR_Stim[AR_Stim$Stimulus=='noise','Weight'])
cohensD(AR_Stim[AR_Stim$Stimulus=='grating','Weight'], AR_Stim[AR_Stim$Stimulus=='simple','Weight'])
cohensD(AR_Stim[AR_Stim$Stimulus=='grating','Weight'], AR_Stim[AR_Stim$Stimulus=='complex','Weight'])
cohensD(AR_Stim[AR_Stim$Stimulus=='noise','Weight'], AR_Stim[AR_Stim$Stimulus=='simple','Weight'])
cohensD(AR_Stim[AR_Stim$Stimulus=='noise','Weight'], AR_Stim[AR_Stim$Stimulus=='complex','Weight'])
cohensD(AR_Stim[AR_Stim$Stimulus=='simple','Weight'], AR_Stim[AR_Stim$Stimulus=='complex','Weight'])
#for ICR levels
cohensD(AR_ICR[AR_ICR$Stimulus=='two','Weight'], AR_ICR[AR_ICR$Stimulus=='four','Weight'])
cohensD(AR_ICR[AR_ICR$Stimulus=='two','Weight'], AR_ICR[AR_ICR$Stimulus=='mono','Weight'])
cohensD(AR_ICR[AR_ICR$Stimulus=='four','Weight'], AR_ICR[AR_ICR$Stimulus=='mono','Weight'])



