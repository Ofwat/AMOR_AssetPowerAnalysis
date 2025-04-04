################################################################################
#                   Author: Joshua Thompson
#   O__  ----       Email:  joshua.thompson@ofwat.gov.uk
#  c/ /'_ ---
# (*) \(*) --
# ======================== Script  Information =================================
# PURPOSE: Power analysis for asset health sampling 
#
# PROJECT INFORMATION:
#   Name: Power analysis for asset health sampling
#
# HISTORY:----
#   Date		        Remarks
#	-----------	   ---------------------------------------------------------------
#	 30/05/2024    Created script                                   JThompson (JT)
#===============================  Environment Setup  ===========================
#==========================================================================================

# load libraries

library(tidyverse)
library(tidyr)
library(pwr)
library(pwr2)
library(extrafont)

loadfonts()

# set 
output_dir <- "directory"

# asset numbers (hard coded)
asset_data <- tribble(
  ~Asset_Group, ~ANH, ~WSH, ~HDD, ~NES, ~SWB, ~SVE, ~SRN, ~TMS, ~NWT, ~WSX, ~YKY, ~AFW, ~BRL, ~PRT, ~SEW, ~SSC, ~SES,
  "Service reservoirs",       252, 429, 84, 304, 327, 482, 224, 240, 347, 299, 396, 108, 107, 17, 169, 57, 31,
  "Contact tanks",            128, 65, 5, 50, 48, 130, 70, 88, 86, 64, 50, 94, 20, 18, 88, 30, 8,
  "Rapid gravity filters",     64, 33, 3, 25, 24, 65, 35, 44, 43, 32, 25, 47, 10, 9, 44, 15, 4,
  "Clarifiers",               128, 65, 5, 50, 48, 130, 70, 88, 86, 64, 50, 94, 20, 18, 88, 30, 8,
  "Activated sludge",         396, 290, 18, 145, 229, 334, 127, 123, 204, 139, 212, NA, NA, NA, NA, NA, NA,
  "Settlement tanks",         566, 414, 25, 207, 328, 478, 182, 176, 292, 199, 303, NA, NA, NA, NA, NA, NA,
  "Screening",               3393,2484,150,1239,1965,2865,1089,1056,1749,1194,1815,NA, NA, NA, NA, NA, NA,
  "Combined public sewers", 10323,8924, 84,8410,9353,12128,2516,5812,22817,3128,16271,NA, NA, NA, NA, NA, NA,
  "Network pumping stations",6284,2519, 96, 973,1223,4782,3519,5144,2765,2172,2608,NA, NA, NA, NA, NA, NA
)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#####################################################
#####################################################
# Question 1. Sector-wide Asset Health Issue 
# (Proportion of Assets in Health Bands)
# Check whether the proportion of assets across health
# bands deviates from an expected distribution.
#####################################################
#####################################################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#####################################################
# first check using Wald method for normal approximation
#####################################################
# sample size calculation for chi squared test
confidence_level <- 0.95 # 95%
margin_error <- 0.05 # ±5%
p_hat <- 0.5 # assuming here that 50% of assets are in bands 1-2
z_score <- qnorm(1 - (1 - confidence_level) / 2)
sample_size_Wald <- (z_score^2 * p_hat * (1 - p_hat)) / (margin_error^2)
sample_size_Wald <- ceiling(sample_size)  

#####################################################
# next compute effect size with null and alternate 
# distrbution and run power analysis 
#####################################################
# assuming 5 health bands, and an expected proportion (equal distribution)
null_distribution <- c(0.2, 0.2, 0.2, 0.2, 0.2) 
alternate_distribution <- c(0.25,0.25, 0.1666667,0.1666667,0.1666667)
effect_size <- ES.w1(null_distribution,alternate_distribution) 

# power analysis
Q1_samplesize <- tibble()
for(i in 100:1000){
print(paste0("Iteration: ",i))
q1_pwr <- pwr.chisq.test(w = effect_size, df = 4, N = i, sig.level = 0.05, power = NULL)
Q1_samplesize = bind_rows(Q1_samplesize,tibble(effect_size = q1_pwr$w,
                                               sample_size = q1_pwr$N,
                                               significance_lvl = q1_pwr$sig.level,
                                               power = q1_pwr$power, 
                                               method =q1_pwr$method))
}


ggplot(Q1_samplesize %>% mutate(Title = "Q1 - sector-wide asset health power analysis")) +
  geom_line(aes(x = sample_size, y = power), linewidth = 1.5,linetype="solid") +
  facet_wrap(~ Title) +
  labs(
    x = "Sample size (sector)",
    y = "Statistical power") +
  geom_hline(yintercept = 0.8, color = "red", linetype = "dashed",linewidth = 1.75) +
  scale_x_continuous(
    breaks = seq(100, max(Q1_samplesize$sample_size), by = 100),  
    labels = seq(100, max(Q1_samplesize$sample_size), by = 100)
  ) +
  theme_minimal() +
  theme(
    text = element_text(family = "Krub"), 
    axis.text.x = element_text(angle = 0, vjust = 1, hjust = 0.5, colour = "black"),  
    strip.background = element_rect(fill = "black"),
    strip.text = element_text(colour = 'white', face = "bold", size = 12),
    panel.background = element_blank(),
    panel.grid.major.x = element_line(size = 0.5, linetype = 'solid', colour = "lightgrey"),
    panel.grid.major.y = element_line(size = 0.5, linetype = 'solid', colour = "lightgrey"),
    panel.border = element_rect(colour = "black", size = 1, fill = NA),
    axis.text.y = element_text(colour = "black"),
    legend.key = element_blank(),
    legend.background = element_blank(),
    legend.spacing = unit(2, "lines"),
    legend.spacing.x = unit(10, "pt"),
    legend.margin = margin(t = -0.5, l = 0.05, b = 0.05, r = 0.05, unit = 'cm'),
    plot.title = element_text(hjust = 0.5, size = 20, face = "bold")
  )
#ggsave(
#  filename = paste0(output_dir, "/", "WOC_WSI", "_worst_performers_plot.png"),
#  width = 9*1.1, height = 6*1.1
#)

sample_size_ChiSq <- Q1_samplesize %>%
  filter(power >= 0.8) %>%
  slice_head(n = 1) %>% pull(sample_size)

#####################################################
# next compute effect size with null and alternate 
# distrbution and run power analysis 
#####################################################

# pivot data
long_asset <- asset_data %>%
  pivot_longer(-Asset_Group, names_to = "Company", values_to = "Count") %>%
  filter(!is.na(Count))

# calculate proportionate sample sizes per company
# also ensuring that at least 5 are sampled (chi sq issues)
Q1sample_allocation <- long_asset %>%
  group_by(Asset_Group) %>%
  mutate(
    Total = sum(Count),
    Wald_Sample_Size = round((Count / Total) * sample_size_Wald),
    ChiSq_Sample_Size = round((Count / Total) * sample_size_ChiSq),
    Wald_Sample_Size = if_else(Wald_Sample_Size < 5, if_else(5 > Count, Count, 5), Wald_Sample_Size),
    ChiSq_Sample_Size = if_else(ChiSq_Sample_Size < 5, if_else(5 > Count, Count, 5), ChiSq_Sample_Size)) %>%
  ungroup()

#writexl::write_xlsx(Q1sample_allocation, paste0(output_dir,"Q1SampleSizes.xlsx"))


#####################################################
# chi square test 
#####################################################

# when you get the data, put counts in these bins
observed <- c(0, 0, 0, 0, 0)

# null hypothesis distribution 
expected <- rep(sum(observed) / length(observed), length(observed))

# chi squared test
chisq.test(observed, p = expected / sum(expected))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#####################################################
#####################################################
# Question 2. Is the issue consistent across all companies? 
# Check if the average number of assets in each health 
# band is common for all companies / aligns to the average 
# across the sector.
#####################################################
#####################################################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# N.B.- a simple way to determine this would be to run an anova for each
# band and understand if companies  are different. you'll need to run a 
# post-hoc test after the fact like tukeys. could also do a two way anova...

# assume small effect size 
effect_size <- 0.1  
num_groups_all <- 17  # Number of companies
num_groups_wasc <- 11
sample_size <- 100  # Expected number of observations per group

Q2_samplesize_all <- tibble()
for(i in 100:1000){
print(paste0("Iteration: ",i))
# Power analysis for ANOVA
q2_pwr_all <- pwr.anova.test(k = num_groups_all, n = i, f = effect_size, sig.level = 0.05, power = NULL)
Q2_samplesize_all <- bind_rows(Q2_samplesize_all,tibble(effect_size = q2_pwr_all$f,
                                                 sample_size = q2_pwr_all$n,
                                                 significance_lvl = q2_pwr_all$sig.level,
                                                 power = q2_pwr_all$power))
}

Q2_samplesize_wasc <- tibble()
for(i in 100:1000){
  print(paste0("Iteration: ",i))
  # Power analysis for ANOVA
  q2_pwr_wasc <- pwr.anova.test(k = num_groups_wasc, n = i, f = effect_size, sig.level = 0.05, power = NULL)
  Q2_samplesize_wasc <- bind_rows(Q2_samplesize_wasc,tibble(effect_size = q2_pwr_wasc$f,
                                                         sample_size = q2_pwr_wasc$n,
                                                         significance_lvl = q2_pwr_wasc$sig.level,
                                                         power = q2_pwr_wasc$power))
}

sample_size_anova_all <- Q2_samplesize_all %>%
  filter(power >= 0.8) %>%
  slice_head(n = 1) %>% pull(sample_size)*5

sample_size_anova_wascs <- Q2_samplesize_wasc %>%
  filter(power >= 0.8) %>%
  slice_head(n = 1) %>% pull(sample_size)*5

# calculate proportionate sample sizes per company
Q2sample_allocation <- long_asset %>%
  group_by(Asset_Group) %>%
  mutate(Total = sum(Count),
         Anova_Sample_Size = case_when(
           Asset_Group %in% c("Service reservoirs", "Contact tanks", "Rapid gravity filters", "Clarifiers") ~ round((Count / Total) * sample_size_anova_all),
           TRUE ~ round((Count / Total) * sample_size_anova_wascs)),
         Anova_Sample_Size = if_else(Anova_Sample_Size < 5, if_else(5 > Count, Count, 5), Anova_Sample_Size)) %>%
  ungroup()

#writexl::write_xlsx(Q2sample_allocation, paste0(output_dir,"Q1SampleSizes.xlsx"))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#####################################################
#####################################################
# Question 3. Is there a sector-wide issue for cohorts 
# (age & material)? Look at assets with a poor health 
# score and assessing if this is due to age, and 
# construction material.
#####################################################
#####################################################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# assume small effect size
effect_size <- 0.1

# predictors (e.g., age + assume 4 material levels = 5 - 1 dummy = 4 total predictors)
num_predictors <- 4

# power analysis
logit_power <- pwr.f2.test(u = num_predictors, f2 = effect_size,
                            sig.level = 0.05, power = 0.8)

# calculate min sample size
logit_sample_size <- ceiling(power_result$v + num_predictors + 1)

# double it to be conservative 
logit_sample_size = logit_sample_size*2

Q3sample_allocation <- long_asset %>%
  group_by(Asset_Group) %>%
  mutate(Total = sum(Count),
         Logit_Sample_Size = round((Count / Total) * logit_sample_size),
         Logit_Sample_Size = if_else(Logit_Sample_Size < 5, if_else(5 > Count, Count, 5), Logit_Sample_Size)) %>%
  ungroup()

#### example logit model
#asset_scores$PoorHealth <- ifelse(assets$AssetBand <= 2, 1, 0)
#asset_logit_model <- glm(PoorHealth ~ Age + Material, data = asset_scores, family = binomial)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#####################################################
#####################################################
# Question 4. Is the cohort issue consistent across 
# all companies? Want to understand how cohorts of 
# asset classes are represented consistently across 
# the sector and are not skewed by company.
#####################################################
#####################################################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# assume small effect size
effect_size <- 0.1

num_predictors <- 1 + 3 + 16

# power analysis
logit_power_company <- pwr.f2.test(u = num_predictors, f2 = effect_size, sig.level = 0.05, power = 0.8)

# calculate min sample size
logit_sample_size_company <- ceiling(logit_power_company$v + num_predictors + 1)

# double it to be conservative
logit_sample_size_company <- logit_sample_size_company * 2

Q4sample_allocation <- long_asset %>%
  group_by(Asset_Group) %>%
  mutate(
    Total = sum(Count),
    Logit_Sample_Size_Company = round((Count / Total) * logit_sample_size_company),
    Logit_Sample_Size_Company = if_else(Logit_Sample_Size_Company < 5,
                                if_else(5 > Count, Count, 5),
                                Logit_Sample_Size_Company)
  ) %>%
  ungroup()

#### example for logit model
#model_base <- glm(PoorHealth ~ Age + Material, data = asset_scores, family = binomial)
#model_Wcompany <- glm(PoorHealth ~ Age + Material + Company, data = asset_scores, family = binomial)
#anova(model_base, model_Wcompany, test = "Chisq")