# Package(s) loading...

library(emmeans)
library(nlme)
library(PhotoGEA)


data2nd <- PhotoGEA::read_gasex_file("C:\\Users\\DELL\\Documents\\Git in R\\Proj_Dan\\Data\\diurnal second set.xlsx")
data2nd <- data.frame(data2nd$main_data)
library(tidyverse)
# Step 2: Preprocess the Data
# Rename columns, create block identifiers, and calculate mean values for each group
data2nd <- data2nd %>%
  rename(tod = obs) %>%  #"tod" (time of day)
  mutate(
    # Create block identifiers based on plot IDs
    block = as.factor(case_when(
      grepl("17", plot.id) ~ "1",
      grepl("20", plot.id) ~ "2",
      grepl("25", plot.id) ~ "3"
    ))
  ) %>%
  # Group data by time of day, plot ID, and block to summarize measurements
  group_by(tod, plot.id, block) %>%
  summarize(
    A = mean(A, na.rm = TRUE),                # Mean photosynthesis (A)
    gsw = mean(gsw, na.rm = TRUE),            # Mean stomatal conductance
    Qin = mean(Qin, na.rm = TRUE),            # Mean light intensity
    Tleaf = mean(TleafEB, na.rm = TRUE),      # Mean leaf temperature
    VPDleaf = mean(VPDleaf, na.rm = TRUE),    # Mean vapor pressure deficit
    PhiPS2 = mean(PhiPS2, na.rm = TRUE),      # Mean quantum yield
    .groups = "drop"                          # Drop grouping after summarizing
  ) %>%
  # Create treatment groups based on plot ID suffix
  mutate(
    treatment = case_when(
      grepl("S$", plot.id) ~ "heat",
      grepl("N$", plot.id) ~ "VPD",
      grepl("C$", plot.id) ~ "control"
    )
  ) %>%
  relocate(treatment, .after = plot.id)  

# Preview the processed data
head(data2nd, 10)
view(data2nd)
# Ensure block, treatment, and tod are appropriately formatted
data2nd <- data2nd %>%
  mutate(
    block = as.factor(block),                # Convert block to a factor
    treatment = as.factor(treatment),        # Convert treatment to a factor
    tod = as.numeric(as.character(tod))      # Convert time of day to numeric
  )

# Step 4: Initial Model Fitting
# Fit a linear mixed-effects model with random intercepts for block/treatment
# and a correlation structure for temporal dependency

gsw_2nd_m1 <- lme(
  gsw ~ treatment * tod,                     # Fixed effects: treatment, time of day, and their interaction
  random = ~ 1 | block/treatment,          # Random intercepts for block and treatment within block
  correlation = corAR1(form = ~ tod | block/treatment),  # Temporal correlation structure
  weights = varIdent(form = ~ 1 | tod),    # Allow heteroscedasticity across time of day
  data = data2nd
)

# Summarize the model to inspect fixed and random effects
summary(gsw_2nd_m1)

gsw_2nd_m2 <- lme(
  gsw ~ treatment * tod,                      
  random = ~ 1 | block/treatment,           
  correlation = corAR1(form = ~ tod | block/treatment),   
 # weights = varIdent(form = ~ 1 | tod),     
  data = data2nd
)

# Summarize the model to inspect fixed and random effects
summary(gsw_2nd_m2)
anova(gsw_2nd_m2)

anova(gsw_2nd_m1,gsw_2nd_m2) # The simpler model reduced the AIC (which is good),
                              # yet there is no significant difference between the 
                              # more complex and simpler model
                            # Moving on with the simpler model

emmeans(gsw_2nd_m2, pairwise ~ treatment)



###################################################################################


Qin_2nd_m1 <- lme(
  Qin ~ treatment * tod,                     
  random = ~ 1 | block/treatment,           
  correlation = corAR1(form = ~ tod | block/treatment),   
  weights = varIdent(form = ~ 1 | tod),     
  data = data2nd
)

summary(Qin_2nd_m1)
    # Correlation looks weak (0.16), but not zero
    # compared to a simper model without Temporal correlation structure

Qin_2nd_m2 <- lme(
  Qin ~ treatment * tod,                     
  random = ~ 1 | block/treatment,           
  #correlation = corAR1(form = ~ tod | block/treatment),   
  weights = varIdent(form = ~ 1 | tod),     
  data = data2nd
)

summary(Qin_2nd_m2)
anova(Qin_2nd_m2)

anova(Qin_2nd_m1, Qin_2nd_m2)  # AIC is reduced but its too small
      # I'd fall back to the more complex model

# Trying to make the model simpler by taking out heteroscedasticity across time made the model run into errors.
# can't explain why

emmeans(Qin_2nd_m1, pairwise ~ treatment)

######################################################################################


Tleaf_2nd_m1 <- lme(
  Tleaf ~ treatment * tod,                     
  random = ~ 1 | block/treatment,           
  correlation = corAR1(form = ~ tod | block/treatment),   
  weights = varIdent(form = ~ 1 | tod),     
  data = data2nd
)

summary(Tleaf_2nd_m1)
anova(Tleaf_2nd_m1)
# Correation structure= 0.2 # quite good, but lets test with a against a simpler model

Tleaf_2nd_m2 <- lme(
  Tleaf ~ treatment * tod,                     
  random = ~ 1 | block/treatment,           
 # correlation = corAR1(form = ~ tod | block/treatment),   
  weights = varIdent(form = ~ 1 | tod),     
  data = data2nd
)

summary(Tleaf_2nd_m2)
anova(Tleaf_2nd_m2)
anova(Tleaf_2nd_m1, Tleaf_2nd_m2) # simpler model is better


Tleaf_2nd_m3 <- lme(
  Tleaf ~ treatment * tod,                     
  random = ~ 1 | block/treatment,           
  # correlation = corAR1(form = ~ tod | block/treatment),   
 # weights = varIdent(form = ~ 1 | tod),     
  data = data2nd
)

summary(Tleaf_2nd_m3)
anova(Tleaf_2nd_m3) # Interaction term is not significant

anova(Tleaf_2nd_m1, Tleaf_2nd_m2,Tleaf_2nd_m3)

# Tleaf_2nd_m3 is the best

Tleaf_2nd_m4<-lme(
  Tleaf ~ treatment + tod,   # Taken off the interaction                  
  random = ~ 1 | block/treatment,           
  # correlation = corAR1(form = ~ tod | block/treatment),   
  # weights = varIdent(form = ~ 1 | tod),     
  data = data2nd
)

summary(Tleaf_2nd_m4)
anova(Tleaf_2nd_m4)

anova(Tleaf_2nd_m3,Tleaf_2nd_m4)
  # Tleaf_2nd_m4 model without interaction term is better and simpler

emmeans(Tleaf_2nd_m4, pairwise ~ treatment)

#####################################################################################

VPDleaf_2nd_m1 <- lme(
  VPDleaf ~ treatment * tod,                     
  random = ~ 1 | block/treatment,           
  correlation = corAR1(form = ~ tod | block/treatment),   
  weights = varIdent(form = ~ 1 | tod),     
  data = data2nd
)  # model failed to converge

VPDleaf_2nd_m2 <- lme(
  VPDleaf ~ treatment + tod,                     
  random = ~ 1 | block/treatment,           
  correlation = corAR1(form = ~ tod | block/treatment),   
  weights = varIdent(form = ~ 1 | tod),     
  data = data2nd
)  # Worked well--Converged!!

summary(VPDleaf_2nd_m2)
anova(VPDleaf_2nd_m2)

# Correlation Structure looks large enough (-0.53)
# interaction term is significant

VPDleaf_2nd_m3 <- lme(
  VPDleaf ~ treatment + tod,                     
  random = ~ 1 | block/treatment,           
 # correlation = corAR1(form = ~ tod | block/treatment),   
  #weights = varIdent(form = ~ 1 | tod),     
  data = data2nd
) 

anova(VPDleaf_2nd_m3)
anova(VPDleaf_2nd_m2, VPDleaf_2nd_m3)

# VPDleaf_2nd_m3 is better

emmeans(VPDleaf_2nd_m3, pairwise ~ treatment)

#################################################################################

PhiPS2_2nd_m1 <- lme(
  PhiPS2 ~ treatment * tod,                     
  random = ~ 1 | block/treatment,           
  correlation = corAR1(form = ~ tod | block/treatment),   
  weights = varIdent(form = ~ 1 | tod),     
  data = data2nd
)  

summary(PhiPS2_2nd_m1)
anova(PhiPS2_2nd_m1) # interaction is insignificant 

PhiPS2_2nd_m2 <- lme(
  PhiPS2 ~ treatment + tod,                     
  random = ~ 1 | block/treatment,           
  correlation = corAR1(form = ~ tod | block/treatment),   
  weights = varIdent(form = ~ 1 | tod),     
  data = data2nd
)  # Failed to converge

PhiPS2_2nd_m3 <- lme(
  PhiPS2 ~ treatment * tod,                     
  random = ~ 1 | block/treatment,           
 # correlation = corAR1(form = ~ tod | block/treatment),   
  # weights = varIdent(form = ~ 1 | tod),     
  data = data2nd
)  

summary(PhiPS2_2nd_m3)

anova(PhiPS2_2nd_m1, PhiPS2_2nd_m3) 
# PhiPS2_2nd_m1 which has a correlation structure is better

emmeans(PhiPS2_2nd_m1, pairwise~treatment)





