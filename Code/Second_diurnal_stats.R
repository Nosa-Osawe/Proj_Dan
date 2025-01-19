# Package(s) loading...

library(emmeans)
library(nlme)
library(PhotoGEA)
library(tidyverse)

data2nd <- PhotoGEA::read_gasex_file("C:\\Users\\DELL\\Documents\\Git in R\\Proj_Dan\\Data\\diurnal second set.xlsx")
data2nd <- data.frame(data2nd$main_data)

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
    block = as.factor(block),                 
    treatment = as.factor(treatment),         
    tod = as.numeric(as.character(tod))       
  )

##############################################################################3
gsw_2nd_m1 <- lme(
  gsw ~ treatment * tod,                      
  random = ~ 1 | block/treatment,           
  correlation = corAR1(form = ~ tod | block/treatment),  
  weights = varIdent(form = ~ 1 | tod),     
  data = data2nd
)

# Summarize the model to inspect fixed and random effects
summary(gsw_2nd_m1)
anova(gsw_2nd_m1) # Interaction is significant
                  # Correlation structure is high (-0.59)

gsw_2nd_m2 <- lme(
  gsw ~ treatment * tod,                      
  random = ~ 1 | block/treatment,           
  correlation = corAR1(form = ~ tod | block/treatment),   
 # weights = varIdent(form = ~ 1 | tod),     
  data = data2nd
)

# Summarize the model to inspect fixed and random effects
summary(gsw_2nd_m2)


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
anova(Qin_2nd_m1) # interaction term is significant

    # Correlation looks weak (0.16), but not zero

Qin_2nd_m2 <- lme(
  Qin ~ treatment * tod,                     
  random = ~ 1 | block/treatment,           
  #correlation = corAR1(form = ~ tod | block/treatment),   
  weights = varIdent(form = ~ 1 | tod),     
  data = data2nd
)

summary(Qin_2nd_m2)

anova(Qin_2nd_m1, Qin_2nd_m2)  # Go with simpler model: Qin_2nd_m2
      
Qin_2nd_m3 <- lme(
  Qin ~ treatment * tod,                     
  random = ~ 1 | block/treatment,           
  #correlation = corAR1(form = ~ tod | block/treatment),   
 # weights = varIdent(form = ~ 1 | tod),     
  data = data2nd
)
summary(Qin_2nd_m3)
anova(Qin_2nd_m2, Qin_2nd_m3) # Go with more complex model: Qin_2nd_m2

emmeans(Qin_2nd_m2, pairwise ~ treatment)

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
            # Interaction term is sig.

Tleaf_2nd_m2 <- lme(
  Tleaf ~ treatment * tod,                     
  random = ~ 1 | block/treatment,           
 # correlation = corAR1(form = ~ tod | block/treatment),   
  weights = varIdent(form = ~ 1 | tod),     
  data = data2nd
)

summary(Tleaf_2nd_m2)
anova(Tleaf_2nd_m1, Tleaf_2nd_m2) # simpler model is better


Tleaf_2nd_m3 <- lme(
  Tleaf ~ treatment * tod,                     
  random = ~ 1 | block/treatment,           
  # correlation = corAR1(form = ~ tod | block/treatment),   
 # weights = varIdent(form = ~ 1 | tod),     
  data = data2nd
)

summary(Tleaf_2nd_m3)
anova(Tleaf_2nd_m3) # interaction term is n.s

anova(Tleaf_2nd_m2,Tleaf_2nd_m3)

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
  correlation = corAR1(form = ~ tod | block/treatment),   
  # weights = varIdent(form = ~ 1 | tod),     
  data = data2nd
) 

summary(VPDleaf_2nd_m3)
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
anova(PhiPS2_2nd_m1) # interaction is n.s. 

PhiPS2_2nd_m2 <- lme(
  PhiPS2 ~ treatment + tod,                     
  random = ~ 1 | block/treatment,           
  correlation = corAR1(form = ~ tod | block/treatment),   
  weights = varIdent(form = ~ 1 | tod),     
  data = data2nd
)  # Failed to converge

PhiPS2_2nd_m3 <- lme(
  PhiPS2 ~ treatment + tod,                     
  random = ~ 1 | block/treatment,           
 # correlation = corAR1(form = ~ tod | block/treatment),   
  weights = varIdent(form = ~ 1 | tod),     
  data = data2nd
)  

summary(PhiPS2_2nd_m3)
anova(PhiPS2_2nd_m3)

anova(PhiPS2_2nd_m1, PhiPS2_2nd_m3) 
# Go with the simpler model: PhiPS2_2nd_m3

emmeans(PhiPS2_2nd_m3, pairwise~treatment)

##################################################################

A_2nd_m1 <- lme(
  A ~ treatment * tod,                      
  random = ~ 1 | block/treatment,           
  correlation = corAR1(form = ~ tod | block/treatment),  
  weights = varIdent(form = ~ 1 | tod),     
  data = data2nd
)

summary(A_2nd_m1)
anova(A_2nd_m1)
            # Interaction is n.s.

A_2nd_m2 <- lme(
  A ~ treatment + tod,                      
  random = ~ 1 | block/treatment,           
  correlation = corAR1(form = ~ tod | block/treatment),  
  weights = varIdent(form = ~ 1 | tod),     
  data = data2nd
)
summary(A_2nd_m2) # correlation structure is strong here

anova(A_2nd_m1, A_2nd_m2) # go with simpler model: A_2nd_m2

A_2nd_m3 <- lme(
  A ~ treatment + tod,                      
  random = ~ 1 | block/treatment,           
  correlation = corAR1(form = ~ tod | block/treatment),  
 # weights = varIdent(form = ~ 1 | tod),     
  data = data2nd
)

summary(A_2nd_m3)
anova(A_2nd_m2, A_2nd_m3) # go with more complex model: A_2nd_m2

emmeans(A_2nd_m2, pairwise~treatment)


