library(tidyverse)
library(emmeans)
library(nlme)
library(PhotoGEA)


data3rd <- PhotoGEA::read_gasex_file("C:\\Users\\DELL\\Documents\\Git in R\\Proj_Dan\\Data\\diurnal third set.xlsx")
data3rd <- data.frame(data3rd$main_data)

# Step 2: Preprocess the Data
# Rename columns, create block identifiers, and calculate mean values for each group
data3rd <- data3rd %>%
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
head(data3rd, 10)
 
# Ensure block, treatment, and tod are appropriately formatted
data3rd <- data3rd %>%
  mutate(
    block = as.factor(block),                # Convert block to a factor
    treatment = as.factor(treatment),        # Convert treatment to a factor
    tod = as.numeric(as.character(tod))      # Convert time of day to numeric
  )

################################################################################

A_3rd_m1 <- lme(
  A ~ treatment * tod,                     
  random = ~ 1 | block/treatment,           
  correlation = corAR1(form = ~ tod | block/treatment),   
  weights = varIdent(form = ~ 1 | tod),     
  data = data3rd
)   

summary(A_3rd_m1)
anova(A_3rd_m1)
    # the interaction effect is not significant

A_3rd_m2 <- lme(
  A ~ treatment + tod,                     
  random = ~ 1 | block/treatment,           
  correlation = corAR1(form = ~ tod | block/treatment),   
  weights = varIdent(form = ~ 1 | tod),     
  data = data3rd
)   

summary(A_3rd_m2)
anova(A_3rd_m2)

anova(A_3rd_m1, A_3rd_m2) # Go with simpler model

A_3rd_m3 <- lme(
  A ~ treatment + tod,                     
  random = ~ 1 | block/treatment,           
 # correlation = corAR1(form = ~ tod | block/treatment),   
  weights = varIdent(form = ~ 1 | tod),     
  data = data3rd
)   

summary(A_3rd_m3)
anova(A_3rd_m2, A_3rd_m3) # Go with simpler model: A_3rd_m3

A_3rd_m4 <- lme(
  A ~ treatment + tod,                     
  random = ~ 1 | block/treatment,           
  # correlation = corAR1(form = ~ tod | block/treatment),   
 # weights = varIdent(form = ~ 1 | tod),     
  data = data3rd
)   

summary(A_3rd_m4)

anova(A_3rd_m3, A_3rd_m4) # Go with more complex model: A_3rd_m3

# Finally
emmeans(A_3rd_m3, pairwise~treatment)$emmeans
anova(A_3rd_m3)

####################################################################################


gsw_3rd_m1 <- lme(
  gsw ~ treatment * tod,                     
  random = ~ 1 | block/treatment,           
  correlation = corAR1(form = ~ tod | block/treatment),   
  weights = varIdent(form = ~ 1 | tod),     
  data = data3rd
)   

summary(gsw_3rd_m1)
anova(gsw_3rd_m1) # interaction term looks n.s.
                  # compare with a simpler model

gsw_3rd_m2 <- lme(
  gsw ~ treatment + tod,                     
  random = ~ 1 | block/treatment,           
  correlation = corAR1(form = ~ tod | block/treatment),   
  weights = varIdent(form = ~ 1 | tod),     
  data = data3rd
)   

summary(gsw_3rd_m2)

anova(gsw_3rd_m1, gsw_3rd_m2)
          # gsw_3rd_m2, the simpler model, is better

gsw_3rd_m3 <- lme(
  gsw ~ treatment + tod,                     
  random = ~ 1 | block/treatment,           
 # correlation = corAR1(form = ~ tod | block/treatment),   
  weights = varIdent(form = ~ 1 | tod),     
  data = data3rd
)   

summary(gsw_3rd_m3)
anova(gsw_3rd_m2, gsw_3rd_m3)  # go with the simpler model:gsw_3rd_m3

gsw_3rd_m4 <- lme(
  gsw ~ treatment + tod,                     
  random = ~ 1 | block/treatment,           
  # correlation = corAR1(form = ~ tod | block/treatment),   
 # weights = varIdent(form = ~ 1 | tod),     
  data = data3rd
)   

summary(gsw_3rd_m4)
anova(gsw_3rd_m3, gsw_3rd_m4) # go with the simpler model: gsw_3rd_m4

# Finally
emmeans(gsw_3rd_m4, pairwise~treatment)$emmeans
anova(gsw_3rd_m4)

###############################################################################

Qin_3rd_m1 <- lme(
  Qin ~ treatment * tod,                     
  random = ~ 1 | block/treatment,           
  correlation = corAR1(form = ~ tod | block/treatment),   
  weights = varIdent(form = ~ 1 | tod),     
  data = data3rd
)   

summary(Qin_3rd_m1)
anova(Qin_3rd_m1)# Interaction effect is n.s

Qin_3rd_m2 <- lme(
  Qin ~ treatment + tod,                     
  random = ~ 1 | block/treatment,           
  correlation = corAR1(form = ~ tod | block/treatment),   
  weights = varIdent(form = ~ 1 | tod),     
  data = data3rd
)   

summary(Qin_3rd_m2)
anova(Qin_3rd_m1, Qin_3rd_m2) # go with simpler model: Qin_3rd_m2
                              # correlation structure looks very low. 


Qin_3rd_m3 <- lme(
  Qin ~ treatment + tod,                     
  random = ~ 1 | block/treatment,           
#  correlation = corAR1(form = ~ tod | block/treatment),   
  weights = varIdent(form = ~ 1 | tod),     
  data = data3rd
)   # Failed to converge properly

Qin_3rd_m4 <- lme(
  Qin ~ treatment + tod,                     
  random = ~ 1 | block/treatment,           
  #  correlation = corAR1(form = ~ tod | block/treatment),   
  #weights = varIdent(form = ~ 1 | tod),     
  data = data3rd
)  

summary(Qin_3rd_m4)
anova(Qin_3rd_m2, Qin_3rd_m4) # go with the more complex model: Qin_3rd_m2

# Finally
emmeans(Qin_3rd_m2, pairwise~treatment)$emmeans
anova(Qin_3rd_m2)

###############################################################################

Tleaf_3rd_m1 <- lme(
  Tleaf ~ treatment * tod,                     
  random = ~ 1 | block/treatment,           
  correlation = corAR1(form = ~ tod | block/treatment),   
  weights = varIdent(form = ~ 1 | tod),     
  data = data3rd
)   

summary(Tleaf_3rd_m1)
anova(Tleaf_3rd_m1) # interaction is n.s.


Tleaf_3rd_m2 <- lme(
  Tleaf ~ treatment + tod,                     
  random = ~ 1 | block/treatment,           
  correlation = corAR1(form = ~ tod | block/treatment),   
  weights = varIdent(form = ~ 1 | tod),     
  data = data3rd
)   # Failed to converge

Tleaf_3rd_m3 <- lme(
  Tleaf ~ treatment + tod,                     
  random = ~ 1 | block/treatment,           
 # correlation = corAR1(form = ~ tod | block/treatment),   
   weights = varIdent(form = ~ 1 | tod),     
  data = data3rd
)  

summary(Tleaf_3rd_m3)
anova(Tleaf_3rd_m1, Tleaf_3rd_m3) # go with simpler model: Tleaf_3rd_m3

Tleaf_3rd_m4 <- lme(
  Tleaf ~ treatment + tod,                     
  random = ~ 1 | block/treatment,           
  # correlation = corAR1(form = ~ tod | block/treatment),   
  # weights = varIdent(form = ~ 1 | tod),     
  data = data3rd
)  

summary(Tleaf_3rd_m4)
anova(Tleaf_3rd_m3, Tleaf_3rd_m4) # go with the more complex model: Tleaf_3rd_m3

emmeans(Tleaf_3rd_m3, pairwise~treatment)$emmeans
anova(Tleaf_3rd_m3)

################################################################################

VPDleaf_3rd_m1 <- lme(
  VPDleaf ~ treatment * tod,                     
  random = ~ 1 | block/treatment,           
  correlation = corAR1(form = ~ tod | block/treatment),   
  weights = varIdent(form = ~ 1 | tod),     
  data = data3rd
)   # Failed to converge properly

VPDleaf_3rd_m2 <- lme(
  VPDleaf ~ treatment + tod,                     
  random = ~ 1 | block/treatment,           
  correlation = corAR1(form = ~ tod | block/treatment),   
  weights = varIdent(form = ~ 1 | tod),     
  data = data3rd
)   # Good!

summary(VPDleaf_3rd_m2)

VPDleaf_3rd_m3 <- lme(
  VPDleaf ~ treatment + tod,                     
  random = ~ 1 | block/treatment,           
  #correlation = corAR1(form = ~ tod | block/treatment),   
  weights = varIdent(form = ~ 1 | tod),     
  data = data3rd
)  
summary(VPDleaf_3rd_m3)

anova(VPDleaf_3rd_m2, VPDleaf_3rd_m3) # Go with simpler model: VPDleaf_3rd_m3


VPDleaf_3rd_m4 <- lme(
  VPDleaf ~ treatment + tod,                     
  random = ~ 1 | block/treatment,           
  #correlation = corAR1(form = ~ tod | block/treatment),   
 # weights = varIdent(form = ~ 1 | tod),     
  data = data3rd
)  

summary(VPDleaf_3rd_m4)
anova(VPDleaf_3rd_m3, VPDleaf_3rd_m4) # Go with simpler model: VPDleaf_3rd_m4

# Finally
emmeans(VPDleaf_3rd_m4, pairwise~treatment)$emmeans
anova(VPDleaf_3rd_m4)

###############################################################################


PhiPS2_3rd_m1 <- lme(
  PhiPS2 ~ treatment * tod,                     
  random = ~ 1 | block/treatment,           
  correlation = corAR1(form = ~ tod | block/treatment),   
  weights = varIdent(form = ~ 1 | tod),     
  data = data3rd
) 

summary(PhiPS2_3rd_m1)
anova(PhiPS2_3rd_m1) # Interaction is n.s.

PhiPS2_3rd_m2 <- lme(
  PhiPS2 ~ treatment + tod,                     
  random = ~ 1 | block/treatment,           
  correlation = corAR1(form = ~ tod | block/treatment),   
  weights = varIdent(form = ~ 1 | tod),     
  data = data3rd
) 

summary(PhiPS2_3rd_m2)
anova(PhiPS2_3rd_m1, PhiPS2_3rd_m2) # Go with simpler model: PhiPS2_3rd_m2

PhiPS2_3rd_m3 <- lme(
  PhiPS2 ~ treatment + tod,                     
  random = ~ 1 | block/treatment,           
#  correlation = corAR1(form = ~ tod | block/treatment),   
  weights = varIdent(form = ~ 1 | tod),     
  data = data3rd
) 

summary(PhiPS2_3rd_m3)
anova(PhiPS2_3rd_m2, PhiPS2_3rd_m3) # go with simpler model: PhiPS2_3rd_m3

PhiPS2_3rd_m4 <- lme(
  PhiPS2 ~ treatment + tod,                     
  random = ~ 1 | block/treatment,           
  #  correlation = corAR1(form = ~ tod | block/treatment),   
  # weights = varIdent(form = ~ 1 | tod),     
  data = data3rd
) 

summary(PhiPS2_3rd_m4)
anova(PhiPS2_3rd_m3, PhiPS2_3rd_m4) # go with the more complex model: PhiPS2_3rd_m3

# Finally
emmeans(PhiPS2_3rd_m3, pairwise~treatment)$emmeans
anova(PhiPS2_3rd_m3)


