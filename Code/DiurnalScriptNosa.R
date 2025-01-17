# Step 1: Load and Read Data
#file_path <-create ur own specified file path
# Read the data from specified file path using the PhotoGEA package
data1 <- PhotoGEA::read_gasex_file(file_path)
data <- data.frame(data1$main_data)

# Step 2: Preprocess the Data
# Rename columns, create block identifiers, and calculate mean values for each group
data <- data %>%
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
head(data, 10)
# Ensure block, treatment, and tod are appropriately formatted
data <- data %>%
  mutate(
    block = as.factor(block),                # Convert block to a factor
    treatment = as.factor(treatment),        # Convert treatment to a factor
    tod = as.numeric(as.character(tod))      # Convert time of day to numeric
  )

# Step 4: Initial Model Fitting
# Fit a linear mixed-effects model with random intercepts for block/treatment
# and a correlation structure for temporal dependency
library(nlme)
model <- lme(
  A ~ treatment * tod,                     # Fixed effects: treatment, time of day, and their interaction
  random = ~ 1 | block/treatment,          # Random intercepts for block and treatment within block
  correlation = corAR1(form = ~ tod | block/treatment),  # Temporal correlation structure
  weights = varIdent(form = ~ 1 | tod),    # Allow heteroscedasticity across time of day
  data = data
)

# Summarize the model to inspect fixed and random effects
summary(model)

# Observation: The corAR1 parameter (temporal correlation) is effectively 0, suggesting no meaningful temporal correlation.

# Step 5: Simplify the Model (Remove corAR1)
simplified_model <- lme(
  A ~ treatment * tod,                     # Fixed effects
  random = ~ 1 | block/treatment,          # Random intercepts for block and treatment within block
  weights = varIdent(form = ~ 1 | tod),    # Heteroscedasticity across time points
  data = data
)

# Compare AIC values of the models to justify simplification
AIC(model, simplified_model)

# Summarize the simplified model
summary(simplified_model)

# Step 6: Perform ANOVA on the Simplified Model
# Evaluate the significance of fixed effects in the simplified model
anova(simplified_model)

# Observation: If tod and treatment:tod interaction are not significant, further simplify the model.

# Step 7: Final Model
# Remove non-significant terms (tod and treatment:tod)
final_model <- lme(
  A ~ treatment,                        # Fixed effects: treatment only
  random = ~ 1 | block/treatment,       # Random intercepts
  weights = varIdent(form = ~ 1 | tod), # Heteroscedasticity by time of day
  data = data
)

# Summarize the final model to inspect fixed and random effects
summary(final_model)

# Compare AIC values of the simplified and final models to ensure further simplification is justified
AIC(simplified_model, final_model)

# Step 8: Interpretation and Visualization
# Pairwise comparisons for treatment using the final model
library(emmeans)
emmeans(final_model, pairwise ~ treatment)

#I hope you gfet the idea at this point
#Let me know if you run into any trouble. I dont expect it to be overly cumbersome
