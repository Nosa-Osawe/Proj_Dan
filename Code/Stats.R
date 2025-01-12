require(tidyverse)
require(agricolae)

plant_data <- read.csv("C:\\Users\\DELL\\Documents\\Git in R\\Proj_Dan\\Data\\tissue nutrient.csv")
View(plant_data)

colnames(plant_data)

#----- Convert the N, S, P, K, Mg, Ca, and Na to mg/kg by multiplying by 100

plant_data <- plant_data %>% 
  mutate(N= N*100,
         S=S*100,
         P=P*100,
         K=K*100,
         MG=MG*100,
         CA=CA*100,
         NA.=NA.*100) %>% 
  rename(Ca=CA,
         Mg=MG,
         Na=NA.)   # Careful running more than once!!!!

#-- Use the nutrient yield values by
#---multiplying the concentration of each element by the tissue biomass yield

plant_Nutr.yield<-plant_data %>% 
  mutate(across(
    where(is.numeric) & !c("biomass"),
    ~ . * biomass,
    .names = "{.col}"
  ))

head(plant_Nutr.yield)

# - using data in R8 and Seed

R8_seed<- plant_Nutr.yield %>% 
  filter(Field.Id=="seed", 
         phase=="R8") %>% 
  select(-c(Field.Id,phase)) %>% 
  as.data.frame()

head(R8_seed)

# sample analysis
summary(N_anova_lm <- lm(N~id, data=R8_seed))
summary(N_anova <- aov(N~id, data=R8_seed))
HSD.test(N_anova, 
         trt = c("id"), group = TRUE)$groups



lm_function <- function(dependent_var, data) {

  formula <- as.formula(paste(dependent_var, "~ id"))
  lm_summary <- summary(lm_model <- lm(formula, data = data))
  aov_summary <- summary(aov_model <- aov(formula, data = data))
  hsd_groups <- agricolae::HSD.test(aov_model, trt = c("id"), group = TRUE)$groups
  
  return(list(
    lm_summary = lm_summary,
    aov_summary = aov_summary,
    hsd_groups = hsd_groups
  ))
}

print(result.N <- lm_function("N",R8_seed))
print(result.S <- lm_function("N", R8_seed))
print(result.P<-lm_function("P", R8_seed))
print(result.K<-lm_function("K", R8_seed))
print(result.Mg<-lm_function("Mg", R8_seed))
print(result.Ca <- lm_function("Ca", R8_seed))
print(result.Na <- lm_function("Na", R8_seed))
print(result.B <- lm_function("B", R8_seed))
print(result.Zn <- lm_function("Zn", R8_seed))
print(result.Mn <- lm_function("Mn", R8_seed))
print(result.Fe <- lm_function("Fe", R8_seed))
print(result.Cu <- lm_function("Cu", R8_seed))
print(result.Al.with.ppm.uom <- lm_function("Al.with.ppm.uom",
                                            R8_seed))

########################PLOTS###################################################

# ---- make the plots here
colnames(R8_seed)

R8_seed_summary <- R8_seed %>% 
  select(-biomass) %>%          # I Do not need biomass anymore
  rename(Al="Al.with.ppm.uom") %>% 
  group_by(id) %>% 
  summarise(across(where(is.numeric), mean, na.rm = TRUE)) %>% 
  as.data.frame()

view(R8_seed_summary)




R8_seed_relative <- R8_seed_summary %>%
  pivot_longer(
    cols = -id, # All columns except `id`
    names_to = "Nutrients", 
    values_to = "Value"   
  ) %>% 
  pivot_wider(
    names_from = id,
    values_from = Value
  ) %>% 
  rename(vpd="vpd ") %>% # this spacing in the name gave issues, now resolved
  mutate(Heat= (((heat-control)/control))*100) %>% 
  mutate(VPD= (((vpd -control)/control))*100) %>% 
  select(-control,-heat, -vpd) %>% 
  pivot_longer(
    cols = -Nutrients,
    names_to = "id",
    values_to = "Mean_Relative_change"
  ) %>% 
  as.data.frame()



view(R8_seed_relative)

#--- based of a quick skim of your codes, I suggest you want some elements on one plot, 
#--- and others in another

R8_seed_relative %>% 
  filter(Nutrients %in% c("N", "S", "P", "K", "Mg", "Ca", "Na")) %>% 
ggplot(aes(x = Nutrients, y = Mean_Relative_change, fill = id)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  labs(
    title = "Relative Change of Nutrients Compared to Control",
    x = "Nutrient",
    y = "Relative Change (%)",
    fill = "I.D"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5)
  ) +
  scale_fill_brewer(palette = "Set2")


R8_seed_relative %>% 
  filter(!Nutrients %in% c("N", "S", "P", "K", "Mg", "Ca", "Na")) %>% 
  ggplot(aes(x = Nutrients, y = Mean_Relative_change, fill = id)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  labs(
    title = "Relative Change of Nutrients Compared to Control",
    x = "Nutrient",
    y = "Relative Change (%)",
    fill = "I.d"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5),
  ) +
  scale_fill_brewer(palette = "Set2")


