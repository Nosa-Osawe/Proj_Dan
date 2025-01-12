# Load the data
data <-read.csv ("C:\\Users\\dsani2\\OneDrive - University of Illinois - Urbana\\Desktop\\Research_Data\\Data\\2024\\tissue nutrient analysis\\tissue nutrient.csv")

r8_data_micro <- data %>%
  filter(phase == "R8", Field.Id == "seed") %>%
  pivot_longer(
    cols = c("N", "S", "P", "K", "MG", "CA"),
    names_to = "nutrient",
    values_to = "value"
  ) %>%
  select(id, nutrient, value) %>%
  group_by(id, nutrient) %>%
  summarise(mean_value = mean(value, na.rm = TRUE), .groups = "drop")

# Calculate relative change
data_relative <- r8_data_micro %>%
  group_by(nutrient) %>%
  mutate(
    control_value = mean_value[id == "control"], # Ensure mean_value is used here
    relative_change = ifelse(id != "control", (mean_value - control_value) / control_value * 100, NA)
  ) %>%
  ungroup() %>%
  select(-control_value)%>%
  filter (!is.na(relative_change))
# Print the resulting data frame
print(data_relative)


library(ggplot2)

# Create the graph
ggplot(data_relative, aes(x = nutrient, y = relative_change, fill = id)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  labs(
    title = "Relative Change of Nutrients Compared to Control",
    x = "Nutrient",
    y = "Relative Change (%)"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5)
  ) +
  scale_fill_brewer(palette = "Set2")
```
Now I will do the same with the macro nutrients:
  ```{r}
data <- read.csv("C:\\Users\\dsani2\\OneDrive - University of Illinois - Urbana\\Desktop\\Research_Data\\Data\\2024\\tissue nutrient analysis\\tissue nutrient.csv")

# Process the data for microelements
r8_data_micro <- data %>%
  filter(phase == "R8", Field.Id == "seed") %>%
  pivot_longer(
    cols = c("B", "Zn", "Mn", "Fe", "Cu"), # Update for selected elements
    names_to = "nutrient",
    values_to = "value"
  ) %>%
  select(id, nutrient, value) %>%
  group_by(id, nutrient) %>%
  summarise(mean_value = mean(value, na.rm = TRUE), .groups = "drop")

# Calculate relative change
data_relative <- r8_data_micro %>%
  group_by(nutrient) %>%
  mutate(
    control_value = mean_value[id == "control"], # Ensure mean_value is used here
    relative_change = ifelse(id != "control", (mean_value - control_value) / control_value * 100, NA)
  ) %>%
  ungroup() %>%
  select(-control_value) %>%
  filter(!is.na(relative_change))

# Print the resulting data frame
print(data_relative)

# Load ggplot2 library
library(ggplot2)

# Create the graph
ggplot(data_relative, aes(x = nutrient, y = relative_change, fill = id)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  labs(
    title = "Relative Change of Microelements Compared to Control",
    x = "Nutrient",
    y = "Relative Change (%)"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5)
  ) +
  scale_fill_brewer(palette = "Set2") 