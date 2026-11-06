library(dplyr)
library(tidyverse)

df <- read.csv("data/cleanned/school_data_joined.csv") %>%
  mutate(percent_immunized = percent_immunized * 100)

ggplot(df, aes(x = total_students, y = percent_immunized)) +
  geom_point(alpha = 0.7) +
  
  # 1. Error beam + trend line just like Claire's (method = "lm", default se = TRUE)
  geom_smooth(method = "lm",color = "blue", fill = "grey60", alpha = 0.3) +
  
  # 2. Herd immunity target line
  geom_hline(aes(yintercept = 95, linetype = "95% Herd Immunity Target"), color = "red") +
  
  scale_linetype_manual(name = "Target", values = c("95% Herd Immunity Target" = "dashed")) +
  labs(
    title = "School Size vs. Vaccine Coverage",
    x = "Total Students", 
    y = "Percent Immunized (%)"
  ) +
  theme_minimal()