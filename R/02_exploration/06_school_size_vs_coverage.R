library(dplyr)  # data wrangling
library(tidyverse)  # data wrangling and plots

df <- read.csv("data/processed/school_data_joined.csv") %>%  # read the table
  mutate(percent_immunized = percent_immunized * 100)  # show coverage as a percent

ggplot(df, aes(x = total_students, y = percent_immunized)) +  # start the plot
  geom_point(alpha = 0.7) +  # one mark per school

  geom_smooth(method = "lm",color = "blue", fill = "grey60", alpha = 0.3) +  # fitted trend

  geom_hline(aes(yintercept = 95, linetype = "95% Herd Immunity Target"), color = "red") +  # reference line

  scale_linetype_manual(name = "Target", values = c("95% Herd Immunity Target" = "dashed")) +  # set the legend colours or line types
  labs(  # titles and axis labels
    title = "School Size vs. Vaccine Coverage",
    x = "Total Students",
    y = "Percent Immunized (%)"
  ) +
  theme_minimal()  # plain plot theme
