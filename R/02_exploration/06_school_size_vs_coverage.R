library(dplyr)  # data wrangling
library(tidyverse)  # data wrangling and plots

df <- read.csv("data/processed/school_data_joined.csv") %>%  # read the table
  mutate(percent_immunized = percent_immunized * 100)  # show coverage as a percent

fit <- lm(percent_immunized ~ total_students, data = df)  # straight line of coverage against school size
coef_ci <- confint(fit, level = 0.95)  # 95% confidence interval for intercept and slope
coef_table <- tibble(
  term = c("intercept", "slope_per_student"),
  estimate = coef(fit),
  ci_low = coef_ci[, 1],
  ci_high = coef_ci[, 2]
)

predicted <- as.data.frame(predict(fit, newdata = df, interval = "confidence", level = 0.95))  # 95% CI for the mean coverage
mean_ci <- bind_cols(df, predicted) %>%
  arrange(total_students) %>%  # sort so the confidence band is drawn left to right
  transmute(
    school_name,
    total_students,
    percent_immunized,
    fitted = fit,  # mean coverage on the line
    ci_low = lwr,  # lower end of the 95% confidence interval
    ci_high = upr  # upper end of the 95% confidence interval
  )

print(coef_table)
print(mean_ci)

dir.create("output/tables", showWarnings = FALSE)
write.csv(coef_table, "output/tables/school_size_coverage_coefficient_ci.csv", row.names = FALSE)  # save the table
write.csv(mean_ci, "output/tables/school_size_coverage_mean_ci.csv", row.names = FALSE)  # save the table

print(
ggplot(df, aes(x = total_students, y = percent_immunized)) +  # start the plot
  geom_point(alpha = 0.7) +  # one mark per school
  geom_ribbon(
    data = mean_ci,
    aes(x = total_students, ymin = ci_low, ymax = ci_high),
    fill = "grey60",
    alpha = 0.3,
    inherit.aes = FALSE
  ) +  # 95% confidence interval around the mean line
  geom_line(
    data = mean_ci,
    aes(x = total_students, y = fitted),
    color = "blue",
    linewidth = 1,
    inherit.aes = FALSE
  ) +  # fitted trend
  geom_hline(aes(yintercept = 95, linetype = "95% Herd Immunity Target"), color = "red") +  # reference line
  scale_linetype_manual(name = "Target", values = c("95% Herd Immunity Target" = "dashed")) +  # set the legend colours or line types
  labs(
    title = "School Size vs. Vaccine Coverage",
    subtitle = "Grey band is the 95% confidence interval for the mean coverage",
    x = "Total Students",
    y = "Percent Immunized (%)"
  ) +
  theme_minimal()  # plain plot theme
)
