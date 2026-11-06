

library(tidyverse)


#  Cleaning  and joining school data


df <- read_csv(
  "data/schools_vax_geocoded_25_26_withdates.csv", 
  col_types = cols(total_students = col_character()),
  show_col_types = FALSE
)

# Picking needed columns and fix student numbers
df_1 <- df %>%
  select(school_name, type, total_students, measles_exposure, percent_immunized, 
         date1, date2, date3, n_cases, as_of, county) %>%
  mutate(
    # Converting text numbers to real numbers
    total_students = case_when(
      total_students == "<10" ~ 9L,
      total_students == ""    ~ NA_integer_,
      TRUE                    ~ suppressWarnings(as.integer(total_students))
    ),
    # Filling in missing student counts with the average count
    total_students = if_else(
      is.na(total_students),
      as.integer(round(mean(total_students, na.rm = TRUE))),
      total_students
    ),
    # Changing missing case numbers to 0
    n_cases = if_else(is.na(n_cases), 0, n_cases)
  )

# Keep only schools in Spartanburg county
df_filtered <- df_1 %>%
  filter(county == "spartanburg")

# Read school ID data
df_left <- read_csv("data/cleanned/schools_data.csv", show_col_types = FALSE) %>%
  select(school_name, unit_Id)

# Merging datasets and remove rows with missing IDs
df_final <- right_join(df_left, df_filtered, by = "school_name", suffix = c("_left", "_right")) %>%
  filter(!is.na(unit_Id))

write_csv(df_final, "data/cleanned/school_data_joined.csv")

###############

#Plotting measles cases over time

# Read weekly measles case counts
df_cases <- read_csv("data/SC_cases_by_week.csv", show_col_types = FALSE)


df_cases$date <- as.Date(df_cases$date, format = "%m/%d/%Y")


write_csv(df_cases, "data/cleanned/cases_counts.csv")

# Plot total cumulative cases over time
plot(
  df_cases$date, df_cases$cum_cases,
  type = "l",
  col = "steelblue",
  lwd = 2,
  xlab = "Date",
  ylab = "Number of Cases",
  main = "Measles Cases",
  ylim = c(0, max(df_cases$cum_cases, na.rm = TRUE))
)

# Add weekly new cases line to the plot
lines(df_cases$date, df_cases$inc_cases, col = "firebrick", lwd = 2)

# Add legend to the top left
legend(
  "topleft",
  legend = c("Cumulative Cases", "Weekly New Cases"),
  col = c("steelblue", "firebrick"),
  lwd = 2,
  bty = "n"
)






