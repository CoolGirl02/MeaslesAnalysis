library(tidyverse)  # data wrangling and plots

df <- read_csv(  # read the table
  "data/raw/schools_vax_geocoded_25_26_withdates.csv",
  col_types = cols(total_students = col_character()),
  show_col_types = FALSE
)

df_1 <- df %>%
  select(school_name, type, total_students, measles_exposure, percent_immunized,
         date1, date2, date3, n_cases, as_of, county) %>%
  mutate(

    total_students = case_when(
      total_students == "<10" ~ 9L,  # treat enrolment under 10 as 9
      total_students == ""    ~ NA_integer_,
      TRUE                    ~ suppressWarnings(as.integer(total_students))
    ),

    total_students = if_else(
      is.na(total_students),
      as.integer(round(mean(total_students, na.rm = TRUE))),
      total_students
    ),

    n_cases = if_else(is.na(n_cases), 0, n_cases)  # missing case counts become zero
  )

df_filtered <- df_1 %>%  # Spartanburg schools only
  filter(county == "spartanburg")  # keep Spartanburg schools

df_left <- read_csv("data/processed/schools_data.csv", show_col_types = FALSE) %>%  # read the table
  select(school_name, unit_Id)

df_final <- right_join(df_left, df_filtered, by = "school_name", suffix = c("_left", "_right")) %>%  # schools that have both coverage and a network id
  filter(!is.na(unit_Id))

write_csv(df_final, "data/processed/school_data_joined.csv")  # save the table

df_cases <- read_csv("data/raw/SC_cases_by_week.csv", show_col_types = FALSE)  # weekly incident and cumulative cases

df_cases$date <- as.Date(df_cases$date, format = "%m/%d/%Y")  # parse the dates

write_csv(df_cases, "data/processed/cases_counts.csv")  # save the table

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

lines(df_cases$date, df_cases$inc_cases, col = "firebrick", lwd = 2)

legend(
  "topleft",
  legend = c("Cumulative Cases", "Weekly New Cases"),
  col = c("steelblue", "firebrick"),
  lwd = 2,
  bty = "n"
)
