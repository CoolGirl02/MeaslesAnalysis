# Measles transmission across Spartanburg County schools

Master's project on measles spread among schools in Spartanburg County, South Carolina. The analysis joins school-level MMR coverage and enrollment to a contact network, then simulates a metapopulation SEIR model (susceptible, exposed, infectious, recovered). The simulations ask how outbreak size changes with vaccination coverage, school size, where the index case starts, and how doses or movement restrictions are allocated.

Open `MeaslesAnalysis.Rproj` and run scripts from the project root. Paths are relative to that root.

## Layout

| Path | What it holds |
| --- | --- |
| `data/raw/` | Source tables: geocoded school vaccination records, the pairwise connectivity file, its variable dictionary, and weekly case counts |
| `data/processed/` | Analysis tables built from the raw files (89 Spartanburg schools and their edges) |
| `R/01_data/` | Build the school list and network, then join coverage, enrollment, and cases |
| `R/02_exploration/` | Coverage summaries, network centrality, the geographic network, and school size versus coverage |
| `R/03_models/` | SEIR simulations, from within-school transmission to between-school coupling |
| `R/04_sensitivity/` | Sensitivity analyses: coverage, hubs versus peripheral schools, vaccine allocation, and travel restrictions |
| `figures/` | Saved plots, grouped by exploration, network, simulation, sensitivity, and supplement |
| `output/tables/` | Result tables written by the sensitivity scripts |
| `archive/` | Draft scripts and untitled plot iterations. These are not the pipeline |

## Data

`data/raw/schools_vax_geocoded_25_26_withdates.csv` is the school vaccination extract (enrollment, percent immunized, measles exposure, reported cases, coordinates, county).

`data/raw/spartanburg_connectivity_data.csv` is a school-to-school edge list. Each row has the two schools, `total_contacts`, and `connect_risk`. Column definitions are in `data/raw/connectivity_data_dictionary.csv`.

`data/raw/SC_cases_by_week.csv` is weekly incident and cumulative measles cases (`date`, `inc_cases`, `cum_cases`).

Processed files used by the models:

- `schools_data.csv` — one row per school, taken from the connectivity file
- `network_df.csv` — undirected edges (`unit_id_A`, `unit_id_B`, `total_contacts`, `connect_risk`)
- `school_data_joined.csv` — Spartanburg schools with a network id, enrollment, MMR coverage, and case fields (89 schools)
- `cases_counts.csv` — the weekly case series with dates parsed
- `schools_by_quadrant.csv` and `school_network_roles_weighted.csv` — derived geography and network-role tables

`02_join_school_and_case_data.R` keeps Spartanburg County, treats enrollment coded `<10` as 9, fills missing enrollment with the county mean, and sets missing case counts to 0. Coverage is stored as a proportion in the joined file. Model scripts still rescale values above 1, so a percentage coded as `95` is read as `0.95`.

## How to run

Install the packages the scripts load: `tidyverse`, `deSolve`, `igraph`, `ggraph`, and `scales`.

Run the data scripts first, in order, from the project root:

1. `R/01_data/01_build_network.R` writes `schools_data.csv` and `network_df.csv`.
2. `R/01_data/02_join_school_and_case_data.R` writes `school_data_joined.csv` and `cases_counts.csv`.

Then the exploration and model scripts can be run on their own. They read `data/processed/`.

| Script | Role |
| --- | --- |
| `03_coverage_eda.R` | MMR summary, herd-immunity threshold, map, and quadrant plots. Saves into `figures/exploration/` |
| `04_network_characteristics.R` | Degree, strength, betweenness, closeness, and eigenvector centrality, plus a geographic network plot |
| `05_geographic_network.R` | Geographic network coloured by MMR coverage and exposure |
| `06_school_size_vs_coverage.R` | Enrollment versus percent immunized, with a 95% coverage reference line |
| `07_seir_within_schools.R` | Transmission inside each school |
| `08_seir_within_and_between.R` | Within-school transmission plus coupling along network edges |
| `09_seir_contacts.R` | Same two-scale model weighted by recorded contacts |
| `10_seir_global.R` | County-wide SEIR on the school network |
| `11_seir_ordered.R` | School curves ordered for comparison |
| `12_seir_baseline.R` | Baseline network SEIR used as the starting point for later runs |
| `13_sensitivity_analyses.R` | Main sensitivity set: coverage, hub versus peripheral seeding, lower connectivity, targeted versus random vaccination, allocation by school size, reactive closure, ring versus mass vaccination, and phased travel restrictions |
| `14_sensitivity_coverage.R` | Coverage and within-school transmission sensitivity |
| `15_vaccine_allocation.R` | Size-based allocation against uniform coverage |
| `16_network_seir_allocation.R` | Network SEIR plus the allocation comparison |
| `17_network_seir_allocation_revised.R` | Revised allocation run. Writes the vaccine-allocation tables and `figures/simulation/measles_epidemic_curves.png` |
| `18_vaccine_allocation_revised.R` | Threshold sensitivity for the same allocation question |
| `19_network_seir_framework.R` and `20_network_seir_framework_revised.R` | Hub versus peripheral seeding on the network |
| `22_multi_hub_vs_peripheral.R` | Three hub schools versus three peripheral schools |

Later numbers revise earlier ones. For a single pass through the sensitivity questions, use `13_sensitivity_analyses.R`. For the allocation tables in `output/tables/`, use `17_network_seir_allocation_revised.R`. For hub-versus-peripheral tables, use `20_network_seir_framework_revised.R` and `22_multi_hub_vs_peripheral.R`.

`21_travel_restrictions.R` does not reload the CSVs. It expects a prepared school table `df` (with `total_students` and `vax_rate`) and an adjacency matrix `A` in the same school order, which the network SEIR scripts create. It compares no reduction, a 50% reduction, and a 90% reduction in inter-school connectivity from day 14, and writes the sensitivity-3 tables and plots.

Baseline transmission settings in the simulation scripts are an infectious period of about 7 days, an incubation period of about 10 days, within-school R0 of 15, and a small between-school coupling (often 5%). Check the parameter block at the top of the script you run; later sensitivity scripts change coverage and coupling on purpose.

New tables are written to `output/tables/`. New sensitivity and epidemic figures are written under `figures/sensitivity/` and `figures/simulation/`.

## Archive

`archive/scripts/` holds fragments that continue an interactive session (`2.R`, `aggregateplot.R`, `more_seir_fragment.R`, `travel_restrictions_fragment.R`, and the `recoreco*` files) and earlier drafts of the network build, exploration, and allocation code. `archive/figures/` holds working plots that were not part of a named figure. Nothing in `archive/` is required to reproduce the tables in `output/tables/`.

Identical copies were removed while sorting the project: the duplicate connectivity extract `spartanburg_connectivity_data(1).csv`, and plot files that matched another file byte for byte.
