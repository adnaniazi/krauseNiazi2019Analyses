# Import all function from the packages
Jmisc::sourceAll('R/')
source("code/load_packages_for_drake.R")
source("code/generate_reports.R")

# make individual plans
data_download_plan = drake::code_to_plan("code/download_data.R")
data_consolidation_plan = drake::code_to_plan("code/consolidate_data.R")

# consolidate all plans into one
plan = drake::bind_plans(data_download_plan,
                         data_consolidation_plan,
                         reports_plan)
drake::drake_config(plan)
