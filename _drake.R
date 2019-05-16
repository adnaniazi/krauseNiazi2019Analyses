# Import all function from the packages
Jmisc::sourceAll('R/')
source("analyses/load_packages_for_drake.R")

# make individual plans
data_download_plan = drake::code_to_plan("analyses/download_data.R")
data_consolidation_plan = drake::code_to_plan("analyses/consolidate_data.R")

# consolidate all plans into one
plan = drake::bind_plans(data_download_plan,
                         data_consolidation_plan)
drake::drake_config(plan)
