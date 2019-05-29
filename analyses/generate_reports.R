reports_plan = drake::drake_plan(
  rna_kr_report = rmarkdown::render(knitr_in("reports/krause_niazi_et_al_rna_analysis.Rmd")),
  dna_kr_report = rmarkdown::render(knitr_in("reports/krause_niazi_et_al_dna_analysis.Rmd")),
  rna_wo_report = rmarkdown::render(knitr_in("reports/workman_et_al_rna_analysis.Rmd"))
)
  