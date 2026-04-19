# Produce plots and tables for report

# Before:
# After:

library(icesTAF)

mkdir("report")
mkdir("report/figs")

rmarkdown::render(
  input="report_01_srr_prod.Rmd",
  output_file = "report/report_01_srr_prod.html"
)


