# 2026_pok.27.3a46_refpts
2026_pok.27.3a46_refpts


```
library(icesTAF)
## icesTAF::clean()  # Remove working directories (force re-install everything).
taf.bootstrap(clean = TRUE)
sourceAll()
```


## Structure of the code

Four main scripts are used: data.R, model.R, output.R, and report.R. In most cases, these scripts source other sub-scripts in a defined sequence. The script hierarchy and a brief description of their is as follows:

- **data.R**
  - **data_01_objects.R** - Creates FLStock from SAM fit (stockassessment package). Passes objects with old Eqsim settings and old reference point values.
- **model.R**
  - **model_01_settings.R** - Defines settings for the Fscans with mse package
  - **model_02_propagate_stock.R** - Uses information on the SAM estimated process residuals and historical uncertainty to generate starting FLStocks for each iteration 
  - **model_03_srr_fits.R** - Fits segmented regression to each of the starting FLStock realizations (from **model_02**), with fixed breakpoint at each iterations unique **empirical Blim**.
  - **model_04_condition_mse** - Condition objects for use in mse projection.
  - **model_05_run_mse** - Runs 3 scenarios:
      - Fscan with assessment error, and without harvest control rule (HCR) to determine MSY reference points (Fmsy, FmsyUpper, FmsyLower, Bmsy, etc.)
      - Single projection without assessment error, without HCR, and at Ftarget equal to Fmsy, to determine Fmsy Btrigger
      - Fscan with assessment error, and with HCR to determine precautionary reference point, Fpa.
- **output.R**
  - **output_01_refpts_summary.R** - Produces summary figures of Fscans: 1. landings & SSB, 2. risk probability.
- **report.R** - To generate final report with Rmarkdown (TO COME). 



## To do:

Further documentation:

 - Move setting of recruitment resampling in STF to **model_01_settings.R** from **model_05_run_mse.R** (`recyrs = 10`) 
 
 
