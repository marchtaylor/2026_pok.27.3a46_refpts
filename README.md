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
 - Switch from `parallel::parLapply` to `mse::mps` to run Fscan?
 
## Other relevant information
 - Rationale for using truncated time series for fitting of SRR, but full time series for calculating empirical Blim: "The benchmark concluded that there is no evidence of impaired recruitment within the range of observed SSB values. Therefore, it was decided to set Blim=Bloss=130090 tonnes (Type 5 stock). At the same time the benchmark concluded that the current low productivity regime will likely continue at least for the next years. Therefore, it was decided to fit a segmented regression with a breakpoint at Bloss of the full time-series (to avoid losing information from previous parts of the time-series), but to truncate the SSB and recruitment time-series to the period from 2002 to the last assessment year, in order to take into account the lower level of productivity (Fig-ure 4.4.4). This avoids too optimistic Fp05 values that would be considerably higher if a segmented regression fitted to the full time-series were used in the simulations to estimate reference points (because of nearly doubled average recruitment estimates in absence of truncation: Figure 4.4.4, left panel). Truncation to estimate Bloss for the sake of consistency was discussed, but would – in addition to loosing possibly relevant information from former stable periods – be in contradic-tion with guidelines as this would set Blim based on recent SSB (final year) in a downward trend context (ICES, 2021a)." (WGBGAD, ICES 2024)


## References
 - ICES, 2024. Benchmark workshop on selected haddock and saithe stocks (WKBGAD). ICES Scientific Reports. https://doi.org/10.17895/ICES.PUB.25002470
 - ICES, 2021. ICES fisheries management reference points for category 1 and 2 stocks. ICES Technical Guidelines. https://doi.org/10.17895/ICES.ADVICE.7891
 




