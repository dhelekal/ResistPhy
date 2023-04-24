In order to run single analysis notebooks `Simu_rbf.rmd` and `gono90.rmd` set `run_mcmc <- TRUE` on line 26, then render the notebooks using `Rmarkdown`.
In order to replicate parameter recovery simulation analysis, first execute the batch simulation scripts using `Rscript -e run_multi_det_sd01.R 1 50` and `Rscript -e run_multi_det_sd01.R 1 50`.
The first argument corresponds to the simulation index to start from, inclusive, and the second is how many simulations to run. E.g., `Rscript -e run_multi_det_sd01.R 1 50` runs simulations 1 to 50 inclusive, whereas `Rscript -e run_multi_det_sd01.R 2 1` runs simulation for index 2 only.
This can be used for splitting the job into batches.

Once all simulations have been run, render the notebooks `Simu_multi_v2_rbf_sd01.rmd` and `Simu_multi_v2_rbf_sd01.rmd` using `Rmarkdown`.
