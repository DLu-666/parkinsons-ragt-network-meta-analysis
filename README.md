# Robot-assisted gait training for Parkinson's disease

This repository contains the data and code for our systematic review and network meta-analysis. We compared robot-assisted gait training with other rehabilitation programmes and usual care in people with Parkinson's disease. The review covered motor impairment, balance, mobility, freezing of gait, quality of life, and walking performance. The database search ended on 5 March 2026.

We have included the extracted study data, the scripts used for the analyses, and the result tables. If you only want to look at the findings, you can open the CSV files directly. You do not need to run the code.

## Where to start

- `source_index.csv` lists the study reports used in the review.
- `Baseline characteristics.csv` describes the participants and treatment programmes.
- `primary_endpoints.csv` contains the data used in the main analyses. `all_timepoints.csv` also includes measurements taken at other time points.
- `comparative_effects.csv` contains the treatment comparisons and their credible intervals.
- `ROBUST_RCT_by_outcome.csv` contains the risk-of-bias assessments. `certainty_assessment.csv` contains the certainty assessments, with the reasons for each judgement.
- `Certainty assessment methods.md` explains how the certainty assessments were made.

The eAddenda accompanying the article give the source pages, data conversions, and further details needed to read these files alongside the study reports. The original articles are not included in this repository.

## Running the analyses

Download the repository and work in a copy of the folder. The scripts write output files, so keeping the supplied results in a separate copy makes it easier to compare them with a new run.

The analyses used R 4.5.2, JAGS, and the R packages `rjags`, `coda`, and `metafor`. `sessionInfo.txt` records the software environment. Run these commands from the repository folder, in this order:

```sh
Rscript validate_model.R
Rscript run_analyses.R
Rscript network_checks.R
Rscript make_diagnostics.R
```

The scripts contain the model settings and random seeds. The saved model objects are not included; they are created when the analyses run. On Windows, set `LANG=C` and `LC_ALL=C` if the system locale causes problems reading UTF-8 files.

To recalculate the certainty assessment from the supplied results, install Python 3 and NumPy, then run:

```sh
python assess_certainty.py
```

This step uses the study contributions supplied in `study_conductance.csv`. To calculate those contributions again after running the main analyses, first run `Rscript export_conductance.R`, then run the Python command above.
