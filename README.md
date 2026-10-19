# Data and code accompanying the systematic review

The corrected efficacy inputs, model scripts, saved numerical results and ROBUST domain judgments are unchanged. Certainty was assessed for all 47 primary comparisons using the CINeMA framework. All comparisons were rated Very low. These are completed manual framework assessments supported by reproducible calculations, not outputs from the official CINeMA web application or the local App's diagnostic-accuracy automatic engine. The local App's unchanged matrix preparation and generic grading helper were checked against the completed assessments; its display accepts the supplied CINeMA app input.csv.

## Reproduce the certainty assessment

Use a separate copy of this folder. With Python 3 and NumPy installed, run `python assess_certainty.py`. The script reads the supplied numerical results, trial ledger, ROBUST judgments and study_conductance.csv. It regenerates certainty_assessment.csv (all six domains, applied penalties and reasons), certainty_study_relevance.csv, study_contributions.csv and CINeMA app input.csv. The supplied conductances came from the saved primary fits; no efficacy model was refitted for this assessment. Certainty assessment methods.md records the frozen rules. The SMD interval of -0.2 to +0.2 was chosen during revision, is not a validated clinical threshold, and was not prospectively specified. Contributions are a plug-in approximation using posterior-mean heterogeneity, not a posterior distribution of contributions.

To regenerate conductances after reproducing the primary fits, run `Rscript export_conductance.R` from this folder, then run the Python script. The R script reads fits/main__<outcome>.rds and preserves the multi-arm covariance. Existing fit objects are retained in the local analysis archive and are not included in this package.

## Reproduce the efficacy analyses

Use R 4.5.2 with rjags/JAGS, coda and metafor and run, in order: `Rscript validate_model.R`, `Rscript run_analyses.R`, `Rscript network_checks.R`, and `Rscript make_diagnostics.R`. Use LANG=C and LC_ALL=C where the Windows locale otherwise prevents UTF-8 input. The scripts record fixed seeds and model settings. Do not run in the sole copy of the supplied results. No proprietary App code is distributed. Source reports remain identifiable by source_index.csv; original PDFs are not redistributed.

ROBUST overall labels were checked against the existing overall-risk algorithm using unchanged domain inputs. The research team completed manual verification of the extracted data, domain judgments and certainty ratings. AI assistance is disclosed in the manuscript. The final database search date was 5 March 2026. Certainty must be reassessed if the evidence set changes. Original submitted CSV files remain unchanged in the local archive.
