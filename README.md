# Semantic Mnemonic Similarity Task (sMST) - Analysis pipeline
The featquery extracted fMRI data and the preprocessed behavioural data merged and analyzed in this pipeline.

Preprocessed data is uploaded to OpenNeuro, at: https://openneuro.org/datasets/ds007275
Experimental scripts available at: https://github.com/IlyesAlex/smst2_mr_experiment

## Content
### data
  - behavioural: the grand dataset combining preprocessed data of all participants from 'smst2_mr_experiment' repository.
  - covariates: all covariates administered during the two sessions 1 mont apart for the whole 39 non-excluded sample
  - fmri: percent signal changes extracted via featquery from session 1 of the sMST for the manuscript
  - stimuli: the whole stimulus set with descriptive and psycholingustic covariates (the file includes human-corrected AI-translation to English - note, that some specific connotations are 'lost in translation')

### derivatives
All figures created by the analysis pipeline, as well as table outputs from Matlab are located here.

### scripts
The main data cleaning, analysis and figure-making pipeline contained in smst_mr1_analysis.Rmd.
Additional analysis of curve fits in MATLAB are in the *.m files
