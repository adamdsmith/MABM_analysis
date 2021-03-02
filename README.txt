This file describes the general structure and content of the MABM analysis and 
power analysis directory/repository. 

If you're familiar with R/RStudio, we recommend opening the `MABM_analysis.Rproj` 
file in RStudio as the starting point. Whether in RStudio or not, the next place
to look for guidance on the general structure of the analysis and power analysis
is the `MABM_analysis_overview.Rmd` file. This file contains the code (or sources
additional code located in the `R` subdirectory) that extracted, cleaned, and 
tidied the raw data from the MABM database, extracted and constructed habitat and
survey covariates, fit and evaluated the GLMM models for each of the five focal
bat species, and used the parameters from those fits to conduct an extensive
power analysis related to evaluating the survey effort necessary to detect 
several levels of population decline.

.
├── MABM_analysis.Rproj           # RStudio R Project file. If using RStudio, open this first; if not, see MABM_analysis_overview.Rmd  
├── MABM_analysis_overview.html   # html file that can be opened for general flow of MABM data analysis  
├── MABM_analysis_overview.Rmd    # RMarkdown file that contains code for reproducing/directing entire analysis  
├── Output                        # Data products derived using code in MABM_analysis_overview.Rmd and `R` directory  
|    ├── Models                    # Species-specific GLMM fits (models can be loaded into R without re-running code)  
|    └── Raw_DB_outputs            # Necessary data exports from MABM Access Database (alternative to storing database with repo)  
├── R                             # R source code; see individual files for context/comments  
└── README.txt  
