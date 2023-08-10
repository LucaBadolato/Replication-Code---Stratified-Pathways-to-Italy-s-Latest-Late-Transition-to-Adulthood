# Replication material for *Stratified Pathways to Italy’s “Latest-Late” Transition to Adulthood*

**Date of the last update**: 2023-07-26

**Title**: Stratified Pathways to Italy’s “Latest-Late” Transition to Adulthood

**Author**: Badolato Luca (badolato.3@osu.edu), Department of Sociology and Institute for Population Research, The Ohio State University

**Journal:** Advances in Life Course Research

**DOI**: [https://doi.org/10.1016/j.alcr.2023.100563](https://doi.org/10.1016/j.alcr.2023.100563)

**Abstract**:
During the last few decades, Western societies have undergone substantial social and demographic changes, and the transition to adulthood progressively moved from an early, contracted, and simple pattern to a late, protracted, and complex one. These trends have been extensively analyzed under the Second Demographic Transition framework, emphasizing the role of individual agency and ideational change. A growing parallel literature underlines social stratification, the gender revolution, and contextual opportunities as driving forces. This paper builds on this emerging literature to analyze trends of the transition to adulthood in Italy, a salient social and demographic context among the “lowest-low” fertility countries. Drawing from the European Social Survey 2018 data, we use Sequence Analysis to compute a taxonomy of ideal types of transition to adulthood and analyze their evolution across cohorts. Our analyses show that the emergence of a late and protracted transition to adulthood, associated with “lowest-low” fertility levels, is stratified by gender and socioeconomic background. We contribute to the growing literature on the social stratification of life course trajectories and the relevance of contextual opportunities and constraints by analyzing the transition to adulthood in a low-opportunity context from a longitudinal, stratified perspective.

## Folder structure:

**Note:** The raw European Social Survey (ESS) data used in this study are publicly available and can be downloaded [here](https://ess-search.nsd.no/en/study/bdc7c350-1029-4cb3-9d5e-53f668b8fa74). 

* **data_cleaning.do** - Script to clean the data, select the variable of interest, and save the cleaned data as *data.dta*. 
   
* **sequence_analysis.R** - Script to run the Sequence Analysis and Clustering to replicate Figures 1, 2, 3, A1, A2, A3 and Tables 4, A1, A4. Takes as input *data.dta* and returns as output *data_cluster_optimal_matching.dta* and *data_cluster_optimal_matching_gender_specific.dta*.

* **regressions_analyses.do** - Script to compute the regression analyses to replicate Figures 4, A4 and Tables 5, A2, A3, A5, A6. Takes as input *data.dta*, *data_cluster_optimal_matching.dta*, and *data_cluster_optimal_matching_gender_specific.dta*. 
