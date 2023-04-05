<br/> **LiD genetic deteminants study.** <br/><br/>

In this remote, you will find:  


1.Workflows we used on the analysis under **./workflows/**

- Clinical data QC and model comparison in PDBP cohorts was performed in Terra platform following: [0. PDBP-clinical_qc_modelscompare.ipynb](hhttps://github.com/AMCalejandro/LID-CPH/blob/main/workflows/0.%20PDBP-clinical_qc_modelscompare.ipynb)  
- Clinical data QC for TP, OPDC, PPMI, and PD-STAT cohorts explained in Methods was performed on HPC cluster following: [1.clinical_qc.Rmd](https://github.com/AMCalejandro/LID-CPH/blob/main/workflows/1.clinical_qc.Rmd)  
- All Cohort Summary statistics after clinical and genetic QC were derived on our HPC cluster following [2.cohortSummary.Rmd](https://github.com/AMCalejandro/LID-CPH/blob/main/workflows/2.cohortSummary.Rmd)  
- Model comparison and study of clinical features associated with LID for TP, OPDC, PPMI, and PD-STAT cohorts was performed on HPC cluster following: [3. model_compare.Rmd](https://github.com/AMCalejandro/LID-CPH/blob/main/workflows/3.%20model_compare.Rmd)  
- Power calculation was performed on HPC cluster following: [4.power_estimate.Rmd](https://github.com/AMCalejandro/LID-CPH/blob/main/workflows/4.power_estimate.Rmd)  
- CPH regression models  for each cohort were derived on our HPC following: [5. run_cph_parallel.Rmd](https://github.com/AMCalejandro/LID-CPH/blob/main/workflows/5.%20run_cph_parallel.Rmd)  
- Forest plots for the GWAS sginificant lead SNP  was performed on HPC cluster following: [6. forest_plots.Rmd](https://github.com/AMCalejandro/LID-CPH/blob/main/workflows/6.%20forest_plots.Rmd)  
- KM curves for the GWAS sginificant lead SNP  as well as for known clinical features reported on literature to be associated with LID was performed on HPC cluster following: [7. kmcurves.Rmd](https://github.com/AMCalejandro/LID-CPH/blob/main/workflows/7.%20kmcurves.Rmd)  
- Sensitivity analyses were performed on HPC cluster following: [8. sensitivityAnalysis.Rmd](https://github.com/AMCalejandro/LID-CPH/blob/main/workflows/8.%20sensitivityAnalysis.Rmd)  
- Functional annotation and fine-mapping was was performed on HPC cluster following: [9. echolocator_run.Rmd](https://github.com/AMCalejandro/LID-CPH/blob/main/workflows/9.%20echolocator_run.Rmd)  
- Colocalization analyses against eQTL datasets were performed on HPC cluster following: [10.1. coloc_psychenc_etqtgen.Rmd](https://github.com/AMCalejandro/LID-CPH/blob/main/workflows/10.1.%20coloc_psychenc_etqtgen.Rmd)  & [10.2. coloc_metabrain.Rmd](https://github.com/AMCalejandro/LID-CPH/blob/main/workflows/10.2.%20coloc_metabrain.Rmd)  
- Analysis to find SNP independently associated with LID sruvival probability was perfoemd on HPC cluster following: [11. runcojo.Rmd](https://github.com/AMCalejandro/LID-CPH/blob/main/workflows/11.%20runcojo.Rmd)  
- Candidate variants analysis was performed on HPC cluster following: [12. candidate_genevariants.Rmd](https://github.com/AMCalejandro/LID-CPH/blob/main/workflows/12.%20candidate_genevariants.Rmd)  
- PRS and ROC curves were generated on HPC cluster following: [13. compute_prs.Rmd](https://github.com/AMCalejandro/LID-CPH/blob/main/workflows/13.%20compute_prs.Rmd)  - Baseline features selection on based on a stepwise regression algorithm implemented on top of logistic model regression were generated on HPC cluster following: [14. stepwiseSelection.Rmd](https://github.com/AMCalejandro/LID-CPH/blob/main/workflows/14.%20stepwiseSelection.Rmd)  
- Colocalization analyses against eQTL datasets were performed on HPC cluster following: [15.1. lidvsnonlid_motorcognitive.Rmdd](https://github.com/AMCalejandro/LID-CPH/blob/main/workflows/15.1.%20lidvsnonlid_motorcognitive.Rmd)  & [15.2. anxiety_tpd_ppmi.Rmd](https://github.com/AMCalejandro/LID-CPH/blob/main/workflows/15.2.%20anxiety_tpd_ppmi.Rmd)
<br/><br/>

2.Functions we created to perform the analysis under **./utils/**  
This directory hosts all functions we have developed during the LID study length.  
A README is added explaining what each R script contains.
<br/><br/>

3.Code to perform genetic data QC under **./genetic_qc/**  
This is the bash script we used to perform the QC at the genetic level.  
A README is added exaplaining what each bash script does.
<br/><br/>

4.Scripts to perform meta-analysis under **./metal/**  
Here we keep the metal software script we used to run meta-analysis on cohort specific GWAS CPH regression model results.
<br/><br/>

Cite the code: (TO BE ADDED AFTER DEPLOYING TO ZENODO)
