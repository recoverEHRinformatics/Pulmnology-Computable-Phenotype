# Incident & Worsening Pulm-PASC

This repository holds code unique to the investigation of the RECOVER initiative publication
-  *Condition-specific PASC: Development and evaluation of a pulmonary computable phenotype in the RECOVER PCORnet EHR cohort*

Other repositories host code for upstream data processes. These include:

- https://github.com/calvin-zcx/pasc_phenotype (EHR preprocessing, iptw based PASC identification specific to pulmonary function)
- https://github.com/recoverEHRinformatics/recover_pasc_implementation (simplified cross-system PASC identification logic)

Table of repository contents:

1. `methods_and_display_items.Rmd`
		Sections:
			- Tables for Incidence > PPV Calculation, from baseline fev/fvc ratio to followup
			- Tables for Worsening > PPV Calculation, code defined GINA/GOLD score from raw columnar data
			- Tables for Worsening > Table I, code defined classification of worsening via a number of pathways (mutually exclusive and greedy approach)
2. `RECOVERY Pulmonary CP Code Lists_8.22.22-NDC-crosswalk.xlxs`
	- Contains ICD and RxNorm:to:NDC codes for relevant pulmonary diagnoses and medications
3. `RECOVER Pulmonary COPD&Asthma CPs_8.19.22.docx`
	- Three page specification of the computatble phenotype
4. `Recover_pulm_incidence.py`
	- Helper functions used in upstream cohort selection processes (see other repos above)
 

