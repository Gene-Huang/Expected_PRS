# Expected_PRS

Here is the R codes for implementing our work on "The expected polygenic risk score (ePRS) framework: an equitable metric for quantifying polygenetic risk via modeling of ancestral makeup" 
(check out our preprint on medRxiv: https://medrxiv.org/cgi/content/short/2024.03.05.24303738v1)

There are three folders that contain all the primary code for the illustration purposes. (1) “compute_global_ePRS” folder contains the code and functions for computing ePRS; (2) “simulate_local_ancsestry_genetic_data” folder contains the code for generating local ancestry patterns and allele counts data for admixed populations; (3) “simulation_studies” folder contains the code for primary simulation studies. 

[compute_global_ePRS folder] The main R function for computing ePRS is in “compute_global_ePRS.R”. “example.R” provides an example for computing global ePRS (using three .csv example files). More detailed descriptions of the functions are summarized in README.docx file.

[simulate_local_ancsestry_genetic_data folder] The “example.R” provides an example for generating global ancestry proportions, local ancestry patterns, and allele counts for individuals of admixed population. 

[simulation_studies folder] The main R function for simulation studies is in “main_simulation_code.R”. In this code, we provide an example for simulating PRSs, simulating unknow confoundings, and conducting PRS-outcome association analyses. More detailed descriptions of the functions are summarized in README.docx file.


