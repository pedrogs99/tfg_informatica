This directory contains all the necessary folders:

-genomes: Contains the database of phages from MillarLab and will contain 
	  the genbank files of all the genomes to be analyzed

-motifs: Should contain motif files of the motifs to analyze

-results: Destination to save results and the analysis of these results

-main.py: Main script that finds hits of a TF in the genomes of a set of phages

-analyzeResults.py: Script that feeds off main.py results and extracts stats of interest

-config.json: Specifies parameters for the code to run (e.g clade of interest)

-myenv.yml: Conda environment file (package ete3 might be an issue when importing, may be necessary to install manually)