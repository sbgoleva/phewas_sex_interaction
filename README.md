# phewas_sex_interaction

the sex interaction script may be used by sourcing the script while in R using the following command while in the directory where the script is located:
source("PhePheWAS_sex_interaction_analysis.R")

If using this approach, you must change the variables at the top of the script before calling it. 

Alternatively, you may comment out the variables at the top and call them into the R workspace before sourcing the script. 

Following files must be present:
1. phecodes <- PheCode table
2. pheno_file <- path to Phenotype file w/ subject ID column and phenotype column (naming does not matter and extraneous columns may be present)
3. demo_file <- Demographics file w/ subject ID column (column naming does not matter and extraneous columns may be present)
(phenotype file & demographics file may be the same file)

Other variables to specify:
pheno.name <- "Phenotype Name"
working.directory <- "path/to/directory/for/output/"
phewas.covariates.to.use <- c("specify","all","extra","columns","to","be","used","as","regression","covariates")
sex_strat_phewas.covariates.to.use <- c("specify","all","extra","columns","to","be","used","as","regression","covariates")
#sex_strat covariates will likely just exclude the Sex column
demo.subject.col.name <- "ID column name"
demo.pheno.col.name <- "Pheno Column Name"
sex.col.name <- "column name specifying sex"
male.denotation <- "Male" OR 1 OR "M", etc. 
male.denotation <- "FeMale" OR 0 OR "F", etc. 
biobank.name <- "Name/of/site/of/analysis/to/avoid/confusion/during/collaborations"
phecodes <- "path/to/R/workspace/with/phecode/table/as/variable/named/phecode.Rdata"
load(phecodes) #reads in phenotype table as variable "phecode." Manually create if need be
#phecode=read.table(phenotype/table.txt, stringsAsFactors=FALSE, header=TRUE, comment.char="", quote="")

###libraries loaded###
#ensure these packages are installed before starting
#library(PheWAS)
#library(devtools)
#library(MASS)
#library(plyr)
#library(ggplot2)
########


Script outputs:
1. 6 PheWAS files:
  PDF file with Manhattan plot of PheWAS results for 1) All subjects, 2) Females, and 3) Males
  Txt file with tabular results from PheWAS for 1) All Subjects, 2) Females, and 3) Males
2. Pheno_phewas_Biobank_sex_interaction_results.txt:
  File with each PheCode that was significant in EITHER males or females, the sex interaction P-value for the phenotype, Bonferonni significance, and additional
    descriptive information: OR, UCI, LCI in males and females
  Bonferonni value = 0.05/# of independent tests
3. pheno_categories_sex_interactions_biobank.txt:
  For each Phecode that was significant in EITHER males or females, describes how many of them were:
    1. Female biased
    2. Male biased
    3. Had no significant sex interaction
  This file is separated by category. 
  Can be paired with python plotting script to visualize the sex bias of each category. 
4. pheno_biobank_forest_plot_sex_diff_interactionp_glm.pdf
  Forest plot of each phenotype with a significant sex interaction. Shows OR and 95% CI in males and females. 
