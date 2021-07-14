#module --ignore-cache load "LLVM/.7.0.1" GCC/8.2.0 OpenMPI/3.1.4 R/3.6.0
#R
.libPaths("~/R/rlib-3.6.0")


# #####to change######
# pheno.name = "PTSD narrow"
# pheno_file = "/data/davis_lab/golevasb/Projects/PTSD_GWAS/20190801_MEGA.GRID.RACE.GEN.batch.PCs.covariates.ptsd.all.case.control.txt"
# demo_file = "/data/davis_lab/shared/phenotype_data/biovu/delivered_data/SD_wide_pull_PheWAS_covariates/20200518_MEGA.GRID.RACE.ETH.GEN.batch.PCs.age.BMI.covariates_EU.filt1.txt"
# working.directory = "/data/davis_lab/golevasb/Projects/PTSD_GWAS/PheWAS/sex_interaction"
# phewas.covariates.to.use = c("is.male", "RACE", "ETHNICITY", "PC1", "PC2", "PC3", "PC4","median_age_of_record")
# sex_strat_phewas.covariates.to.use = c("RACE", "ETHNICITY", "PC1", "PC2", "PC3", "PC4","median_age_of_record")
# demo.subject.col.name = "GRID"
# demo.pheno.col.name = "PTSD_narrow"
# sex.col.name = "is.male"
# male.denotation = 1
# female.denotation = 0
# biobank.name = "VUMC"
# phecodes = "/data/davis_lab/sealockj/projects/scripts/sex_strat_prs_phewas/mega_geno_set_phecode_table_02062020.Rdata"
# load(phecodes) #reads in phenotype table as variable "phecode." Manually create if need be
# #phecode=read.table(phenotype/table)
# ###############

###libraries to load###
library(PheWAS)
library(devtools)
library(MASS)
library(plyr)
library(ggplot2)
########



if(!dir.exists(working.directory)){
  dir.create(file.path(working.directory))
}
setwd(file.path(working.directory))



demographics=read.table(demo_file, stringsAsFactors = FALSE, header = TRUE, quote="", comment.char="")
demographics=demographics[append(demo.subject.col.name,(unlist(strsplit(phewas.covariates.to.use, ","))))]

all_pheno=read.table(pheno_file, stringsAsFactors = FALSE, header = TRUE, quote="", comment.char="")
all_pheno=all_pheno[c(demo.subject.col.name,demo.pheno.col.name)]
colnames(all_pheno) = c(demo.subject.col.name, "Phenotype")

#create sex stratified pheno files 
gender_key = demographics[c(demo.subject.col.name,sex.col.name)]

all_pheno_gender = merge(gender_key, all_pheno, by=demo.subject.col.name, all.x = FALSE)
pheno_males = all_pheno_gender[all_pheno_gender[[sex.col.name]]==male.denotation,]
pheno_females = all_pheno_gender[all_pheno_gender[[sex.col.name]]==female.denotation,]


##Running phewases in all, male and female##
###if already run phewas and sex-strat phewas, can comment out code below and read tables in:
###if you do this, each results must have column names: 'phenotype'(phecode #), 'OR', 'lci', 'uci', and 'p'
# f_results_d = read.table("/data/davis_lab/golevasb/Projects/PNES_PRS/PhePheWAS/PNES_without_epilepsy_phephewas/new_mega/new_pheno_file/eeg_cpt_and_keyword/PNES without concurrent epilepsy_females.txt", comment.char="", quote="", stringsAsFactors=FALSE, header=TRUE, fill = TRUE, sep='\t',colClasses = "character")
# m_results_d = read.table("/data/davis_lab/golevasb/Projects/PNES_PRS/PhePheWAS/PNES_without_epilepsy_phephewas/new_mega/new_pheno_file/eeg_cpt_and_keyword/PNES without concurrent epilepsy_males.txt", comment.char="", quote="", stringsAsFactors=FALSE, header=TRUE, fill = TRUE, sep='\t',colClasses = "character")
# results_d = read.table("/data/davis_lab/golevasb/Projects/PNES_PRS/PhePheWAS/PNES_without_epilepsy_phephewas/new_mega/new_pheno_file/eeg_cpt_and_keyword/PNES without concurrent epilepsy_all.txt", comment.char="", quote="", stringsAsFactors=FALSE, header=TRUE, fill = TRUE, sep='\t', colClasses = "character")

###### Run phewas on all ########

chart.title = paste0('PheWAS of ',biobank.name,'-EHR ',pheno.name)

all_pheno$Phenotype = scale(all_pheno$Phenotype)

results=phewas(phenotypes=phecode,genotypes=all_pheno, covariates=demographics[,append(demo.subject.col.name,(unlist(strsplit(phewas.covariates.to.use, ","))))], additive.genotypes=F, min.records=30, significance.threshold=c('bonferroni', 'fdr'), cores=4)

results[,1]=gsub("X", "", results[,1])

results_d = addPhecodeInfo(results)

results_d$lci = exp(results_d$beta - (1.96*results_d$SE))
results_d$uci = exp(results_d$beta + (1.96*results_d$SE))

write.table(results_d[order(results_d$p), ], file=file.path(working.directory, paste0(pheno.name, '_all', '.txt')), col.names=T, row.names=F, quote=F, sep='\t')

#label top 20 hits
x = 25
xth_most_extreme = sort(results_d$p)[x]
pdf(file=file.path(working.directory, paste0(pheno.name, '_all', '.pdf')))
p=phewasManhattan(results, OR.direction=T, annotate.level=xth_most_extreme, title=paste0(chart.title), base.labels = FALSE, annotate.size = 4)
print(p)
dev.off()

#### males ####
chart.title = paste0('PheWAS of ',biobank.name,'-EHR males with ',pheno.name)

pheno_males$Phenotype = scale(pheno_males$Phenotype)

m_results=phewas(phenotypes=phecode[phecode[[demo.subject.col.name]] %in% pheno_males[[demo.subject.col.name]],], genotypes=pheno_males, covariates=demographics[phecode[[demo.subject.col.name]] %in% pheno_males[[demo.subject.col.name]],append(demo.subject.col.name,(unlist(strsplit(sex_strat_phewas.covariates.to.use, ","))))], additive.genotypes=F, min.records=25, significance.threshold=c('bonferroni', 'fdr'), cores=4)

m_results[,1]=gsub("X", "", m_results[,1])

m_results_d = addPhecodeInfo(m_results)

m_results_d$lci = exp(m_results_d$beta - (1.96*m_results_d$SE))
m_results_d$uci = exp(m_results_d$beta + (1.96*m_results_d$SE))


write.table(m_results_d[order(m_results_d$p), ], file=file.path(working.directory, paste0(pheno.name, '_males', '.txt')), col.names=T, row.names=F, quote=F, sep='\t')

x = 20
xth_most_extreme = sort(m_results_d$p)[x]
pdf(file=file.path(working.directory, paste0(pheno.name, '_males', '.pdf')))
p=phewasManhattan(m_results, OR.direction=T, annotate.level=xth_most_extreme, title=paste0(chart.title), base.labels = FALSE, annotate.size = 4)
#saveRDS(p, "all_plot.rds")
print(p)
dev.off()


#### females ####

chart.title = paste0('PheWAS of ',biobank.name,'-EHR females with ',pheno.name)

pheno_females$Phenotype = scale(pheno_females$Phenotype)

f_results=phewas(phenotypes=phecode[phecode[[demo.subject.col.name]] %in% pheno_females[[demo.subject.col.name]],],genotypes=pheno_females, covariates=demographics[phecode[[demo.subject.col.name]] %in% pheno_females[[demo.subject.col.name]],append(demo.subject.col.name,(unlist(strsplit(sex_strat_phewas.covariates.to.use, ","))))], additive.genotypes=F, min.records=25, significance.threshold=c('bonferroni', 'fdr'), cores=4)

f_results[,1]=gsub("X", "", f_results[,1])

f_results_d = addPhecodeInfo(f_results)

f_results_d$lci = exp(f_results_d$beta - (1.96*f_results_d$SE))
f_results_d$uci = exp(f_results_d$beta + (1.96*f_results_d$SE))


write.table(f_results_d[order(f_results_d$p), ], file=file.path(working.directory, paste0(pheno.name, '_females', '.txt')), col.names=T, row.names=F, quote=F, sep='\t')


x = 20
xth_most_extreme = sort(f_results_d$p)[x]
pdf(file=file.path(working.directory, paste0(pheno.name, '_females', '.pdf')))
p=phewasManhattan(f_results, OR.direction=T, annotate.level=xth_most_extreme, title=paste0(chart.title), base.labels = FALSE, annotate.size = 4)
print(p)
dev.off()


print("done")






######### sex interaction analysis ############

###if already run phewas and sex-strat phewas, can ignore code above and read tables in:
###if you do this, each results must have column names: 'phenotype'(phecode #), 'OR', 'lci', 'uci', and 'p'
# f_results_d = read.table("/data/davis_lab/golevasb/Projects/PNES_PRS/PhePheWAS/PNES_without_epilepsy_phephewas/new_mega/new_pheno_file/eeg_cpt_and_keyword/PNES without concurrent epilepsy_females.txt", comment.char="", quote="", stringsAsFactors=FALSE, header=TRUE, fill = TRUE, sep='\t',colClasses = "character")
# m_results_d = read.table("/data/davis_lab/golevasb/Projects/PNES_PRS/PhePheWAS/PNES_without_epilepsy_phephewas/new_mega/new_pheno_file/eeg_cpt_and_keyword/PNES without concurrent epilepsy_males.txt", comment.char="", quote="", stringsAsFactors=FALSE, header=TRUE, fill = TRUE, sep='\t',colClasses = "character")
# results_d = read.table("/data/davis_lab/golevasb/Projects/PNES_PRS/PhePheWAS/PNES_without_epilepsy_phephewas/new_mega/new_pheno_file/eeg_cpt_and_keyword/PNES without concurrent epilepsy_all.txt", comment.char="", quote="", stringsAsFactors=FALSE, header=TRUE, fill = TRUE, sep='\t', colClasses = "character")

female_results2 = f_results_d[f_results_d$bonferroni=='TRUE' & !is.na(f_results_d$OR), c("phenotype")]
male_results2 = m_results_d[m_results_d$bonferroni=='TRUE' & !is.na(m_results_d$OR), c("phenotype")]
all_results2 = results_d[results_d$bonferroni=='TRUE' & !is.na(results_d$OR), c("phenotype")]
all_to_examine = append(all_results2, female_results2)
all_to_examine = append(all_to_examine, male_results2)
all_to_examine = unique(all_to_examine)


phenotypes_sig  = phecode[c(demo.subject.col.name,all_to_examine)]

####Running logistic regression on results
demog = demographics[,append(demo.subject.col.name,(unlist(strsplit(phewas.covariates.to.use, ","))))]
demo_pheno = merge(all_pheno, demog, by = demo.subject.col.name)
demo_pheno = merge(demo_pheno, phenotypes_sig, by=demo.subject.col.name, all.x=FALSE, all.y=FALSE)

#function to apply logistic regression
interaction_term_regression = function(var) {
  tryCatch({
    formula_to_use = formula(paste0(paste0("var ~ Phenotype + ",sex.col.name, " + "),paste(unlist(sex_strat_phewas.covariates.to.use), collapse=" + "),paste0(" + Phenotype * ",sex.col.name)))
    regression = glm(formula_to_use, data=demo_pheno, family="binomial")
    coef = coef(summary(regression))
    return(coef[nrow(coef),ncol(coef)])
  }, error=function(e){cat("ERROR", conditionMessage(e),"\n")})}

#get results of interaction linear regression, but only take phenotype columns (after the first few identifying columns)
interaction_term_pval = lapply(demo_pheno[(length(phewas.covariates.to.use)+3):ncol(demo_pheno)],interaction_term_regression)
##since results are in form of a list, convert to a data frame
reg_results = ldply (interaction_term_pval, data.frame)
colnames(reg_results) = c("PheCode",  "Interaction_P")
#take "X" out of results PheCodes, resulting in only numerical phecodes
reg_results[,1]=gsub("X", "", reg_results[,1])
#map phecodes to phenotypes and groups
phecode_mapped = addPhecodeInfo(reg_results)
phecode_mapped$Significant = phecode_mapped$Interaction_P < (0.05/(ncol(demo_pheno)-(length(phewas.covariates.to.use)+3)))


##add info for males and females OR and CI
female_phewas = f_results_d[f_results_d$note != ' [Error: non-varying phenotype or genotype]',]
female_phewas = female_phewas[c('phenotype', 'OR', 'lci', 'uci', 'p','description')]
female_phewas_plotting=female_phewas[female_phewas$phenotype %in% phecode_mapped$PheCode[phecode_mapped$Significant==TRUE],]
female_phewas_plotting$Sex="Female"
female_phewas$description=NULL
colnames(female_phewas) = c("PheCode", 'OR_female', 'lci_female', 'uci_female', "p_female")

male_phewas = m_results_d[m_results_d$note != ' [Error: non-varying phenotype or genotype]',]
male_phewas = male_phewas[c('phenotype', 'OR', 'lci', 'uci', 'p','description')]
male_phewas_plotting=male_phewas[male_phewas$phenotype %in% phecode_mapped$PheCode[phecode_mapped$Significant==TRUE],]
male_phewas_plotting$Sex="Male"
male_phewas$description=NULL
colnames(male_phewas) = c("PheCode", 'OR_male', 'lci_male', 'uci_male', "p_male")

phecode_mapped_plotting=rbind(female_phewas_plotting,male_phewas_plotting)

phecode_mapped = merge(phecode_mapped, female_phewas, by='PheCode', all.x=TRUE, all.y=FALSE)
phecode_mapped = merge(phecode_mapped, male_phewas, by='PheCode', all.x=TRUE, all.y=FALSE)

phecode_mapped_order = phecode_mapped[order(phecode_mapped$Interaction_P),]
write.table(phecode_mapped_order, paste0(pheno.name,"_phewas_",biobank.name,"_sd_wide_sex_interaction_results.txt"), col.names=T, row.names=F, sep="\t", quote=F)


##df with all that had no significant sex interaction
no_interaction = phecode_mapped[phecode_mapped$Significant==FALSE,]

###determine number of each category that is skewed female, male, or is neutral
categories_no_interaction = as.data.frame(table(as.factor(no_interaction$group)))
colnames(categories_no_interaction) = c('Group', 'Number_no_interaction')

##df with female skewed sex interactions
if(nrow(phecode_mapped[phecode_mapped$Significant==TRUE & (phecode_mapped$OR_female>phecode_mapped$OR_male),])>0){
  female_interaction = phecode_mapped[phecode_mapped$Significant==TRUE & (phecode_mapped$OR_female>phecode_mapped$OR_male),]
  categories_female_interaction = as.data.frame(table(as.factor(female_interaction$group)))
  colnames(categories_female_interaction) = c('Group', 'Number_female_interaction')
  categories = merge(categories_no_interaction, categories_female_interaction, by='Group', all.x = TRUE, all.y=TRUE)
}else{
  categories=categories_no_interaction
  categories$Number_female_interaction = 0
}

##df with male skewed sex interactions
if(nrow(phecode_mapped[phecode_mapped$Significant==TRUE & (phecode_mapped$OR_female<phecode_mapped$OR_male),])>0){
  male_interaction = phecode_mapped[phecode_mapped$Significant==TRUE & (phecode_mapped$OR_female<phecode_mapped$OR_male),]
  categories_male_interaction = as.data.frame(table(as.factor(male_interaction$group)))
  colnames(categories_male_interaction) = c('Group', 'Number_male_interaction')
  categories = merge(categories, categories_male_interaction, by='Group', all.x = TRUE, all.y=TRUE)
}else{
  categories$Number_male_interaction = 0
}
categories[is.na(categories)] = 0
write.table(categories, paste0(pheno.name,"_categories_sex_interactions",biobank.name,".txt"), col.names=T, row.names=F, sep="\t", quote=F)


pdf(paste0(pheno.name,"_",biobank.name,"_forest_plot_sex_diff_interactionp_glm.pdf"))
p = ggplot(data=phecode_mapped_plotting,
           aes(x = Sex,y = OR, ymin = lci, ymax = uci))+
  geom_pointrange(aes(col=Sex))+
  geom_hline(aes(fill=Sex),yintercept =1, linetype=2)+
  xlab('Phenotype')+ ylab("Odds Ratio (95% Confidence Interval)")+
  geom_errorbar(aes(ymin=lci, ymax=uci,col=Sex),width=0.2,cex=0.5)+
  facet_wrap(~description,strip.position="top",nrow=nrow(phecode_mapped[phecode_mapped$Significant==TRUE,]),scales = "free_y") +
  theme(plot.title=element_text(size=8,face="bold"),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x=element_text(face="bold"),
        axis.title=element_text(size=10,face="bold"),
        strip.text.y = element_text(hjust=0,vjust = 1,angle=180,face="bold"))+
  coord_flip()
p
dev.off()
# 
