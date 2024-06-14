#!/bin/bash

#Software: PLINK

#First starting with a logistic regression using APOE haplotypes as interaction covariates.


plink   --assoc fisher   --bfile final-vars_IPDGC_hc   --covar covar_APOE_plink._PCA2.txt  --covar-number 2,3,5-11  --exclude plink_APOE_EXCLUDE.txt  --interaction   --logistic genotypic  --out interactiontest_6.05_1  --parameters 1,2,3,4,5,6,7,8,9,10,11,16,17


#To explore all potential interactions between APOE and SNCA, expanded the range to include all SNPs in APOE and SNCA, plus 100kb outside the region ends on both genes.
#First, removed SNPs that were in LD with an r^2 of 0.5

#Used PLINK's epistatis regression-based method to model SNP x SNP interactions

plink --file IPDGC_GeneGeneInter_LD_Pruned --pheno pheno_IPDGC.txt --covar covar_APOE.txt --epistasis --epi1 .05 --threads 64 --all-pheno --out IPDGC_GeneGeneInter_LD_Pruned_Epistasis 


#In contrast to SNP-level approaches that have limited power, gene-level testing can help detect relevant and significant interactions. 
#Gene based tests can be classified into 2 categories 
#1) Tests that consider multiple markers in a gene as part of a joint model and 
#2) tests that combine marker-based test statistics or p-values into a gene-based equivalent. 


#Both categories of gene-based tests were run using the GeneGeneInteR package in R.
#GeneGeneInteR is part of the BioConductor Package
#The output is limited to statistically significant results only. In the analysis so far, we did not find any.

#The R script is in GeneGeneInteR folder > Scripts. There are 10 methods used--you can use only a handful at a time because they are computationally intensive. 


#Software: R--GeneGeneInteR Package
#Some of the methods, especially KCCA -kernal canonical correlation analysis will take ALOT of time. 40 cores with 128GB memory on the cluster took ~12 hours. Some methods may crash 

module load nixpkgs/16.09 gcc/7.3.0 r/3.6.1


if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("GeneGeneInteR")

Rscript GeneGeneInteR_LD_Pruned.R


#Lastly, a model-based Multifactor Dimensionality Reduction method was implemented as a final screen strategy to detect any significant 2-way SNP x SNP and 3-way  SNP x SNP x SNP interactions using the open-source MB-MDR software.
#You can download it here: http://www.statgen.ulg.ac.be/software_mbmdr.html OR find it in the MB-MDR > Files > mbmdr-4.4.1-linux-64bits.out . The manual is included in the folder.
#MB-MDR has parallel processing capabilities but I have not been able to implement them. For our sample size and the parameters below, CPU running time was ~6 hours. Altering 


#Software: MB-MDR

#This creates the mbmdr file --here it is mbmdrFile7

./mbmdr-4.4.1-linux-64bits.out --plink2mbmdr --binary -ped 'IPDGC_GeneGeneInter_LD_Pruned.ped' -map 'IPDGC_GeneGeneInter_LD_Pruned.map' -phe 'covar_APOE_MBMDR_PCA.phe.' -o 'mbmdrFile7' -tr 'trFile7'

#SNPxSNP
./mbmdr-4.4.1-linux-64bits.out --binary -n 10000 -p 1000 -mt gammaMAXT -x .05 -rc RESIDUALS  -ac 14 -o mbmdr4_LD_PCA_0505 -o2 mbmdr4_LD_PCA_0505_model -v LONG -pb NORMAL 'mbmdrFile7'

#SNPxSNPxSNP
./mbmdr-4.4.1-linux-64bits.out --binary -n 10000 -p 1000 -mt gammaMAXT -x .05 -rc RESIDUALS  -d 3D -ac 14 -o mbmdr7_LD_PCA_0505 -o2 mbmdr7_LD_PCA_0505_model -pb NORMAL 'mbmdrFile7'