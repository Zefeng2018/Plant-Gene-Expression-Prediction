# Plant Gene Expression Prediction (PGEP)


Function: This script is to predict gene expression level based on chromain modification or transcription factor binding infomation. It was purely written by R. 

## **Workflow of the script**
![image](https://github.com/Zefeng2018/Plant-Gene-Expression-Prediction/raw/master/workflow.png)


## **Dependence of R package:**  
	
	optparse;  
	e1071;  
	randomForest;  
	bamsignals;  
	GenomicFeatures;  

## **Usage:**  

	Rscript PGEP.R -h

## **Options:**

	-b CHARACTER, --bams=CHARACTER
		Input bam files directory

	-g CHARACTER, --gtf=CHARACTER
		Input the gtf file

	-e CHARACTER, --expression=CHARACTER
		Input gene expression file

	-m CHARACTER, --method=CHARACTER
		Input predictive model (LR: linear regression; RF: random forest; SVR: support vector regresison)

	-o OUTPUT, --output=OUTPUT
		output directory or prefix

	-h, --help
		Show this help message and exit
		
## **Command line example:**
    
    Rscript PGEP.R  -b bam_file_directory/ -g genome.gtf -e gene_expression.txt -o outfile_name -m LR
    
## **:exclamation:Input request**

_The alignment files (.bam files) contained in the file directory (-b options) shoud be named as follows:_
1. Histone modifications data from ChIP-Seq should be named as, e.g., H3K4me3.bam.
2. DHS data should be named as e.g., DNaseI.bam.
3. Transcription factors data from ChIP-Seq should be named as TF's name, e.g., AT1G22640.bam.

Notes: all the .bam files shoud be indexed using samtools under same directory before using this script.


## **Output file 1 : Predicted and original measured gene expression levels (log2)**

    predicted	    measured	    method
    0.965	2.438	LR
    5.107	8.860	LR
    3.4475	3.070	LR
    ...
    3.422	2.377	LR
    4.914	5.545	LR
    4.770	3.806	LR

## **Output file2 : PCCs from 10-fold cross validation**

    0.79	LR
    0.80	LR
    0.81	LR
    0.79	LR
    0.81	LR
    0.79	LR
    0.79	LR
    0.80	LR
    0.80	LR
    0.79	LR
## **Other files**
1. Predict_matrix.RData: this file is the data matrix (an R data object) prepared beforehand that can be directly used as the input of prediction models.
2. Arabidopsis_rice_orthologs_list.txt：This file contains the orthologous gene paris between Arabidopsis and rice.
3. example_data is the directory containing some examples of alignment files.
4. The workflow of ChIP-Seq data analysis can be found in https://zefeng2018.github.io/Plant-Gene-Expression-Prediction/

