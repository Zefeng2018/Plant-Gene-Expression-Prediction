# Plant Gene Expression Prediction (PGEP)


Function: This script is to predicte gene expression level based on chromain modification or transcription factor binding infomation. It was purely written by R.

**Dependence:**  
	
	optparse;  
	e1071;  
	randomForest;  
	bamsignals;  
	GenomicFeatures;  

**Usage:**  

	Rscript PGEP.R -h

**Options:**

	-b CHARACTER, --bams=CHARACTER
		Input bam files directory

	-g CHARACTER, --gtf=CHARACTER
		Input the gtf file

	-e CHARACTER, --expression=CHARACTER
		Input gene expression

	-m CHARACTER, --method=CHARACTER
		Input predictive model (LR: linear regression; RF: random forest; SVR: support vector regresison)

	-o OUTPUT, --output=OUTPUT
		output directory or prefix

	-h, --help
		Show this help message and exit
		
**Example:**
    
    Rscript PGEP.R  -b bam_file_director/ -g genome.gtf -e gene_expression.txt -o outfile_name -m LR
    
**Output file 1 (example): Predictive and orignial mesured gene expression levels (log2)**
predicted	measured	method
0.965	2.438	LR
5.107	8.860	LR
3.4475	3.070	LR
...
3.42	2.377	LR
4.91424972121951	5.54505404319212	LR
4.77015987244539	3.80631374502264

** Output file2 (example): PCC for 10-fold cross validation)**
0.79	LR
0.8	LR
0.81	LR
0.79	LR
0.81	LR
0.79	LR
0.79	LR
0.8	LR
0.8	LR
0.79	LR
