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

