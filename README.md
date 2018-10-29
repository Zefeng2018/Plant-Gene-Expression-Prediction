# Plant-Gene-Expression-Prediction

Prediction gene expression based on chromain modification or transcription factor


Usage: Rscript PGE.R -h

Options:
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

