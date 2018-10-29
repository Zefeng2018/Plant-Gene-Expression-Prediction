#!/usr/bin/env Rscript

# 1. Main Functions and Steps:

# Functions： Read alignment files directions contaning  alignment files (*.bam file), gtf file to convert them into read counts matrix where each
# row strand for a gene, each column is a histone moidification

# example
# USAGE:  anova.r -i bam_files_director -o output file name

# parameters
# Options
# -i/--input  bam file directers
# -g/--gtf  gtf file
# -o/--output output name

options(warn = -1)

# 2. library dependence insepct, if not, install.
package_list <- c("optparse","ggplot2","e1071","randomForest","bamsignals","GenomicFeatures")

for(p in package_list){
  if(!suppressWarnings(suppressMessages(require(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))){
    install.packages(p, repos="http://cran.r-project.org")
    suppressWarnings(suppressMessages(library(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
  }
}
# another install source 
if (FALSE){
  # Bioconductor install
  source("https://bioconductor.org/biocLite.R")
  biocLite(c("bamsignals"))
}

# Clean enviroment object
rm(list=ls()) 

# Load essential packages
library(optparse)
library(ggplot2)
library(e1071)
library(randomForest)
library(bamsignals)
library(GenomicFeatures)

# Parse commond lines
if (TRUE){
  option_list <- list(
    make_option(c("-b", "--bams"), type = "character", default = NULL,
                help = "Input bam files directory", metavar = "character"),
    
    make_option(c("-g", "--gtf", type = "character",default = NULL),
                help = "Input the gtf file", metavar = 'character'),
    
    make_option(c("-e", "--expression", type = "character",default = NULL),
                help = "Input gene expression", metavar = 'character'),
    
    make_option(c("-m", "--method", type = "character",default = NULL),
                help = "Input predictive model (LR: linear regression; RF: random forest; SVR: support vector regresison)", metavar = 'character'),
    
    make_option(c("-o", "--output"), type = "character", default = NULL,
                help = "output directory or prefix")
    
    
  )
  opt_parser = OptionParser(option_list=option_list)
  opts = parse_args(opt_parser)
  # Show and verify your input information 
  print(paste("The director containing bam files is", opts$bams,  sep = ""))
  print (paste("The gtf file is ", opts$gtf,sep = ""))
  print (paste("The gene expression is ", opts$expression, sep = ""))
  print (paste("The predictive model is ", opts$method, sep = ""))
  print(paste("The output file prefix is ", opts$output, sep = ""))
}

## verify the parameters existence 
if (is.null(opts$b)){
  print_help(opt_parser)
  stop("Please input the bam file director ", call.=FALSE)
}

if (is.null(opts$g)){
  print_help(opt_parser)
  stop("Please input gtf file ", call.=FALSE)
}

if (is.null(opts$e)){
  print_help(opt_parser)
  stop("Please input gene expression file ", call.=FALSE)
}

if (is.null(opts$m)){
  print_help(opt_parser)
  stop("Please input predictive model ", call.=FALSE)
}

if (is.null(opts$o)){
  print_help(opt_parser)
  stop("Please input out file prefix ", call.=FALSE)
}

### main function
## get genes from gtf

Getgenes<-function(gtf_file){
  require(GenomicFeatures)
  tr<-makeTxDbFromGFF(gtf_file)
  genes<-genes(tr)
  return(genes)
}

genes<-Getgenes(opts$gtf)

###### calculate normalized tag counts of X(ij) for gene i, hitone mark j  in give region for gene  
#######

Normalized_tag_counts_genes<-function(region_gr,bampath){
  require(bamsignals)
  bc<-bamCount(bampath, region_gr, mapqual = 0, paired.end = c("ignore"), tlenFilter = NULL,verbose = FALSE)
  bc<-(bc*1e6*1000/sum(bc))/width(region_gr) # fpkm
  
  
  df<-data.frame(rc=bc)
  rownames(df)<-names(region_gr)
  colnames(df)<-Get_marks_from_filenames(bampath)
  return(df)
}


Get_marks_from_filenames<-function(peak_file_name){
  basename<-basename(peak_file_name)
  mark_name<-unlist(strsplit(basename,"\\."))[1]
  return(mark_name)
}

##### initiate empty  df and get all histone noralized read info

histone_df<-data.frame(matrix(nrow = length(genes)))
rownames(histone_df)<-names(genes)


file_list<-dir(path = opts$bam,full.names = TRUE,pattern = ".bam")
for (file in file_list){
  if (!endsWith(file,".bai")&!endsWith(file,".input.bam")){
    marks<-Get_marks_from_filenames(file)
    message(paste(marks," is processing...",sep = ""))
    if (marks=="DNaseI"){
      histone_df<-cbind(histone_df,Normalized_tag_counts_genes(promoters(genes,upstream = 1000,downstream = 0),file))  #gene promoters(TSS 1KB)
    }
    else{
      histone_df<-cbind(histone_df,Normalized_tag_counts_genes(genes,file))  # gene_body
    }
  }
}

histone_df<-histone_df[-1] ## all gene 

genes_expression<-read.table(opts$expression,header = TRUE,stringsAsFactors = FALSE) # (best rna-seq)

## common genes in both gtf and gene expression
common_genes<-intersect(rownames(histone_df),rownames(genes_expression))
genes_expression<-genes_expression[common_genes,]

histone_df2<-histone_df[rownames(genes_expression),]  # genes with expression used in following analysis

#apply(histone_df2,2,function(x) cor(log(genes_expression$leaf1+1),log(x+1)))


### cross validation for d2 data (add the pesudo counts derived from d1 data)
d2_histone<- log(histone_df2+1,2)
d2_exp<-data.frame(leaf=log((genes_expression$leaf1+genes_expression$leaf2+genes_expression$leaf3)/3+1,2))
d2_data<- cbind(d2_histone,d2_exp) #***####used to cross validation  (triaining and test data)

#10 cvs
nrFolds <- 10
folds <- rep_len(1:nrFolds, nrow(d2_data))
folds <- sample(folds, nrow(d2_data))  # https://stats.stackexchange.com/questions/61090/how-to-split-a-data-set-to-do-10-fold-cross-validation


######################################################################
#linear regresison
if (opts$method == "LR"){
cross_validations_cor<-c()  # record each cor
predicted_ten_fold<-c()
measured_ten_fold<-c()

for(k in 1:nrFolds) {
  fold <- which(folds == k)
  data.train <- d2_data[-fold,]
  data.test <- d2_data[fold,]
  #message(dim(data.train),"**",dim(data.test)) # 90% used to train, 10% used to validate
  
  model<-lm(formula = leaf ~ .,data = data.train) 
  predictions<-predict(model,data.test)
  #message(cor(predictions,data.test$leaf))
  predicted_ten_fold<-c(predicted_ten_fold,predictions)
  measured_ten_fold<-c(measured_ten_fold,data.test$leaf)

  cross_validations_cor<-c(cross_validations_cor,round(cor(predictions,data.test$leaf),2))
}

### output prediction accuracy

data<-data.frame(predicted=predicted_ten_fold, measured=measured_ten_fold, method=opts$method)
data1 <-data.frame(accuracy=cross_validations_cor,method=opts$method)

write.table(data,file=opts$output,row.names = FALSE,
            col.names=!file.exists(opts$output),
            sep="\t",
            quote=FALSE,
            append = TRUE)

write.table(data1,file=paste(opts$output,".cv",sep=""),row.names = FALSE,
            col.names=!file.exists(opts$output),
            sep="\t",
            quote=FALSE,
            append = TRUE)

}

 ### rnadom foresrt

if (opts$method == "RF"){
cross_validations_cor<-c()  # record each cor
predicted_ten_fold<-c()
measured_ten_fold<-c()

for(k in 1:nrFolds) {
  fold <- which(folds == k)
  data.train <- d2_data[-fold,]
  data.test <- d2_data[fold,]
  #message(dim(data.train),"**",dim(data.test)) # 90% used to train, 10% used to validate
  model<-randomForest(leaf ~ ., data=data.train,importance = TRUE, proximity = FALSE, ntree = 200)
  predictions<-predict(model,data.test)
  predicted_ten_fold<-c(predicted_ten_fold,predictions)
  measured_ten_fold<-c(measured_ten_fold,data.test$leaf)
  cross_validations_cor<-c(cross_validations_cor,round(cor(predictions,data.test$leaf),2))
}
data<-data.frame(predicted=predicted_ten_fold, measured=measured_ten_fold, method=opts$method)
data1 <-data.frame(accuracy=cross_validations_cor,method=opts$method)

write.table(data,file=opts$output,
            row.names = FALSE,
            col.names=!file.exists(opts$output),
            sep="\t",
            quote=FALSE,
            append = TRUE)

write.table(data1,file=paste(opts$output,".cv",sep=""),row.names = FALSE,
            col.names=!file.exists(opts$output),
            sep="\t",
            quote=FALSE,
            append = TRUE)

}

if (opts$method == "SVR"){
cross_validations_cor<-c()  # record each cor
predicted_ten_fold<-c()
measured_ten_fold<-c()
library(e1071)

for(k in 1:nrFolds) {
  fold <- which(folds == k)
  data.train <- d2_data[-fold,]
  data.test <- d2_data[fold,]
  #message(dim(data.train),"**",dim(data.test)) # 90% used to train, 10% used to validate
  
  model<-svm(leaf ~ ., data=data.train)
  predictions<-predict(model,data.test)
  #message(cor(predictions,data.test$leaf))
  predicted_ten_fold<-c(predicted_ten_fold,predictions)
  measured_ten_fold<-c(measured_ten_fold,data.test$leaf)
  
  cross_validations_cor<-c(cross_validations_cor,round(cor(predictions,data.test$leaf),2))
}

data<-data.frame(predicted=predicted_ten_fold, measured=measured_ten_fold,method=opts$method)
data1 <-data.frame(accuracy=cross_validations_cor,method=opts$method)

write.table(data,file=opts$output,
            row.names = FALSE,
            col.names=!file.exists(opts$output),
            sep="\t",
            quote=FALSE,
            append = TRUE)

write.table(data1,file=paste(opts$output,".cv",sep=""),row.names = FALSE,
            col.names=!file.exists(opts$output),
            sep="\t",
            quote=FALSE,
            append = TRUE)
}
