# Setup R error handling to go to stderr
options(show.error.messages=F, error=function(){cat(geterrmessage(),file=stderr());q("no",1,F)})

# We need to not crash galaxy with an UTF8 error on German LC settings.
loc <- Sys.setlocale("LC_MESSAGES", "en_US.UTF-8")

# Import library
library(Biobase)
library(knitr)
library(reshape2)
library(ggplot2)
library(monocle)
options(stringAsfactors = FALSE, useFancyQuotes = FALSE)
# Take in trailing command line arguments
args <- commandArgs(trailingOnly = TRUE)

# get options, using the spec as defined by the enclosed list.
# we read the options from the default: commandArgs(TRUE).
option_specification = matrix(c(
  'cellData', 'c', 2, 'character',
  'phenoData', 'p', 2, 'character',
  'featureData', 'f', 2, 'character',
  'lowerDetectionLimit', 2,'double',
  'expressionFamily', 2, 'character',
  'output', 'o', 2, 'character'
  ), byrow=TRUE, ncol=4);

# Parse options
options = getopt(option_specification);

# Print options to see what is going on
cat("\n cellData: ",options$cellData)
cat("\n phenoData: ",options$phenoData)
cat("\n featureData: ",options$featureData)
cat("\n lowerDetectionLimit: ",options$lowerDetectionLimit)
cat("\n expressionFamily: ",options$expressionFamily)
cat("\n output: ",options$output)

#Monocle holds single cell expresion data in objects of the CellDataSet class
#requires 3 input files:
#1. exprs: a numberic matrix of expression values, where rows are genes, and columns are cells
sample_sheet <- read.csv(options$cellData, sep=",",row.names=1)
#create CellDataSet from the expression level file
#to avoid formatting errors, make sure to convert csv to matrix
expression_data <- as.matrix(sample_sheet)

#2. phenoData: an AnnotatedDataFrame object, where rows are cells, and columns are cell attributes (cell type, culture condition, day captures, etc)
pheno_data <-read.csv(options$phenoData, sep=",",row.names=1)
pheno_data_df <- data.frame(type=pheno_data) #must be data frame object
rownames(pheno_data_df) <= colnames(expression_data) #rownames must match expression data
pd<-new('AnnotatedDataFrame',data=pheno_data_df)

#3. featureData: an AnnotatedDataFrame object, where rows are features (e.g. genes), and columns are gene attributes, such as biotype, gc content, etc.
#This would eventually be changed to get the short name of a gene from a different file - ensembl
#using long gene name as feature data for now
feature_data <- rownames(expression_data) # Get full gene name
feature_data_df <- data.frame(gene=feature_data) # Must be data frame object
rownames(feature_data_df) <- rownames(expression_data) # Rownames must match expression data
fd <- new('AnnotatedDataFrame', data = feature_data_df) 

##Create CellDataSet object
#hsc_mpp object will be used in downstream analysis
#need 'expressionFamily=negbinomial.size()' in order to use estimateDispersions() function

if(options$expressionFamily == "default"){
    expFam = negbinomial.size()
}else if(options$expressionFamily == "negbinomial"){
expFam = negbinomial()
}else if(options$expressionFamily == "tobit"){
expFam = tobit()
}else if(options$expressionFamily == "gaussianff"){
expFam = gaussianff()
    }else{
}

hsc_mpp <- newCellDataSet(expression_data, phenoData = pd, featureData=fd, lowerDetectionLimit=options$lowerDetectionLimit, expressionFamily=expFam) 

##estimate size and dispersion
hsc_mpp <- estimateSizeFactors(hsc_mpp)
hsc_mpp <- estimateDispersions(hsc_mpp)
#Removing 786 outliers
#Warning message:
#Deprecated, use tibble::rownames_to_column() instead. 

#detect genes that have a min expression of 10% and are expressed in more than ten cels
hsc_mpp <- detectGenes(hsc_mpp, min_expr = 0.1)
expressed_genes <- row.names(subset(fData(hsc_mpp), num_cells_expressed >= 10))

##show_mRNA_totals
pData(hsc_mpp)$Total_mRNAs <- Matrix::colSums(exprs(hsc_mpp))
hsc_mpp <- hsc_mpp[,pData(hsc_mpp)$Total_mRNAs < 1e6]
upper_bound <- 10^(mean(log10(pData(hsc_mpp)$Total_mRNAs)) + 2*sd(log10(pData(hsc_mpp)$Total_mRNAs)))
lower_bound <- 10^(mean(log10(pData(hsc_mpp)$Total_mRNAs)) - 2*sd(log10(pData(hsc_mpp)$Total_mRNAs)))

pdf(file="Total_mRNAs.pdf")
qplot(Total_mRNAs, data=pData(hsc_mpp), geom="density") + geom_vline(xintercept=lower_bound) + geom_vline(xintercept=upper_bound)
dev.off()