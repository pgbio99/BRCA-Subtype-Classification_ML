#Normalization of data present in cel file format

setwd("E:/Bioteck_Project/combined result/GSE86374/cel.files")
getwd()

library(oligo)
library(affy)

list.celfiles()
celf_path <- "E:/Bioteck_Project/combined result/GSE86374/cel.files" 

data <- ReadAffy(celfile.path = celf_path)
normalized_data <- affy::rma(data)
exprs_matrix <- exprs(normalized_data)

write.csv(exprs_matrix, file = 'normalized_data_1.csv')

#Normalization of data present in text file format

setwd("E:/Bioteck_Project/combined result/GSE57297/processed files")
getwd()

library(limma)
library(arrayQualityMetrics)

files <- list.files(pattern = "*.txt")
data_list <- lapply(files, read.table, header = TRUE, sep = "\t")

bg_corrected_data_list <- lapply(data_list, function(data)
  {exprs_values <- data[,2]
  bg_corrected_values <- exprs_values-min(exprs_values, na.rm = TRUE)
  data[,2] <- bg_corrected_values
  return(data)})

common_probes <- Reduce(intersect, lapply(bg_corrected_data_list, function(x) x[,1]))
bg_corrected_data_list <- lapply(bg_corrected_data_list, function(data)
  {data[data[,1] %in%common_probes,]})

exprs_matrix_n <- do.call(cbind, lapply(bg_corrected_data_list, function(x)x[,2]))
rownames(exprs_matrix_n) <- bg_corrected_data_list[[1]][,1]
exprs_matrix_log2 <- log2(exprs_matrix_n+1)
exprs_matrix_normalized <- normalizeBetweenArrays(exprs_matrix_log2, method = "quantile")

write.table(exprs_matrix_normalized, file = "normalized_expression_matrix_n.txt", sep = "\t", quote = FALSE, col.names = NA)
