# !/usr/bin/env Rscript
# December 12,2020
# Version 1.4.0
#
#' @title FIREcaller: an R package for detecting frequently interacting regions from Hi-C data
#' @description Seperate function for calling differential FIREs outside of the main FIREcaller function
#' @param data a data frame of normalized cis-interactions, which are computed by FIREcaller, for each sample/replicate.
#' @param targetFile a vector specifying the class of each sample/replicate.
#' @param coef specifies the column number or column name for the interested coefficient or contrast of the linear model.
#' @param adjust.method specifies the method used to adjust the p-values for multiple testing. Possible options include "none", "BH", "BY" and "holm". Default is \code{adjust.method = "BH"}.
#' @param adj.pvalue sets the cutoff value for adjusted p-values where only the interactions with lower adjusted p-values will be returned in the output.
#' @param lfc sets the minimum absolute log2-fold-change required. Only the interactions with larger absolute log2-fold-change value are returned in the output.
#' @import limma
#' @export
InteractionCompare <- function(data, targetFile, coef = 2, adjust.method = "BH", adj.pvalue = 0.05, lfc = log2(2)){
  designMatrix <- model.matrix(~targetFile)
  fit <- lmFit(data, designMatrix)
  efit <- eBayes(fit)
  test_result <- topTable(efit, num=inf, coef=coef, adjust.method = adjust.metho, p.value = adj.pvalue, sort.by="P", lfc=lfc)
  write.table(test_result, file = "Differentially_Interacting_FIREs.txt", sep = "\t", quote = F)   
  return(test_result)
}
  








