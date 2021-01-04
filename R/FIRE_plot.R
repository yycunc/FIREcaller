# !/usr/bin/env Rscript
# December 12,2020
# Version 1.4.0
#
#' @title FIREcaller: an R package for detecting frequently interacting regions from Hi-C data
#' @description This function FIRE_plot() in the FIREcaller package (a user-friendly R package for detecting FIREs from Hi-C data), 
#' Visualized the FIREcaller FIRE and superFIRE results for a sample and plots them in a circos plot by utilizing the circlize function in R, with the allowance of 3 additional annotation data.
#' @usage FIRE_plot(FIRE="",prefix="", gb=c('hg19','gb38','mm9','mm10'),chr.list=c('chr1','chr2',...),outname="",anno_list=c(""),color=c('red','darkgreen','orange','blue','purple'))
#' 
#' @param FIRE The name of the FIREcaller output for identifying the FIREs. With default parameters, this file is 'FIRE_ANALYSIS_40000_200000_poisson.txt'.
#' @param prefix The prefix of the sample used in the "prefix.list" option in the FIREcaller function. 
#' @param gb The name of the genome build for the ideogram. Default is hg19 but options include gb38,mm9, and mm10
#' @param chr.list The list of chromosomes visualized 
#' @param outname A string of the name for a JPEG file of the circos plot is desired. If a jpeg is not desired, set to NA. 
#' @param anno_list The optional list of files (up to 3) to also include into the circos plot, such as TADs, enhancers, or super-enhancers, with each file containing 3 columns (chr/start/end  position)
#' @param color The optional list of colors for the tracks if the default is not preferred. The number of colors specified needs to be >= the number of annotations 
#' @details The process includes  Quantile Normalization (if number of samples > 1 and qqnorm=TRUE), and highlighting significant Fire Scores, and calculating the Super Fires.
#' @return Two sets of files will be returned. The total number of files outputted are 1+ (number of prefixes/samples):
#' @return \itemize{
#'     \item Fire: a single text file is outputted for all the samples and all chromosomes.This file contains the Fire Score, associated ln(pvalue), and an indicator if the region is a FIRE or not with I(pvalues > -ln(0.05)).
#'     \item SuperFire: a text file for each sample with a list of Super Fires and corresponding -log10(pvalue).
#'     }
#' @note The input for this function should be same format as the output of the complement function find_cis_inter()
#' @seealso Paper: https://doi.org/10.1016/j.celrep.2016.10.061  bioarxiv: https://doi.org/10.1101/619288
#' @author Crowley, Cheynna Anne <cacrowle@live.unc.edu>, Yuchen Yang <yyuchen@email.unc.edu>,
#' Ming Hu <afhuming@gmail.com>, Yun Li <yunli@med.unc.edu>
#' @references Cheynna Crowley, Yuchen Yang, Ming Hu, Yun Li. FIREcaller: an R package for detecting frequently interacting regions from Hi-C data
#' @import circlize
#' @import data.table
#' @examples
#'
#' # the name of the file for the  normalized cis-interactions
#' norm_file<-"Hippo_norm.txt"
#'
#' # define the alpha cut off for a significant p-value
#' alpha<-0.05
#' 
#' #specify if quantile normalization across samples are needed (qqnorm=TRUE). If the normalization procedure does within and across, then set qqnorm=FALSE.
#' qqnorm<-FALSE
#' 
#' #option for user to specify an outname for the FIRE analysis. Defaulted at FIRE_ANALYSIS.txt
#' outname='FIRE_ANALYSIS.txt'
#' 
#' # run the function
#' FIRE_calc(norm_file, alpha=0.05,qqnorm=TRUE,outname='FIRE_ANALYSIS.txt')
#' 
#' @export

FIRE_plot<-function(FIRE="FIRE_ANALYSIS_40000_200000_poisson.txt",prefix, gb='hg19',chr.list,outname='FIREcaller_circos.jpg',anno_list,color=c('red','orange','darkorange','blue','purple')){
  #FIRE Data format
  FIRE_data<-fread(FIRE)
  FIRE_data.1<-FIRE_data[,c("chr","start","end",paste0(prefix,"_indicator"))]
  colnames(FIRE_data.1)<-c('chr','start','end','indicator')
  FIRE_data.2<-FIRE_data.1[FIRE_data.1$chr %in% chr.list & FIRE_data.1$indicator==1,]
  FIRE_data.3<-FIRE_data.2[,-4]
  
  #superFIRE
  SF<-fread(paste0("super_FIRE_call_",prefix,".txt"))
  SF.1<-SF[SF$chr %in% chr.list,c(1:3)]
  
  if(is.na(outname)==FALSE){
    jpeg(outname)
  }
  #run ciros function;
  if(length(anno_list)==0){
    circos.clear()
    circos.initializeWithIdeogram(species=gb,chromosome.index = chr.list)
    circos.genomicDensity(FIRE_data.3, col = color[1], track.height = 0.1)
    circos.genomicDensity(SF, col = color[2], track.height = 0.1)
  }

  if(length(anno_list)==1){
    temp1<-fread(anno_list)
    circos.clear()
    circos.initializeWithIdeogram(species=gb,chromosome.index = chr.list)
    circos.genomicDensity(FIRE_data.3, col = color[1], track.height = 0.1)
    circos.genomicDensity(SF, col = color[2], track.height = 0.1)
    circos.genomicDensity(temp1, col = color[3], track.height = 0.1) 
  }
  
  if(length(anno_list)==2){
    temp1<-fread(anno_list[1])
    temp2<-fread(anno_list[1])
    circos.clear()
    circos.initializeWithIdeogram(species=gb,chromosome.index = chr.list)
    circos.genomicDensity(FIRE_data.3, col = color[1], track.height = 0.1)
    circos.genomicDensity(SF, col = color[2], track.height = 0.1)
    circos.genomicDensity(anno_list, col = color[3], track.height = 0.1) 
    circos.genomicDensity(temp2, col = color[4], track.height = 0.1)  
  }
  
  if(length(anno_list)==3){
    temp1<-fread(anno_list[1])
    temp2<-fread(anno_list[1])
    temp3<-fread(anno_list[1])
    circos.clear()
    circos.initializeWithIdeogram(species=gb,chromosome.index = chr.list)
    circos.genomicDensity(FIRE_data.3, col = color[1], track.height = 0.1)
    circos.genomicDensity(SF, col = color[2], track.height = 0.1)
    circos.genomicDensity(anno_list, col = color[3], track.height = 0.1) 
    circos.genomicDensity(temp2, col = color[4], track.height = 0.1)  
    circos.genomicDensity(temp3, col = color[5], track.height = 0.1) 
  }
       
  if(is.na(outname)==FALSE){
    dev.off()
  }     
}



 




