# !/usr/bin/env Rscript
# March 10,2021
# Version 1.4.1
#
#' @title FIREcaller: an R package for detecting frequently interacting regions from Hi-C data
#' @description This function FIREcaller() in the FIREcaller package ( user-friendly R package for detecting FIREs from Hi-C data), 
#' For default parameters: FIREcaller takes raw Hi-C NxN contact matrix as input, performs within-sample and cross-sample normalization via
#' HiCNormCis and quantile normalization respectively, and outputs FIRE scores, FIREs and super-FIREs.
#' Input is either an NxN contact matrices in the form of a .gz file, or data in .cool and .hic format.
#' @usage FIREcaller(file.list=(...), gb=c("hg19","GRCh38","mm9","mm10"), map_file="", binsize=c(10000,20000,40000),upper_cis=200000, normalized=c("TRUE","FALSE"),filter=c("TRUE","FALSE"),rm_mhc=c("TRUE","FALSE"),rm_EBL=c("TRUE","FALSE"), rm_perc=0.25, dist=c('poisson','nb'),alpha=0.05,diff_fires=c('TRUE','FALSE')) 
#' @param file.list a list of files used for FIREcaller. If in .gz format, the naming convention is ${prefix}.chr${chr#}.gz. If in .cool format, the naming convention is ${prefix}.cool.
#' @param gb a string that defines the genome build type. If missing, an error message is returned.
#' @param map_file a string that defines the name of the mappability file specific to the samples genome build, restriction enzyme, and resolution. Only contains chromosomes you want to input. See read me for format.
#' @param binsize a numeric value for the binsize. Default is 40000 (40Kb) with other options being 10Kb or 20Kb.
#' @param upper_cis a bound for the cis-interactions calculation. The default is 200000 (200Kb).
#' @param normalized  a logical value for whether the input matrices are ALREADY normalized. If TRUE, the normalization procedures are skipped. Default=FALSE.
#' @param rm_perc is the percentage of "bad-bins" in a cis-interaction calculation to filter. Default is 0.25 (25\% filtered)
#' @param rm_mhc a logical value indicating whether to remove the MHC region of the sample. Default is "TRUE".
#' @param rm_EBL a logical value indicating whether to remove the ENCODE blacklist regions of the sample. Default is "TRUE".
#' @param dist is the distribution specification for the HiCNormCis normalization and FIREscore calculation. The default is Poisson.
#' @param alpha is the type 1 error for the p-value cut off. Default is 0.05.
#' @param diff_fires a logical value for whether to include the differential FIRE analysis. Samples need to have {_rep1} or {_rep2} differences in the file name.Default=FALSE.
#' @details The process includes calculating the raw fire scores, filtering (with the option of removing the MHC region and ENCODE black list), HiCNormCis, Quantile Normalization (if number of samples > 1), highlighting significant Fire Scores, and calculating the SuperFires.
#' @return Two sets of files will be returned at default. The total number of files outputted are 1+ (number of prefixes/samples):
#' @return \itemize{
#'     \item Fire: a single text file is outputted for all the samples and all chromosomes.This file contains the Fire Score, associated ln(pvalue), and an indicator if the region is a FIRE or not with I(pvalues > -ln(0.05)).
#'     \item SuperFire: a text file for each sample with a list of Super Fires and corresponding -log10(pvalue).
#'     }
#' @note   
#' Mappability files are available https://yunliweb.its.unc.edu/FIREcaller/
#' @seealso Paper: https://doi.org/10.1016/j.celrep.2016.10.061  ; https://doi.org/10.1101/619288
#' @author Crowley, Cheynna Anne <cacrowle@live.unc.edu>, Yuchen Yang <yyuchen@email.unc.edu>,
#' Ming Hu <afhuming@gmail.com>, Yun Li <yunli@med.unc.edu>
#' @references Cheynna Crowley, Yuchen Yang, Ming Hu, Yun Li. FIREcaller: an R package for detecting frequently interacting regions from Hi-C data
#' @import MASS
#' @import preprocessCore
#' @import data.table
#' @import stringr
#' @import circlize
#' @import limma
#' @import HiCcompare
#' @importFrom stats glm
#' @importFrom stats pnorm
#' @importFrom stats sd
#' @importFrom utils write.table
#' @importFrom utils read.table
#' @examples
#' # set working directory: the location of the NxN matrices and mappability files
#' setwd('~/Documents/Schmitt_Hippo_40KB_input/')
#'
#' # define the filename following if in .gz format, the naming convention is ${prefix}.chr${chr#}.gz. If in .cool format, the naming convention is ${prefix}.cool
#' file.list <- c(paste0('Hippo_chr',1:22,'.gz'))
#'
#' # define the genome build
#' gb<-'hg19'
#' 
#' # define the name of the mappability file
#' map_file<-'Hind3_hg19_40Kb_encodeBL_F_GC_M_auto.txt.gz'
#'
#' # define the binsize. Default=40000 (40Kb). Other recommended bin sizes are 10000 (10Kb) and 20000 (20Kb).
#' binsize<-40000
#' 
#' # define the upper bound of the cis-interactions; default=200000 (200Kb); if not a multiple of the bin, then takes the ceiling;
#' upper_cis<-200000
#' 
#' # define if the input matrix is ALREADY normalized. Default=FALSE. If true, it will skip within-normalization step.
#' normalized<-FALSE
#' 
#' # define whether to remove MHC region; Default=TRUE
#' rm_mhc <- TRUE
#' 
#' # define whether to remove ENCODE blacklist region; Default=TRUE
#' rm_EBL<- TRUE
#'
#' # define the percentage to problematic bins allowed in the cis-interaction calculation (0-1); Default is 25%.
#' rm_perc<-0.25
#' 
#' # define whether the distribution should be poisson or negative binomial; Default=Poisson.
#' dist<-'poisson'
#'
#' # define the alpha cut off for a significant p-value. Default=0.05.
#' alpha<-0.05
#' 
#' # define if a circos plot should be created of FIREs and super-FIREs; Default=FALSE
#' plots<-FALSE
#' 
#' # specify if differential fires should be calculated between 2 samples and atleast 2 replicates per sample; Defaul=FALSE
#' diff_fires<-FALSE
#'
#' # run the function
#' FIREcaller(file.list, gb, map_file,binsize=40000, upper_cis=200000,normalized=FALSE, rm_mhc = TRUE,rm_EBL=TRUE, rm_perc=0.25, dist='poisson',alpha=0.05, plots=FALSE,diff_fires=FALSE)
#' 
#' @export

FIREcaller <- function(file.list,gb, map_file,binsize=40000, upper_cis=200000,normalized=FALSE, rm_mhc = TRUE,rm_EBL=TRUE, rm_perc=0.25, dist='poisson',alpha=0.05, plots=FALSE,diff_fires=FALSE){

  #reference files
  options(scipen = 999)
  ref_file <-as.data.frame(fread(map_file))
  colnames(ref_file) <- c("chr", "start", "end", "F", "GC", "M","EBL")
  res <- as.numeric(binsize)
  
  bin_num <- ceiling(upper_cis/res)
  output_name<-paste0('FIRE_ANALYSIS_', res,"_",upper_cis,'_',dist,'.txt')
  
  #check file type;
  temp<-file.list[1]
  
  if(str_count(temp,".gz")==1){
    prefix.list<-unique(str_split_n(file.list,"_chr",1))
  }
  
  if(str_count(temp,".cool")==1){
    prefix.list<-unique(str_split_n(file.list,".cool",1))
  }
  
  if(str_count(temp,".hic")==1){
    prefix.list<-unique(str_split_n(file.list,".hic",1))
  }

  #(1) calculate cis-interactions 
  #(1A) .gz -- 
  if(str_count(temp,".gz")==1){
    chr_list<-paste0("chr",unique(str_split_n(str_split_n(file.list,"_chr",2),"[.]",1)))
    t <- cis_15KB_200KB(prefix.list,chr_list, bin_num, ref_file)
  }
  
  #(1B) cooler
  if(str_count(temp,".cool")==1){
    t <- cis_all_cool_sample(gb,binsize,bin_num,file.list,ref_file,prefix.list)
  }
  
  #(1C) hic
  if(str_count(temp,".hic")==1){
    t <- cis_all_hic_sample(gb,binsize,bin_num,file.list,ref_file,prefix.list)
  }

  #(2) filtering
    t2 <- filter_count(file = t, rm_mhc, bin_num, gb, res,rm_perc,rm_EBL)

  ###########if normalized data;
  
  if(normalized==FALSE){
    #within sample normalization
    if(dist=='poisson'){t3 <- HiCNormCis(file = t2)} else if (dist=='nb'){t3<-HiCNormCis.2(file=t2)}} else{t3<-t2}
    
    
    #across sample normalization
    if(length(prefix.list) > 1){t4 <- quantile_norm(file = t3)}else{t4 <- t3}
  
    
  
  #final fire call
  final_fire <- Fire_Call(file = t4, prefix.list, alpha, output_name)
  
  #super fire call and plots
  Super_Fire(final_fire, prefix.list,plots,gb)
  
  #differential fires;
  
  if(diff_fires==TRUE){
    diff_interactions(final_fire)
  }
}


#####################################################################################
# Functions for string split 
#####################################################################################

subset_safely <- function(x, index) {
  if (length(x) < index) {
    return(NA_character_)
  }
  x[[index]]
}

str_split_n <- function(string, pattern, n){
  out <- str_split(string, pattern)
  vapply(out, subset_safely, character(1L), index = n)
}


#####################################################################################
# Functions for FIRE caller with NxN raw and normalized contact matrixto prepare for calc_cis()
#####################################################################################

cis_15KB_200KB <- function(prefix.list, chr_list, bin_num, ref_file){
  all <- NULL
  for(j in 1:length(chr_list)){
    chr <- chr_list[j]
    x <- cis_single_chr(prefix.list, chr, bin_num, ref_file)
    all <- rbind(all,x)
  }
  return(all)
}

cis_single_chr <- function(prefix.list, chr, bin_num,ref_file){
  file_list <- paste0(prefix.list, '_', chr, '.gz')
  ref_file_chr <- ref_file[ref_file$chr == chr,]
  all_samples <- as.data.frame(ref_file_chr)

  for(i in 1:length(file_list)){
    matrix <- as.matrix(fread(file_list[i], data.table = FALSE))
    colnames(matrix)<-rownames(matrix)
    sample <- prefix.list[i]
    diag(matrix) <- 0
    x <- calculate_cis(matrix, bin_num, sample, ref_file_chr)
    all_samples <- merge(all_samples, x, by = c("chr", "start", "end", "F", "GC", "M","EBL"), sort = FALSE, all = TRUE)
  }
  
  return(all_samples)
}


#####################################################################################
# Functions for FIRE caller with .hic samples to prepare for calc_cis()
#####################################################################################
cis_all_hic_sample<-function(gb,binsize,bin_num,file.list,ref_file,prefix.list){
  t2<-ref_file
  for(i in 1:length(file.list)){
    file<-file.list[i]
    t<-cis_hic_sample(file,binsize, gb,ref_file,bin_num)
  }
  t2<-cbind(t2,t)
  colnames(t2)<-c("chr","start","end","F","GC","M","EBL",prefix.list)
  return(t2)
}

cis_hic_sample<-function(file,binsize, gb,ref_file,bin_num){
  load(system.file("extdata","chrom_sizes.rda",package="FIREcaller"))
  
  if(gb %in% c('mm9','mm10')){
    chroms <- paste0("chr", c(1:19, "X"))
  }
  if(gb %in% c('hg19','gb38')){
    chroms <- paste0("chr", c(1:22, "X"))
  }
  
  m3<-NULL
  
  name<-str_split_n(file,".hic",1)
  input<-file
  
  for(chrom in chroms){
    name.out<-paste0(name,"_",chrom,".txt")
    #juicer_dir<-system.file("extdata","juicer_tools_1.21.01.jar",package="FIREcaller")
    juicer <- "java -jar juicer_tools_1.21.01.jar"
    juicer_command <- paste(juicer, "dump observed NONE", input, chrom, chrom, "BP", binsize, name.out)
    print(juicer_command)
    system(juicer_command)
    data<-fread(paste0(name,"_",chrom,".txt"))
    colnames(data)<-c('start','end','N')
    long_matrix<-as.data.frame(data)
    
    all_chr<-as.data.frame(all_chr)
    chr_size<-all_chr[all_chr$chr==chrom,grepl(gb, colnames(all_chr))]
    
    
    matrix_size<-ceiling(chr_size/binsize)
    colnames(long_matrix)<-c('row','column','N')
    
    long_matrix$row<-ceiling(long_matrix$row/binsize)
    long_matrix$column<-ceiling(long_matrix$column/binsize)
    
    short_matrix <- matrix(0, matrix_size,matrix_size)
    
    for(i in 1:nrow(long_matrix)){
      short_matrix[long_matrix$row[i],long_matrix$column[i]]=long_matrix$N[i] 
      #short_matrix[long_matrix$column[i],long_matrix$row[i]]=long_matrix$N[i] 
    }
    
    ref_file_chr <- ref_file[ref_file$chr == chrom,]
    
    m2<-calculate_cis(short_matrix, bin_num, name, ref_file_chr)
    m3<-rbind(m3,m2)
    
  }
  return(m3)
}



#####################################################################################
# Functions for FIRE caller with .cool samples to prepare for calculate cis_function()
#####################################################################################

cis_all_cool_sample<-function(gb,binsize,bin_num,file.list,ref_file,prefix.list){
  t2<-ref_file
  for(i in 1:length(file.list)){
    sample<-prefix.list[i]
    dat<-cooler2bedpe(path = paste0(file.list[i]))
    t<-cis_inter_sample(dat,gb,bin_num,ref_file,sample,binsize)
    t2<-t2[t2$chr %in% unique(t$chr),]
    t2<-merge(t2,t,by=c("chr","start","end","F","GC","M","EBL"),all=TRUE)
  }
  colnames(t2)<-c("chr","start","end","F","GC","M","EBL",prefix.list)
  return(t2)
}


cis_inter_sample<-function(dat,gb,bin_num,ref_file,sample,binsize){
  load(system.file("extdata","chrom_sizes.rda",package="FIREcaller"))
  all_chr<-as.data.frame(all_chr)
  chroms.0 <- names(dat[[1]])
  chroms<-chroms.0[chroms.0 %in% c(paste0('chr',1:22),'chrX')]
  m3<-NULL
  
  for(chr in chroms){
    ref_file_chr <- ref_file[ref_file$chr == chr,]
    data_mat <- dat[[1]][[chr]] #separate out matrix from main file;
    
    chr_size<-all_chr[all_chr$chr==chr,grepl(toString(gb), colnames(all_chr))]
    matrix_size<-ceiling(chr_size/binsize)
    m1<-cooler_2_matrix(data_mat,matrix_size,binsize) #create a nxn contact matrix
    m1<-as.matrix(m1)
    m2<-calculate_cis(m1, bin_num, sample, ref_file_chr)
    m3<-rbind(m3,m2)
  }
  
  return(m3)
  
}


cooler_2_matrix<-function(data_mat,matrix_size,binsize){
  long_matrix<-data_mat[,c(1,3,6,7)]
  colnames(long_matrix)<-c('chr','row','column','N')
  
  long_matrix$row<-ceiling(long_matrix$row/binsize)
  long_matrix$column<-ceiling(long_matrix$column/binsize)
  
  short_matrix <- matrix(0, matrix_size,matrix_size)
  
  for(i in 1:nrow(long_matrix)){
    short_matrix[long_matrix$row[i],long_matrix$column[i]]=long_matrix$N[i] 
    #short_matrix[long_matrix$column[i],long_matrix$row[i]]=long_matrix$N[i] 
  }
  return(short_matrix)
}


#####################################################################################
# Functions calculate_cis()
#####################################################################################

calculate_cis <- function(matrix, bin_num, sample, ref_file_chr){
  length <- nrow(matrix)
  x <- vector("integer", length = length)
  for(i in 1:length){
    x[i] <- ifelse(i <= bin_num,
                   sum(matrix[i,(i:(i+bin_num))], matrix[(1:i),i]),
                   ifelse((i+bin_num) > length,
                          sum(matrix[i,(i:length)], matrix[(i-bin_num):i,i]), #case 2
                          sum(matrix[((i-bin_num)):i, i], matrix[i,(i:(i+bin_num))]))) #case 3
  }
  x <- as.data.frame(x)
  colnames(x) <- sample
  final <- cbind(ref_file_chr, x)
  return(final)
}


#####################################################################################
# Filtering
#####################################################################################

filter_count <- function(file, rm_mhc, bin_num, gb, res,rm_perc,rm_EBL){
  options(scipen = 999)

  y <- file
  # filter 1: find bad bins-- F=0, GC=0 or M=0
  y$flag <- 0
  y$flag <- ifelse(y$F == 0 | y$GC == 0 | y$M == 0, 1, 0)

  # filter 2: find neighbor bins
  y2 <- y
  y2$bad_neig <- 0
  for(i in 1:nrow(y2)){
    if((i >= bin_num) & (i + bin_num <= nrow(y2))){
      y2$bad_neig[i] <- sum(y2[c(((i - bin_num):i),(i:(i + bin_num))),"flag"])
    }
    if(i < bin_num){
      y2$bad_neig[i] <- sum(y2[c((1:i),(i:(i + bin_num))),"flag"])
    }
    if((i + bin_num) > nrow(y2)){
      y2$bad_neig[i]<-sum(y2[c(((i - bin_num):i),(i:nrow(y2))),"flag"])
    }
  }

  #find percentage of bad bins of
  y2$perc <- y2$bad_neig/(2 * bin_num)

  #remove the bad bins and if 25 percent of the neighbooring bins are "bad"
  y3 <- y2[y2$flag == 0,]
  y4 <- y3[y3$perc <= rm_perc,]

  # filter 3: remove if Mapp<0.9
  y5 <- y4[y4$M > 0.9,]

  #remove the last two columns
  n1 <- ncol(y5)
  n2 <- n1-2
  y6 <- y5[,-c(n2:n1)]

  #filter 4: MHC regions
  if(rm_mhc == TRUE){
    y7 <- remove_mhc(y6, gb, res)
  } else{
    y7 <- y6
  }
  
  if(rm_EBL == TRUE){
    y8 <-y7[y7$EBL==0,]
  } else{
    y8<-y7
  }
  
  
  return(y8)
}

remove_mhc <- function(y, gb, res){
  if(toupper(gb) == 'HG19'){y2 <- mhc_hg19(y, res)}
  if(toupper(gb) == 'GRCH38'){y2 <- mhc_grch38(y, res)}
  if(toupper(gb) == 'MM9'){y2 <- mhc_mm9(y, res)}
  if(toupper(gb) == 'MM10'){y2 <- mhc_mm10(y, res)}
  return(y2)
}

mhc_hg19 <- function(y, res) {
  y_remove <- (y$chr == "chr6" & y$start > floor(28477797/res)*res & y$end <= ceiling(33448354/res)*res)
  y <- y[!y_remove,]
  return(y)
}

mhc_grch38 <- function(y, res) {
  y_remove <- (y$chr == "chr6" & y$start > floor(28510120/res)*res & y$end <= ceiling(33480577/res)*res)
  y <- y[!y_remove,]
  return(y)
}

mhc_mm9 <- function(y, res) {
  y_remove <- (y$chr == "chr17" & y$start > floor(33888191/res)*res & y$end <= ceiling(35744546/res)*res)
  y_remove2 <- (y$chr == "chr17" & y$start > floor(36230820/res)*res & y$end <= ceiling(38050373/res)*res)
  y <- y[!y_remove,]
  y <- y[!y_remove2,]
  return(y)
}

mhc_mm10 <- function(y,res) {
  y_remove <- (y$chr == "chr17" & y$start > floor(33681276/res)*res & y$end <= ceiling(38548659/res)*res)
  y <- y[!y_remove,]
  return(y)
}


#####################################################################################
# Normalization
#####################################################################################

HiCNormCis <- function(file){
  x <- file[,-7]
  corout <- matrix(0, nrow = ncol(x)-6, ncol = 6)
  cn<-colnames(x[,-c(1:6)])
  FIREscore <- x[,1:3]
  for(id in 7:(ncol(x))){
    y <- x[,id]
    fit <- glm(x[,id] ~ x$F + x$GC + x$M, family="poisson")
    #summary(fit)
    coeff <- round(fit$coeff, 8)
    res <- round(x[,id]/exp(coeff[1] + coeff[2]*x$F + coeff[3]*x$GC + coeff[4]*x$M), 4)
    FIREscore <- cbind(FIREscore, res)
  }

  FIREscore <- as.data.frame(FIREscore)
  colnames(FIREscore)<-c('chr','start','end',cn)
  return(FIREscore)
}

HiCNormCis.2 <- function(file){
  x <- file[,-7]
  corout <- matrix(0, nrow = ncol(x)-6, ncol = 6)
  cn<-colnames(x[,-c(1:6)])
  FIREscore <- x[,1:3]
    for(id in 7:(ncol(x))){
      y <- x[,id]
      fit<-glm.nb(x[,id] ~ x$F + x$GC + x$M)
      #summary(fit)
      coeff <- round(fit$coeff, 8)
      res <- round(x[,id]/exp(coeff[1] + coeff[2]*x$F + coeff[3]*x$GC + coeff[4]*x$M), 4)
      FIREscore <- cbind(FIREscore, res)
    }

  FIREscore <- as.data.frame(FIREscore)
  colnames(FIREscore)<-c('chr','start','end',cn)
  return(FIREscore)
}

quantile_norm <- function(file){
  x <- file
  y <- x[, 4:ncol(x)]
  yqq <- normalize.quantiles(as.matrix(y))
  z <- cbind(x[,1:3], round(yqq,4))
  z <- as.data.frame(z)
  colnames(z) <- colnames(x)
  return(z)
}


#####################################################################################
# FIRE call 
#####################################################################################

Fire_Call <- function(file, prefix.list,alpha,output_name){
  x<- file
  a<-as.numeric(alpha)
  annp <- x[, 1:3] #fire scores
  annf <- x[, 1:3]
  for(id in 4:ncol(x)){
    y <- x[, id]
    mean(y)
    ym  <- mean(y)
    ysd <- sd(y)
    p <- round(-pnorm(y,  mean = ym, sd = ysd, lower.tail = FALSE, log.p = TRUE), 4) #returns log pvalue
    q <- x[ p > -log(a),]
    f <- as.numeric(I(p > -log(a)))
    annp <- cbind(annp, p)
    annf <- cbind(annf, f)
  }

  annp <- as.data.frame(annp) # log(p-values) scores
  annf <- as.data.frame(annf) #indicators for pvalues > -log(0.05)
  fires <- x

  colnames(annp) <- c('chr', 'start', 'end', paste0(prefix.list, '_neg_ln_pval'))
  colnames(annf) <- c('chr', 'start', 'end', paste0(prefix.list, '_indicator'))
  colnames(fires) <- c('chr', 'start', 'end', paste0(prefix.list, '_norm_cis'))

  final0 <- merge(fires, annp, by = c('chr', 'start', 'end'), all=TRUE, sort = FALSE)
  final <- merge(final0, annf, by = c('chr', 'start', 'end'), all = TRUE, sort = FALSE)
  write.table(final, output_name, quote = FALSE, row.names = FALSE)
  
  return(final)
}

#####################################################################################
# Super-FIRE
#####################################################################################

Super_Fire <- function(final_fire, prefix.list,plots,gb){
  final_fire<-as.data.frame(final_fire)
  length_pl <- as.numeric(length(prefix.list))
  NP <- final_fire[, c(1:3,(3+1+length_pl):(3+2*length_pl))]
  colnames(NP) <- c('chr', 'start', 'end', prefix.list)
  ID <- final_fire[, c(1:3,(1+3+2*length_pl):ncol(final_fire))]
  colnames(ID) <- c('chr', 'start', 'end', prefix.list)
  
  
  for(INDEX in 4:(length_pl+3)){
    x0 <- ID[, c(1:3, INDEX) ]
    y0 <- NP[, c(1:3, INDEX) ]
    
    x <- x0[x0[,4]==1,]
    z <- y0[x0[,4]==1,]
    
    final <- NULL
    chr.list<-unique(x0$chr)
    for(chrid in 1:length(chr.list)){
      u <- z[ z[,1] == chr.list[chrid],]
      u2<-as.data.frame(u[1,])
      if(nrow(u)>2){
        out <- NULL
        
        for(i in 2: (nrow(u)) ){
          t<-nrow(u2)
          start <- as.numeric(u2[t, 2])
          end <- as.numeric(u2[t, 3])
          sum <-as.numeric(u2[t,4])
          if(abs(end-as.numeric(as.numeric(u[i, 2])))==0){
            u2[t,1]<-chr.list[chrid]
            u2[t,2]<-start
            u2[t,3]<-as.numeric(as.numeric(u[i, 3]))
            u2[t,4]<-sum+as.numeric(as.numeric(u[i, 4]))
          }
          if(abs(end-as.numeric(as.numeric(u[i, 2])))>0){
            u2[t+1,1]<-chr.list[chrid]
            u2[t+1,2]<-as.numeric(u[i, 2])
            u2[t+1,3]<-as.numeric(u[i, 3])
            u2[t+1,4]<-as.numeric(u[i, 4])
          }
        }
        
        colnames(u2)<-c('chr', 'start', 'end', 'cum_FIRE_score')
        
        final <- rbind(final, u2)
      }
    }
    x <- final
    
    #max(final$cum_FIRE_score)
    y <- x[order(x$cum_FIRE_score),]
    rank<-seq(1,nrow(y),1)
    z <- cbind(rank, y)
    #scale data to 0,1
    z[,6] <- z[,1]/max(z[,1]) #divide each rank by the max
    z[,7] <- z[,5]/max(z[,5]) #divide each score by the max
    
    z[,8] <-  1/sqrt(2)*z[,6] + 1/sqrt(2)*z[,7]
    z[,9] <- -1/sqrt(2)*z[,6] + 1/sqrt(2)*z[,7]
    
    RefPoint <- z[z[,9] == min(z[, 9]), 1] # 1423
    RefValue <- z[RefPoint, 5]
    
    xout <- z[z[,5] >= RefValue,c(2:5)]
    #xout$NegLog10Pvalue <- round(xout$cum_FIRE_score/log(10), 4)
    if(plots==TRUE){
      plot_fire_sf(FIRE=x0,SF=xout,name=colnames(ID)[INDEX],gb)
    }
    
    write.table(xout, file = paste0('super_FIRE_call_', colnames(ID)[INDEX], '.txt'), row.names = F, col.names = T, sep = '\t', quote = F)
  }
}

###################################################################################
# Plot FIREs and SuperFIREs
###################################################################################

plot_fire_sf<-function(FIRE,SF,name,gb){
  list_chr<- paste0('chr',order(unique(FIRE$chr)))
  colnames(FIRE)<-c('chr','start','end','I')
  FIRE.2<-FIRE[FIRE$I==1,c(1:3)]
  SF.2<-SF[,-4]
  circos.clear()
  png(paste0(name,"FIREcaller_circos.png"))
  circos.initializeWithIdeogram(species=gb,chromosome.index = list_chr)
  circos.genomicDensity(FIRE.2, col = 'red', track.height = 0.2)
  circos.genomicDensity(SF.2, col = 'darkblue', track.height = 0.2)
  dev.off()

}



###################################################################################
# Differential FIRE analysis
###################################################################################

diff_interactions<-function(FC_results){
  FC.1<-FC_results
  n<-(ncol(FC.1)-3)/6 #number of samples

  rownames(FC.1) = paste0(FC.1$chr, ":", FC.1$start, "-", FC.1$end) #unique for merging back
  FC.2<-FC.1[ , grepl( "_norm_cis" , names(FC.1) ) ] #keep only the normalized cis-interactions
  FC.3<-FC.2[,order(colnames(FC.2))] #order by alphabetical so different samples are split

  # con-case, because in terms of spelling con is later than case
  targetFile<-c(rep("sample1", n), rep("sample2", n)) # sample 1 is the first in the alphabet
  
  designMatrix<-model.matrix(~targetFile) 
  fit <- lmFit(FC.3, designMatrix)
  efit <- eBayes(fit)
  
  # output the full results
  test_results=topTable(efit, num=Inf, coef=2, adjust.method = "BH", sort.by = "P", lfc=log2(2))
  
  #merge back with the FIREcaller results by row.names
  test_result.2<-merge(FC.1,test_results,by=0, all.y=TRUE) 
  
  #remove row.names and only keep rows with atleast one indicator=1
  test_result.3<- test_result.2[rowSums(test_result.2[,grepl( "_indicator" , names(test_result.2))])>0,-1]

  write.table(test_result.3, file = "Differentially_Interacting_FIREs.txt", sep = "\t", quote = F)                                                                                                       
  }





