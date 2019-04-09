# FIREcaller
FIREcaller: an R package for detecting frequently interacting re-gions from Hi-C data

FIREcaller is maintained by Cheynna Crowley [cacrowle@live.unc.edu] and Yuchen Yang [yyuchen@email.unc.edu]

## News and Updates
April 08, 2018
* Version 0.99.0 released
  + First offical release
  + It can work on Window, Mac and Linux platforms

## Installation
You can install FIREcaller from github with:
```{r install}
install.packages("devtools")

devtools::install_github("yycunc/FIREcaller")
```

## FIREcaller Examples

In this tutorial, we will analyze Hippocampus dataset from Schmitt *et al*., (Cell Reports, 2016), which contains Hi-C contact matrix of 22 autosomes from human Hippocampus. Hi-C input files and the mappability file can be downloaded from [Yun Li Group website](https://yunliweb.its.unc.edu/FIREcaller/download.php).
  

### Setup the library
```{r init}
library("FIREcaller")
```

### FIRE and SuperFIRE calling
#### Set working directory

The directory should contain Hi-C input files of all autosomes and the mappability file required for the analysis.

```{r setup for working directory}
setwd('~/Documents/Schmitt_Hippo_40KB_input/')
```

#### Hi-C input file

The Hi-C input file *prefix.list* is defined according to the naming convention of the NxN matrices. The HiC Input files need to be sample and chromosome-specific NxN contact frequency matrix with no row name and no column name.The files for the contact matrices must have the naming convention "\${prefix}_chr\${number}.gz"

```{r define the prefix.list according to the naming convention of the NxN matrices, warning=FALSE}
prefix.list <- c('Hippo')
```

Here is an example for the required format of Hi-C input files.

```{r an example for the required format of Hi-C input files, warning=FALSE}
Hippo_chr1 <- read.table("Hippo_chr1.gz", header = FALSE)

# A subset contact matrix from 100 ~ 110th rows and 100 ~ 110th columns
Hippo_chr1[100:110,100:110]
```

#### Define the genome build

Users are required to define the genome build type. It can be "hg19" or "GRCh38" for human data, and "mm9" or "mm10" for mouse data. If missing, an error message is returned.

```{r define the genome build, message=FALSE}
gb<-'hg19'
```

#### Define the name of the mappability file

There are some mappability files of different genome build (hg19, GRCh38, mm9 and mm10) and different resolutions (10 kb and 40 kb) available in [YunJiang's website](http://enhancer.sdsc.edu/yunjiang/resources/genomic_features/).

The mappability file needs to be in the format of column names as =c('chr','start', 'end', 'F', 'GC','M'), and the chromosome column needs to be in the format 'chr${number}'.

```{r define the name of the mappability file, message=FALSE}
map_file<-'F_GC_M_HindIII_40KB_hg19.txt.gz'
```

Here is an example for the required format of the mappability file

```{r an example for the required format of the mappability file, message=FALSE}
# The format of mappability file
F_GC_M_HindIII_40KB_hg19 = read.table("F_GC_M_HindIII_40KB_hg19.txt.gz", header = TRUE)
head(F_GC_M_HindIII_40KB_hg19)
```

Users can also use their own mappability file in the same format.

#### Define whether to remove MHC region

The default setting is "TRUE"", that is, to remove the MHC region 

```{r define whether to remove MHC region, message=FALSE}
rm_mhc <- TRUE
```

#### Call FIREs and SuperFIREs

Using FIREcaller function, we call both FIREs and SuperFIREs for 22 autosomes of Hippocampus dataset.
  
```{r call FIRE and SuperFIRE for 22 autosomes of Hippocampus dataset}
FIREcaller(prefix.list, gb, map_file, rm_mhc)
```

In this case, two sets of files will be returned: one for FIREs and the other one for SuperFIREs.

```{r An example for Fire and SuperFire outputs}
# An example for FIRE output
FIRE_output = read.table("FIRE_ANALYSIS_40KB.txt", header = T)
head(FIRE_output)

# An example for SuperFIRE output
SuperFIRE_output = read.table("super_FIRE_call_Hippo.txt", header = T)
head(SuperFIRE_output)
```
