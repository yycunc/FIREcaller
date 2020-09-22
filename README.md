# FIREcaller
FIREcaller: an R package for detecting frequently interacting regions from Hi-C data

FIREcaller is maintained by Cheynna Crowley [cacrowle@live.unc.edu], Yuchen Yang [yyuchen@email.unc.edu] and Yun Li [yun_li@med.unc.edu].

## News and Updates
Sep 20, 2020
* Version 1.30 released
  + Allows user to change the cis-interacting region (200Kb by default)
  + Allows user to pick the regression distribution (Poisson or Negative Binomial; Poisson by default)
  + Allows the user to select the percentage of problematic bins filtered (25% by default)
  + Added filtering of the ENCODE blacklist regions
  + Can do single chromosomes, or full list based on the mappability file.

June 30, 2020
* Version 1.2.0 released
  + Allows the option to include chrX
  + Allows symmetric and upper triangular matrices as input
  
Feb 17, 2020
* Python version released by Jakub Lipiński
  + https://github.com/jakublipinski/python-FIREcaller

Oct 21, 2019
* Version 1.10 released
  + Loosened the criterion for a "bad bin". Previously, a bin is considered to be bad if ANY of its neighboring bins within 200kb is bad. Now loosened to allow <=25% bad neighboring bins.
  + Fixed an error that caused a lack of convergence.

Sep 4, 2019
* Version 1.00 released
  + The Default value of "rm_mhc" is set to "TRUE" instead of "NULL".
  + Fixing an error in SuperFIRE calling.

April 08, 2019
* Version 0.99.0 released
  + First offical release.
  + It can work on Window, Mac and Linux platforms.

## Installation
Users can install FIREcaller from github with:
```{r install}
install.packages("devtools")

devtools::install_github("yycunc/FIREcaller")
```

Or downloaded the package from https://yunliweb.its.unc.edu/FIREcaller/

```{r}
install.packages("~/FIREcaller_1.30.tar.gz",repos=NULL, type="source")

```

Users can also install a python version of FIREcaller from github page
https://github.com/jakublipinski/python-FIREcaller.

## FIREcaller Examples

In this tutorial, we will analyze Hippocampus dataset from Schmitt *et al*., ([Cell Reports, 2016](https://www.cell.com/cell-reports/pdfExtended/S2211-1247(16)31481-4)). It contains 22 Hi-C NxN contact matrix of the autosomal chromosomes. Hi-C input files and the mappability file can be downloaded from [Yun Li Group website](https://yunliweb.its.unc.edu/FIREcaller/download.php).
  

### Setup the library
```{r init}
library("FIREcaller")
```

### FIRE and SuperFIRE calling
#### Set working directory

The directory should contain Hi-C input files of all autosomes and the mappability file required for the analysis.

```{r setup for working directory}
setwd('~/Desktop/FIREcaller_example')
```

#### Hi-C input file

The Hi-C input file *prefix.list* is defined according to the naming convention of the NxN matrices. The HiC Input files need to be sample and chromosome-specific NxN contact frequency matrix with no row name and no column name.The files for the contact matrices must have the naming convention "\${prefix}_chr\${number}.gz".

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

There are some mappability files of different genome build (hg19, GRCh38, mm9 and mm10) and different resolutions (10 kb and 40 kb) available in [Yun Li's's website](https://yunliweb.its.unc.edu/FIREcaller/download.php).

The mappability file needs to be in the format of column names as =c('chr','start', 'end', 'F', 'GC','M'), and the chromosome column needs to be in the format 'chr${number}'. The chromosomes in the file need to directly relate to the chromosomes in the process. For example, the mappability file will only contain information on chromosomes 1 through 22 to do an analysis on chromosome 1 through 22.

```{r define the name of the mappability file, message=FALSE}
map_file<-'Hind3_hg19_40Kb_encodeBL_F_GC_M_auto.txt.gz'
```

Here is an example for the required format of the mappability file.

```{r an example for the required format of the mappability file, message=FALSE}
# The format of mappability file
Hind3_hg19_40Kb_encodeBL_F_GC_M_auto= read.table("Hind3_hg19_40Kb_encodeBL_F_GC_M_auto.txt.gz", header = TRUE)
head(Hind3_hg19_40Kb_encodeBL_F_GC_M_auto)
```

Here is an example of the chromosomes present in the mappability file.

```{r,warning=FALSE}
unique(Hind3_hg19_40Kb_encodeBL_F_GC_M_auto$chr)
```

Users can also use their own mappability file in the same format.

#### Define whether to remove MHC region

The MHC regions defined in the R package are:<br/>

| GB   | CHR  | BP RANGE  |
| ---- |:---: | ---------:|
|hg19  | chr6 |  28477797-33448354                   |
|GRCh38| chr6 | 28510120-33480577                    |
|mm9   |chr17 |33888191-35744546 & 36230820-38050373 |
|mm10  |chr17 | 33681276 & 38548659                  |

The default setting is "TRUE"", that is, to remove the MHC region.

```{r define whether to remove MHC region, message=FALSE}
rm_mhc <- TRUE
```

#### Define whether to remove the ENCODE black list regions

The ENCODE blacklist regions are described [here](https://www.nature.com/articles/s41598-019-45839-z). The EBL variable in the mappability file is an indicator of whether it is a black list region (1) or not (0).

```{r define whether to remove ENCODE black list regions, message=FALSE}
rm_EBL <- TRUE
```

The default setting is "TRUE" to remove the ENCODE blacklist regions.

#### Change the cis-interacting regions
The user has the option to change the cis-interacting region threshold. Default is 200Kb.

```{r define cisinteracting region, message=FALSE}
upper_cis=200000
```

#### Change the regression distribution
The user has the option to change the regression distribution used. The HiCNormCis method uses Poisson ("poisson") distribution, but the user can change to negative binomial ("nb")

```{r change the regression distribution, message=FALSE}
dist='poisson'
```

#### Change the filtering threshold 
The user has the option to change the filtering threshold, where if a cis-interaction is calculated by more than 25% bins that contains a mappability of 0, a GC content of 0 , or a effective fragment length of 0,then it is also filtered.

```{r define filter, message=FALSE}
rm_perc=0.25
```

#### Call FIREs and SuperFIREs

Using FIREcaller function, we call both FIREs and SuperFIREs for 22 autosomes of Hippocampus dataset.
  
```{r call FIRE and SuperFIRE for 22 autosomes of Hippocampus dataset}
FIREcaller(prefix.list, gb, map_file, rm_mhc, rm_EBL, upper_cis, dist,rm_perc)
```

In this case, two sets of files will be returned: one for FIREs and the other one for SuperFIREs. If multiple (*n*) prefix's are specified in the *prefix.list*, there will be one file for FIREs and *n* SuperFIRE files.

```{r An example for Fire and SuperFire outputs}
# An example for FIRE output
FIRE_output = read.table("FIRE_ANALYSIS_40000_200000_poisson.txt", header = T)
head(FIRE_output)

# An example for SuperFIRE output
SuperFIRE_output = read.table("super_FIRE_call_Hippo.txt", header = T)
head(SuperFIRE_output)
```

## Citation
Crowley, C., Yang, Y., Qiu, Y., Hu, B., Lipiński, J., Plewczynski, D., Won, H., Ren, B., Hu, M., Li, Y. FIREcaller: Detecting Frequently Interacting Regions from Hi-C Data. *bioRxiv*, doi: https://doi.org/10.1101/619288.
