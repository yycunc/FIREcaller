# FIREcaller
FIREcaller: an R package for detecting frequently interacting regions from Hi-C data

FIREcaller is maintained by Cheynna Crowley [cacrowle@live.unc.edu], Yuchen Yang [yangych68@mail.sysu.edu.cn] and Yun Li [yun_li@med.unc.edu].

## News and Updates
Nov 3, 2022
* Version 1.42
  + Allows user to perform FIREcaller on data of genome build other than hg19, GRCh38, mm9 and mm10, with a file of chromosome size being provided.
  + Allows user to input their own version of juicer_tools for .hic format data.

Nov 21, 2021
* Version 1.40
  + Provide a R script for users to add the ENCODE blacklist region into mappability file

Dec 17, 2020
* Version 1.40 released
  + Allows users to define bin size
  + Allows users to define whehter the input matrices are ALREADY normalized
  + Allows users to change the alpha cutoff for p-value
  + Allows users to do circos plots for the FIREs and super-FIREs
  + Allows users to perform differential FIRE analysis

Sep 20, 2020
* Version 1.30 released
  + Allows user to change the cis-interacting region (200Kb by default)
  + Allows user to pick the regression distribution (Poisson or Negative Binomial; Poisson by default)
  + Allows the user to select the percentage of problematic bins filtered (25% by default)
  + Added filtering of the ENCODE blacklist regions
  + Can do single chromosomes, or full list based on the mappability file.

## Installation
Users can install FIREcaller from github with:
```{r install}
install.packages("devtools")

devtools::install_github("yycunc/FIREcaller")
```

Or downloaded the package from https://yunliweb.its.unc.edu/FIREcaller/

```{r}
install.packages("~/FIREcaller_1.40.tar.gz",repos=NULL, type="source")

```

Users can also install a python version of FIREcaller from github page
https://github.com/jakublipinski/python-FIREcaller.

## FIREcaller Examples

In this tutorial, we will analyze Hippocampus dataset from Schmitt *et al*., ([Cell Reports, 2016](https://www.cell.com/cell-reports/pdfExtended/S2211-1247(16)31481-4)), which contains Hi-C contact matrix of 22 autosomes from Human Hippocampus. Hi-C input files and the mappability file can be downloaded from [Yun Li Group website](https://yunliweb.its.unc.edu/FIREcaller/download.php).

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

As default, The Hi-C input file *file.list* is defined according to the naming convention of the NxN matrices. The HiC Input files need to be sample and chromosome-specific NxN contact frequency matrix with no row name and no column name.The files for the contact matrices must have the naming convention "{prefix}\_chr{number}.gz". For .hic and .cool files, it should be {sample}.cool or {sample}.hic. For .hic file, juicer is required and needs to be downloaded first, and the version of your juicer_tools needs to be provided later as a parameter.

```{r define the prefix.list according to the naming convention of the NxN matrices, warning=FALSE}
file.list <- paste0('Hippo_chr', 1:22, '.gz')
```

Here is an example for the required format of Hi-C input files.

```{r an example for the required format of Hi-C input files, warning=FALSE}
Hippo_chr1 <- read.table("Hippo_chr1.gz", header = FALSE)

# A subset contact matrix from 100 ~ 110th rows and 100 ~ 110th columns
Hippo_chr1[100:110,100:110]
```

#### Define the genome build

Users are required to define the genome build type. It can be "hg19" or "GRCh38" for human data, "mm9" or "mm10" for mouse data. Other genome build type is also allowed, if so, the chromosome size of the tested species must be provided. If missing, an error message is returned.

```{r define the genome build, message=FALSE}
gb<-'hg19'
```

#### Define the name of the mappability file

Some mappability files of different genome build (hg19, GRCh38, mm9 and mm10) and different resolutions (10 kb and 40 kb) are available in [Yun Li's's website](https://yunliweb.its.unc.edu/FIREcaller/download.php).

The mappability file needs to be in the format of column names as =c('chr','start', 'end', 'F', 'GC','M','EBL'), where, for a certain bin, they correspond to chromosome id, start and end positions, whether filter this bin, GC content and the mappability, whether this bin is overlapped with ENCODE blacklist region or not. The chromosome column needs to be in the format 'chr${number}'. The chromosomes in the file need to directly relate to the chromosomes in the process. For example, the mappability file will only contain information on chromosomes 1 through 22 to do an analysis on chromosome 1 through 22.

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

If you have a six-column mappability file without the EBL column, you can use the script `mappability.R` with blacklist region data to generate a mappability file in the required format (i.e., with the EBL column). It has three options and can be ran from the command line as follows: 

`Rscript mappability.R -i yourFile -b blacklistFile -o outputPrefix`

where yourFile is your six-column mappability file with columns in the order of chr, start, end, F, GC and M, and blacklistFile is the blacklist region file with three columns in the order of chr, start, and end. These two options should be the full file names with extensions (e.g., .txt or .txt.gz). The option outputFile is the prefix of your output file and outputs will be written to a txt file (e.g., MboI_20Kb_el.mm10.txt). The output file is a updated mappability file with seven columns, including chr, start, end, F, GC, M and EBL, where EBL is short for "ENCODE blasklist".

If there is no blacklist region data available, you can also manually add an EBL column with all zeros so that no regions will be filtered out.

#### Define the chromosomes number and the name of chromosome size file
If gb is one of hg19, GRCh38, mm9 and mm10, chromosome size file and the number of chromosome can be omitted. Otherwise user needs to provide a chromosome size file. The file needs to be composed of two columns, the first column is the name of the chromosome, and the second is the corresponding size. The chromosomes in the file need to directly relate to the chromosomes in the process. When chromosome size file is provided, the number of chromosomes can be omitted. In this example, both parameters can be omitted.

```{r define the number of chromosomes and the name of chromosome size file}
nchrom <- 23
chroms_file <- NULL
```

#### Define the version of juicer_tools if input file is in .hic format
If the input file is in .hic format, the full name of juicer_tools needs to be provided. For example, "juicer_tools.2.20.00.ac.jar". For other format of input file, it is not required.

```{r define the juicer_tools_version,warning=FALSE}
juicer_tools_version <- NULL
```

#### Define the bin size
Default is 40000 (40 kb). Other recommended bin size are 10000 (10 kb) and 20000 (20 kb).

```{r define the bin size,warning=FALSE}
binsize <- 40000
```

#### Change the cis-interacting regions
Users have the option to change the cis-interacting region threshold. Default is 200 kb.

```{r change the cis-interacting regions,warning=FALSE}
upper_cis <- 200000
```

#### Define if the input matrices are ALREADY normalized
Default is FALSE. If true, it will skip within-normalization step.

```{r define if the input matrices are ALREADY normalized,warning=FALSE}
normalized <- FALSE
```

#### Define whether to remove MHC region

The MHC regions defined in the R package are:<br/>

| GB   | CHR  | BP RANGE  |
| ---- |:---: | ---------:|
|hg19  | chr6 |  28477797-33448354                   |
|GRCh38| chr6 | 28510120-33480577                    |
|mm9   |chr17 |33888191-35744546 & 36230820-38050373 |
|mm10  |chr17 | 33681276 & 38548659                  |

The default setting is "TRUE"", that is, to remove the MHC region if gb is in c("hg19","GRCh38","mm9","mm10"), "FALSE" if not.

```{r define whether to remove MHC region, message=FALSE}
rm_mhc <- TRUE
```

#### Define whether to remove the ENCODE black list regions

The ENCODE blacklist regions are described [here](https://www.nature.com/articles/s41598-019-45839-z). The EBL variable in the mappability file is an indicator of whether it is a black list region (1) or not (0). The ENCODE blacklist regions can be downloaded here https://sites.google.com/site/anshulkundaje/projects/blacklists. If there is no blacklist region information, please specify this option to "FALSE".

```{r define whether to remove ENCODE black list regions, message=FALSE}
rm_EBL <- TRUE
```

The default setting is "TRUE" to remove the ENCODE blacklist regions if gb is in c("hg19","GRCh38","mm9","mm10"), "FALSE" if not.

#### Change the filtering threshold 
Users are able to change the filtering threshold, where if a cis-interaction is calculated by more than 25% bins that contains a mappability of 0, a GC content of 0 , or a effective fragment length of 0,then it is also filtered.

```{r define filter, message=FALSE}
rm_perc <- 0.25
```

#### Change the regression distribution
Users have the option to change the regression distribution used. The HiCNormCis method uses Poisson ("poisson") distribution, but the user can change to negative binomial ("nb")

```{r change the regression distribution, message=FALSE}
dist <- 'poisson'
```

#### Change the alpha cutoff for the significant p-value
Default is 0.05.

```{r change the alpha cutoff for the significant p-value, message=FALSE}
alpha <- 0.05
```

#### Specify if you would like circos plots for the FIREs and super-FIREs
Default is FALSE. If users would like to incorporate other multi-omics data, see *FIRE_plot* function for details.

```{r specify if you would like circos plots for the FIREs and super-FIREs, message=FALSE}
plots <- FALSE
```

#### Specify if you would like differential FIRE analysis
Default is FALSE. If users would like to perform differential FIRE analysis, see *diff_interactions* function fro details. Differential FIRE analysis can be performed between two samples when each has at least two replicates. Naming convention between replicates needs to be similar. Sample name is ordered alphabetically.

```{r specify if you would like differential FIRE analysis, message=FALSE}
diff_fires <- FALSE
```

#### Call FIREs and super-FIREs

Using FIREcaller function, we call both FIREs and super-FIREs for 22 autosomes of Hippocampus dataset.
  
```{r call FIRE and super-FIRE for 22 autosomes of Hippocampus dataset}
FIREcaller(file.list, gb, map_file, nchrom, chroms_file=NULL, juicer_tools_version=NULL, binsize = 40000, upper_cis = 200000, normalized = FALSE, rm_mhc = TRUE, 
           rm_EBL = TRUE, rm_perc = 0.25, dist = 'poisson', alpha = 0.05, plots = FALSE, diff_fires = FALSE)
```

In this case, two sets of files will be returned: one for FIREs and the other one for super-FIREs. If multiple (*n*) prefix's are specified in the *file.list*, there will be one file for FIREs and *n* superFIRE files.

```{r An example for Fire and SuperFire outputs}
# An example for FIRE output
FIRE_output = read.table("FIRE_ANALYSIS_40000_200000_poisson.txt", header = T)
head(FIRE_output)

# An example for super-FIRE output
SuperFIRE_output = read.table("super_FIRE_call_Hippo.txt", header = T)
head(SuperFIRE_output)
```

## Citation
Crowley, C., Yang, Y.\*, Qiu, Y., Hu, B., Abnousi, A., Lipiński, J., Plewczynski, D., Wu, D., Won, H., Ren, B., Hu, M.\*, Li, Y\*. (2021) FIREcaller: Detecting Frequently Interacting Regions from Hi-C Data. *Computational and Structural Biotechnology Journal*, 19: 355–362.
