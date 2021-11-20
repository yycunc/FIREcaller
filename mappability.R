suppressMessages(library(data.table),library(argparse))
library(data.table)
library(argparse)

setup_parser <- function() {
  parser <- argparse::ArgumentParser()
  parser$add_argument('-i', '--input', help = "input mappability file with 6 columns chr, start, end, F, GC, M", type="character", required = TRUE)
  parser$add_argument('-o', '--output', help = "output mappability file with 7 columns chr, start, end, F, GC, M, EBL", type="character", required = TRUE)
  parser$add_argument('-b', '--blacklist', help = "Encode blacklist regions file", type="character", required = TRUE)
  return (parser)
}

parse_arguments <- function(parser) {
  args <- parser$parse_args()
  message(paste('reading input from directory:', args$input))
  message(paste('output will be written to:', args$output))
  message(paste('blacklist region file specified:', args$blacklist))
  return (args)
}


parser <- setup_parser()
args <- parse_arguments(parser)

map = fread(args$input)
bl = fread(args$blacklist)

if ( ncol(map) != 6 ) {
    stop('The nput file must have 6 columns')
}

if ( !identical(colnames(map), c("chr","start","end","F","GC","M")) ) {
    warning('Please double check if THE input file has columns in this order: chr, start, end, F, GC, M')
}

colnames(map)=c("chr","start","end","F","GC","M")
colnames(bl)=c("chr","start","end")

dat=data.frame()

for (chr_num in unique(map$chr)) { 

  bl_chr=subset(bl, chr==chr_num)
  map_chr=subset(map, chr==chr_num)
  
  x = data.table(bl_chr)
  y = data.table(map_chr)
  
  y$EBL = 0
  setkey(x, start, end)
  setkey(y, start, end)
  index=foverlaps(x, y, type="any", which=TRUE, nomatch=NULL)
  y$EBL[index$yid]=1
  
  dat=rbind(dat, y)
}

write.table(dat, args$output, row.names=F, quote=F, sep="\t")


