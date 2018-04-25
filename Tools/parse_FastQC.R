# -----------------------------------------
# ParseQC script, by R.Gacesa (UMCG, 2018)
# -> parses through FastQC files before and after 
# cleaning by kneaddata, scythe or the like
# and returns QC stats
# -> uses fastqcr library
# -----------------------------------------
library(fastqcr)
library(tidyr)

# debug (remove later)
#setwd('C:/Users/ranko/Documents/UMCG/DAG3_codes/')

# parse CL
args = commandArgs(trailingOnly=TRUE)

# for coding only
#args[1] = './testdata/qc_preclean'
#args[2] = './testdata/qc_postclean'
#args[3] = './testdata/qc_stats.csv'

if (length(args) !=3 ) {
  stop("Call the script with <preQC path> <postQC path> <output file> \n 1 & 2 should be folders, not files, 3rd should be .csv output
       \n example: ./qc_preclean ./qc_postclean ./qc_stats.csv", call.=FALSE)
}
print ('Processing started')
# load fastqc files:
#  -> pre QC
print (' --> aggregating qc data (pre-qc)')
preQC = qc_aggregate(qc.dir = args[1], progressbar = F)
preQcWide <- spread(preQC,module,status)
preQcWide$sample <- gsub(".fq","",preQcWide$sample)
preQcWide$sample <- gsub(".gz","",preQcWide$sample)
colnames(preQcWide) <- c("sample","tot.seq","seq.lt","pct.gc","pct.dup","adpt.cont","basic.stat","overrep.seq",
                      "N.cont","seq.cont","seq.qual","seq.GC","seq.qual2","tile.qual","seq.dupllvl","seq.lt","seqs.dedup")
preQcWide$tot.seq <- as.numeric(preQcWide$tot.seq )
preQcWide$pct.gc <- as.numeric(preQcWide$pct.gc)
preQcWide$pct.dup <- as.numeric(preQcWide$pct.dup)
preQcWide$seq.dedup <- preQcWide$tot.seq * (1.0-preQcWide$pct.dup/100.0)
preQcWide$seq.lt <- as.numeric(preQcWide$seq.lt)

preQcWideA <- preQcWide
preQcWideA$sample <- gsub('_1','',preQcWideA$sample)
preQcWideA$sample <- gsub('_2','',preQcWideA$sample)
preQcWideA <- preQcWideA[,c(1:5)]
preQcWideA2 <- aggregate(data=preQcWideA,FUN = mean,.~sample)
preQcWideA3 <- aggregate(data=preQcWideA,FUN = sum,.~sample)
preQcWideA <- preQcWideA2
preQcWideA$tot.seq <- preQcWideA3$tot.seq
colnames(preQcWideA) <- c("sample","QC.raw.reads","QC.raw.read.lt","QC.raw.gc","QC.raw.duplication")
#  -> post QC
print (' --> aggregating qc data (post-qc)')
postQC = qc_aggregate(qc.dir = args[2], progressbar = F)
postQcWide <- spread(postQC,module,status)
postQcWide$sample <- gsub(".fq","",postQcWide$sample)
postQcWide <- postQcWide[,2:5]
colnames(postQcWide) <- c("QC.cleaned.reads","QC.cleaned.read.lt","QC.cleaned.gc","QC.cleaned.duplication")
# -> merge
print (' --> merging and saving output')
mergedQC <- cbind(preQcWideA,postQcWide)
# save output
write.table(x=mergedQC,file=args[3],row.names = F,col.names = T,sep=',')
print ('DONE!')
