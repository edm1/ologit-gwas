#!/usr/bin/env Rscript
#
# Usage: Rscript ologit_gwas_v1.0.R <geno> <sample> <outfile> <chunk size> <ncores> <phenocol> <covarcols>...
#

# Load libraries
library("MASS")
library("foreach")
library("doMC")

# Get command line args
args = commandArgs(trailingOnly=TRUE)
ingen = args[1]
insample = args[2]
outres = args[3]
chunksize = as.integer(args[4])
ncores = as.integer(args[5])
phenocol = args[6]
covarcols = args[7:length(args)]
# print(insample)

# Set arguments
# ingen = "example_data/data_maf0.01_info0.8_chr21_head1000.add.gz"
# insample = "example_data/all_phenos.sample"
# outres = "out_temp.tsv"
# phenocol = "motorABC_f7_balanceScore"
# covarcols = c("covars_f7_age", "covars_all_sex", "covars_all_pc1", "covars_all_pc2")
# chunksize = 10000
# ncores = 16

# Set number of cores
registerDoMC(ncores)

#
# Load sample file
#

# Load coltypes from top of file
coltypes = read.table(insample, sep="\t", header=T, stringsAsFactors=F, nrows=1)
# Load the rest of the sample
sample = read.table(insample, sep="\t", header=F, stringsAsFactors=F, skip=2)
colnames(sample) = colnames(coltypes)
# Keep only required cols
sample = sample[, c(phenocol, covarcols)]
coltypes = coltypes[, c(phenocol, covarcols)]
# Set discrete vars as factors
sample[, coltypes == "D"] = lapply(sample[, coltypes == "D", drop=F], factor)
# Set continuous vars as numerics
sample[, coltypes == "C"] = lapply(sample[, coltypes == "C", drop=F], as.numeric)

# Get list of complete.cases
to_keep = complete.cases(sample)

# Split into pheno and covars
pheno = sample[to_keep, phenocol]
covars = sample[to_keep, covarcols, drop=F]

# Open outfile and write headers
out_h = file(outres, open="w")
write.table(t(c("SNP_id", "pos", "alleleA", "alleleB", "maf", "pvalue")), file=out_h, sep="\t", row.names=F, quote=F, col.names=F)

#
# Load and process genetic data in chunks
#

# Open geno file
ingeno_h = gzfile(ingen, open="r")

# Read in first chunk for geno
geno = read.table(ingeno_h, nrows=chunksize, sep=" ", header=F, stringsAsFactors=F)

# Repeat processing until there are no chunks left
repeat {
  
  # Apply to each snp
  result = foreach(i=1:nrow(geno), .combine=rbind) %dopar% {
    
    #i=1
    
    # Get snp additive data
    snp_add = as.numeric(geno[i, c(7:ncol(geno))])
    snp_add = snp_add[to_keep]
    
    # Calc maf
    maf = sum(snp_add) / (2 * length(snp_add))
    if (maf > 0.5) { maf = 1 - maf }
    
    # Caluclate assoication pval
    snp_add_p = tryCatch(
    {
      
      # Do regression
      m = polr(pheno ~ . + snp_add, data=covars, Hess=T)
      
      # Calculate P-values
      ctable = coef(summary(m))
      p = pnorm(abs(ctable[, "t value"]), lower.tail = FALSE) * 2
      ctable = cbind(ctable, "p value" = p)
      snp_add_p = ctable["snp_add", "p value"]
      snp_add_p
    },
    
    error = function(e) {
      return(NA)
    }
    )
    
    
    # Out line
    list(geno[i, 3], geno[i, 4], geno[i, 5], geno[i, 6], maf, snp_add_p)
  }
  
  # Write results to outfile
  write.table(result, file=out_h, sep="\t", row.names=F, quote=F, col.names=F)
  
  # Load next chunk of data
  if (nrow(geno) == 0)
    break
  ## process chunk 'data' here, then...
  ## ...read next chunk
  if (nrow(geno) != chunksize)   # last chunk was final chunk
    break
  geno <- tryCatch({
    read.table(ingeno_h, nrows=chunksize, sep=" ", header=F, stringsAsFactors=F)
  }, error=function(err) {
    ## matching condition message only works when message is not translated
    if (identical(conditionMessage(err), "no lines available in input"))
      data.frame()
    else stop(err)
  })
}

# Close input file
close(ingeno_h)  

# Close output file
close(out_h)
