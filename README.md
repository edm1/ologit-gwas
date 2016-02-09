# ologit-gwas
Ordinal logistic regression for GWAS.

- Uses *polr* from R package *MASS* to do ordinal logistic regression
- Runs on multiple cores
- Reads genotypes in chunks to prevent memory problems

## Usage

```
Usage: Usage: Rscript ologit_gwas_v1.0.R <geno> <sample> <outfile> <chunk size> <ncores> <phenocol> <covarcols>...
```

Where
- `geno` - genotype file (format below)
- `sample` - file containing phenotype and covariate data ([format](http://www.stats.ox.ac.uk/~marchini/software/gwas/file_format.html#Sample_File_Format)). Second line variable type for phenotype column should be D
- `outfile` - output file for results
- `chunk size` - Num of genotypes to read at once
- `ncores` - Num of cores to use
- `phenocol` - Phenotype column name
- `covarcols` - Covariate column names (list separated by spaces)

### Genotype format

- Designed to work with expected minor/alternate allele counts
- Genotypes should be in same order as in the sample file (1 per individual)
- No header line
- Genotypes should start at column 7
- The first 6 columns should be:

```CHROM, SNP_id, RS_id, position, ref_allele, alt_allele```


