# Making the .sites files (one per chromosome) for analysis using ARGweaver software

library(tidyverse)
library(bigsnpr)
library(future)
library(future.apply)
library(glue)
db <- DBI::dbConnect(RSQLite::SQLite(), "data/derived/annotations.sqlite3")
snp_db <- tbl(db, "variants")


# Function to make one .sites file, for a given vector of SNPs
make_sites_file <- function(chromosome, snps_to_include){

  marker_IDs <- as.character(bigsnp_attached$map$marker.ID)

  # Limit to chromosome and set of SNPs in the GWAS
  focal_snps <- snp_db %>%
    filter(chr == chromosome & SNP %in% snps_to_include) %>%
    filter(SNP %in% marker_IDs) %>% # Double check everything matches as it should
    select(SNP, position, minor_allele, major_allele) %>% # don't change this line, indexed below as 3:4
    collect(n = Inf) %>% distinct() %>%
    mutate(position = as.numeric(position))

  # Make sure the SNPs in focal_snps are in the same order as in the DGRP genotypes file
  focal_snps <- tibble(SNP = marker_IDs) %>%
    left_join(focal_snps, by = "SNP") %>%
    filter(!is.na(minor_allele))

  column_indexes_focal_SNPs <- which(marker_IDs %in% (focal_snps %>% pull(SNP)))

  genotypes <- (bigsnp_attached$genotypes[, column_indexes_focal_SNPs])
  genotypes[is.na(genotypes)] <- "N"

  plan("multiprocess")
  position_geno <-
    future_lapply(1:ncol(genotypes),
                  function(i){
                    alleles <- as.character(focal_snps[i, 3:4]) # minor then major allele
                    position <- format(focal_snps$position[i], scientific = FALSE)
                    genos_at_locus_i <- genotypes[,i]
                    genos_at_locus_i[genos_at_locus_i == "0"] <- alleles[1]
                    genos_at_locus_i[genos_at_locus_i == "2"] <- alleles[2]
                    c(position, paste0(genos_at_locus_i, collapse = ""))
                  }) %>%
    do.call("rbind", .)

  # Create the .sites file:
  NAMES <- c("NAMES", gsub("line_", "L", bigsnp_attached$fam$family.ID))
  REGION <- c("REGION", "chr", 2, max(focal_snps$position)) # just write 2 for all chromosomes: Argweaver doesn't accept letters as chr names!
  writeLines(
    paste(paste0(NAMES, collapse = "\t"),
          paste0(REGION, collapse = "\t"),
          paste0(apply(position_geno, 1, function(x) paste0(x, collapse = "\t")),
                 collapse = "\n"),
          sep = "\n"),
    con = paste("data/argweaver/sites_files/DGRP_", chromosome, ".sites", sep = ""))
}



# Get list of all SNPs (not indels or MNPs) included in our GWAS
snps_to_include <- tbl(db, "univariate_lmm_results") %>% pull(SNP)
snps_to_include <- snps_to_include[!str_detect(snps_to_include, "DEL")]
snps_to_include <- snps_to_include[!str_detect(snps_to_include, "INS")]
snps_to_include <- snps_to_include[!str_detect(snps_to_include, "MNP")]

# Access the DGRP's bed/bim/fam files, for all the >200 lines
bigsnp <- snp_readBed("data/input/dgrp2.bed", backingfile = tempfile())
bigsnp_attached <- snp_attach(bigsnp)

# Set up a ".sites" file, needed by Argweaver. It's basically an alignment.
lapply(c("2L", "2R", "3L", "3R", "X"), make_sites_file, snps_to_include = snps_to_include)

#########################################################################################################
#########################################################################################################
# Part 2: Write some scripts to run Argweaver on the supercomputer (it uses a lot of memory, but doesn't take too long)
library(tidyverse)
library(glue)

make_slurm_scripts <- function(chr){

  chromosome_recomb_rates <-
    data.frame(
      chr = c("2L", "2R", "3L", "3R", "X"),
      r   = c("2.14e-8", "2.77e-8", "2.2e-8", "1.97e-8", "3.09e-8"),
      stringsAsFactors = FALSE)


  recombination_rate <- chromosome_recomb_rates$r[chromosome_recomb_rates$chr == chr]

  glue(
    "
    #!/bin/bash
    #SBATCH -N 1  # number of nodes
    #SBATCH -n 1  # number of cores
    #SBATCH --mem 50000
    #SBATCH --time 300:00:00
    #SBATCH -o slurm.%N.%j.out # STDOUT
    #SBATCH -e slurm.%N.%j.err # STDERR

    module load Python
    cd /data/projects/punim0243/argweaver
    arg-sample -s sites_files/DGRP_{chr}.sites --resume -N 1e6 -r {recombination_rate} -m 1.8e-8 --ntimes 20 --maxtime 200e3 -c 10 -n 1000 -o output/outfile_{chr}
    arg-extract-tmrca output/outfile_{chr}.%d.smc.gz > output/TMRCA_{chr}.tmrca.txt") %>%
    writeLines(con = glue("data/argweaver/scripts/{chr}.slurm"))
}


lapply(c("2L", "2R", "3L", "3R", "X"), make_slurm_scripts)
