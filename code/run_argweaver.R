library(glue)
library(rslurm)



run_argweaver <- function(chr){

  chromosome_recomb_rates <-
    data.frame(
      chr = c("2L", "2R", "3L", "3R", "X"),
      r   = c("2.14e-8", "2.77e-8", "2.2e-8", "1.97e-8", "3.09e-8"),
      stringsAsFactors = FALSE)

  wd <- "/data/projects/punim0243/argweaver"

  run_command <- function(shell_command){
    cat(system(glue("cd ", wd, "\n", shell_command), intern = TRUE), sep = '\n')
  }

  recombination_rate <- chromosome_recomb_rates$r[chromosome_recomb_rates$chr == chr]

  run_command(glue(
    "arg-sample -s sites_files/DGRP_{chr}.sites --overwrite -N 1e6 -r {recombination_rate} -m 1.8e-8 --ntimes 20 --maxtime 200e3 -c 10 -n 10000 -o output/outfile_{chr}
    arg-extract-tmrca output/outfile_{chr}.%d.smc.gz > output/TMRCA_{chr}.tmrca.txt"))
}

########## TESTING THE RESULTS....
parameters <- data.frame(chr = c("2L", "2R", "3L", "3R", "X"))

xx <- read.delim("data/argweaver/output/TMRCA_2R.tmrca.txt", sep = "\t", header = FALSE)[,-1] %>% as_tibble()
names(xx) <- c("start", "end", "age", "lwr", "upr")

snps <- left_join(tbl(db, "variants") %>% filter(chr == "2R") %>% select(SNP, MAF, position, site.class),
          all_snps, by = "SNP") %>% arrange(position) %>% collect(n=Inf) %>% mutate(age = NA, position = as.numeric(position))

for(i in 1:nrow(xx)){
  poses <- (xx$start[i]):(xx$end[i])
  snps$age[snps$position %in% poses] <- xx$age[i]
  if(i %% 1000 == 0) print(i)
}
##############3


objects_to_export <- ls()

sopt <- list(time = '96:00:00',   # max run time per node in hours
             mem  = '50000')      # 50MB memory per node

sjob <- slurm_apply(
  f = run_argweaver,
  params = parameters,
  add_objects = objects_to_export,
  jobname = "data_61",
  nodes = nrow(parameters),
  cpus_per_node = 1,
  slurm_options = sopt)
