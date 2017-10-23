# Learn errors 

library(dada2)
library(purrr)
library(stringr)

n_threads = snakemake@threads[[1]]
message("Running with ", n_threads, " provided cores")

# --- Setup file names

r1_files = snakemake@input$r1
r2_files = snakemake@input$r2

r1_names = str_split_fixed(basename(r1_files), "_", 2)[,1]
r2_names = str_split_fixed(basename(r2_files), "_", 2)[,1]

names(r1_files) = r1_names
names(r2_files) = r2_names

if (!identical(r1_names, r2_names)) stop("R1 and R2 don't match")

# --- Learn errors

err_r1 = learnErrors(r1_files, multithread = n_threads, nreads = snakemake@config$nreads)
err_r2 = learnErrors(r2_files, multithread = n_threads, nreads = snakemake@config$nreads) 

# --- Save output
save(r1_files, r1_names, err_r1, 
     r2_files, r2_names, err_r2,
     file = snakemake@output[[1]])



