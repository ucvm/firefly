# Infer sequence variants with DADA2

library(dada2)
library(purrr)

n_threads = snakemake@threads[[1]]
message("Running with ", n_threads, " provided cores")

# --- load error rate results

load(snakemake@input[[1]])

# --- Dereplicate, infer sequences, and merge pairs

mergers = map(r1_names, function(x) {
    message("Processing ", x)
    
    # Process forward read
    derep_r1 = derepFastq(r1_files[[x]])
    dd_r1 = dada(derep_r1, err = err_r1, multithread = n_threads)

    # Process reverse read
    derep_r2 = derepFastq(r2_files[[x]])
    dd_r2 = dada(derep_r2, err = err_r2, multithread = n_threads)

    # Merge Pairs
    merger = mergePairs(dd_r1, derep_r1, dd_r2, derep_r2)

    return(merger)
})

names(mergers) = r1_names

saveRDS(mergers, file = snakemake@output[[1]])

