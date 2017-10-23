# Filter and Trim one sample

filt = dada2::filterAndTrim(
    fwd = snakemake@input$r1, filt = snakemake@output$r1,
    rev = snakemake@input$r2, filt.rev = snakemake@output$r2,
    truncLen = snakemake@config$truncLen,
    trimLeft = snakemake@config$trimLeft,
    maxEE = snakemake@config$max_expected_error,
    verbose = TRUE, 
    multithread = snakemake@threads[[1]]
)

saveRDS(filt, snakemake@output$filt_out)

