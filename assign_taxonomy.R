############################################################
#                                                          #
#          Sequence table, chimeras and taxonomy           #
#                                                          #
############################################################


library(dada2)

n_threads = snakemake@threads[[1]]
message("Running with ", n_threads, " provided cores")

mergers = readRDS(snakemake@input[[1]])

message("Making sequence table")
seqtab = makeSequenceTable(mergers)

message("Removing chimeras")
seqtab = removeBimeraDenovo(seqtab, method = "consensus", multithread = n_threads, verbose = TRUE)

message("Assigning taxonomy")
rdp_tax = assignTaxonomy(seqtab, "training/rdp_train_set_16.fa.gz", multithread = n_threads, 
						 verbose = TRUE, tryRC = TRUE)
silva_tax = assignTaxonomy(seqtab, "training/silva_nr_v128_train_set.fa.gz", multithread = n_threads, 
						   verbose = TRUE, tryRC = TRUE)

message("Assigning species")
rdp_species = addSpecies(rdp_tax, "training/rdp_species_assignment_16.fa.gz", 
						 verbose = TRUE, allowMultiple = TRUE)
silva_species = addSpecies(silva_tax, "training/silva_species_assignment_v128.fa.gz", 
						   verbose = TRUE, allowMultiple = TRUE)

# --- Sequences and OTU ids
seqs = Biostrings::DNAStringSet(colnames(seqtab))
otus = paste0("OTU_", 1:ncol(seqtab))
names(seqs) = otus
colnames(seqtab) = otus

rownames(rdp_species) = otus
rownames(rdp_tax) = otus
rownames(silva_species) = otus
rownames(silva_tax) = otus

# --- Write it out
Biostrings::writeXStringSet(seqs, snakemake@output$otus)
save(seqtab, rdp_tax, rdp_species, silva_species, silva_tax, file = snakemake@output$taxonomy)
