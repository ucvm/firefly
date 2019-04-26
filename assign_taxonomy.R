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

train_dir = snakemake@params$train_dir

train = list(
 rdp = file.path(train_dir, "rdp_train_set_16.fa.gz"),
 silva = file.path(train_dir, "silva_nr_v132_train_set.fa.gz")
 )

tax = map(train, ~assignTaxonomy(seqtab, .x, multithread = n_threads, verbose = TRUE, tryRC = TRUE)

message("Assigning species")

spec_train = list(
    rdp = file.path(train_dir, "rdp_species_assignment_16.fa.gz"),
    silva =file.path(train_dir, "silva_species_assignment_v132.fa.gz")
) 

species = map2(tax, spec_train, ~addSpecies(.x, .y, verbose = TRUE, allowMultiple = TRUE)

               
# --- Sequences and OTU ids
            
seqs = Biostrings::DNAStringSet(colnames(seqtab))
otus = paste0("OTU_", 1:ncol(seqtab))
names(seqs) = otus
colnames(seqtab) = otus

rownames(tax$rdp) = otus
rownames(species$rdp) = otus
rownames(tax$silva) = otus
rownames(species$silva) = otus

# --- Write it out
Biostrings::writeXStringSet(seqs, snakemake@output$otus)
save(seqtab, tax, species, file = snakemake@output$taxonomy)
