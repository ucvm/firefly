# Firefly

Dependable pipeline for processing 16S amplicon data


## Pipeline

Custom parameters stored in `config.yaml`

### Preprocessing

* Clip off primers with cutadapt
* Check quality with fastqc, using multiqc to compile report
* Quality filter with dada2::filterAndTrim using max expected error

### Dada2 pipeline

* Learn error with whole dataset: `learn_errors.R`
* De-replicate and infer sequences: `infer_seqs.R`
* Remove bimeras, assign taxonomy (including species): `assign_taxonomy.R`

### Make tree
* Align seqs with `ssu-align`
* Make tree with `FastTree`




