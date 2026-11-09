
dat <- read_tsv("cleandata/MTG_euk_master_obs.tsv.gz") |>
  left_join(read_tsv("cleandata/MTG_euk_master_samples.tsv.gz"), by = "label") |>
  left_join(read_tsv("cleandata/MTG_euk_master_taxa.tsv.gz"), by = "taxid")
