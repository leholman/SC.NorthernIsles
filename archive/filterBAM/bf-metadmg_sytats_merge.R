library(lvplot)
library(showtext)
library(tidyverse)
library(DescTools)
#showtext_auto()

source("~/Library/CloudStorage/GoogleDrive-mikkel.w.pedersen@gmail.com/My Drive/MediterreneanProject/sapropels/get-metadata.R")
source("~/Library/CloudStorage/GoogleDrive-mikkel.w.pedersen@gmail.com/My Drive/MediterreneanProject/sapropels/dmg.R")
source("~/Library/CloudStorage/GoogleDrive-mikkel.w.pedersen@gmail.com/My Drive/MediterreneanProject/mediterranean_R/get_calculate_plot_grid.R")
source("~/Library/CloudStorage/GoogleDrive-mikkel.w.pedersen@gmail.com/My Drive/MediterreneanProject/mediterranean_R/get_penalized_weighted_median_reads.R")
source("~/Library/CloudStorage/GoogleDrive-mikkel.w.pedersen@gmail.com/My Drive/MediterreneanProject/sapropels/perk.R")
source("~/Library/CloudStorage/GoogleDrive-mikkel.w.pedersen@gmail.com/My Drive/MediterreneanProject/sapropels/damage_est_function.R")
source("~/Library/CloudStorage/GoogleDrive-mikkel.w.pedersen@gmail.com/My Drive/MediterreneanProject/sapropels/perk_wrapper.R")
source("~/Library/CloudStorage/GoogleDrive-mikkel.w.pedersen@gmail.com/My Drive/MediterreneanProject/sapropels/perk_wrapper_function.R")
source("~/Library/CloudStorage/GoogleDrive-mikkel.w.pedersen@gmail.com/My Drive/MediterreneanProject/sapropels/get_dmg_decay_fit.R")


# Let's load the cdata
cdata <- metadata_file!!!! 

# Load metaDMG data
holi_data <- read_tsv("~/Library/CloudStorage/GoogleDrive-mikkel.w.pedersen@gmail.com/My Drive/MediterreneanProject/sapropels/tp-damage.tsv.gz") |>
  mutate(
    label = paste(str_extract(label, "(MED)-(\\d+)-(\\d+|NC\\d+|NL\\d+)"), sep = "-")
  ) |>
  inner_join(cdata |> select(label, internal_name_by_dom, specific_feature, label_fig))


# Get species and with at least 100 reads
holi_data_sp_euk_100 <- holi_data |>
  filter(rank == "species") |>
  filter(grepl("Eukaryota", taxa_path)) |>
  filter(nreads >= 100) |>
  rename(tax_name = taxid, n_reads = nreads)

#
holi_data_sp_euk_50 <- holi_data |>
  filter(rank == "species") |>
  filter(grepl("Eukaryota", taxa_path)) |>
  filter(nreads < 100 & nreads >= 50) |>
  rename(tax_name = taxid, n_reads = nreads)

# Let's get the damage fits using CCC
samples <- cdata$label |> unique()
dat100 <- dmg_fwd_CCC(holi_data_sp_euk_100, samples, ci = "asymptotic", nperm = 100, nproc = 24)
dat50 <- dmg_fwd_CCC(holi_data_sp_euk_50, samples, ci = "asymptotic", nperm = 100, nproc = 24)


# Get good fits for the 100 reads
dat100_filt <- dat100 |>
  ungroup() |>
  inner_join(holi_data_sp_euk_100) |>
  mutate(fit = ifelse(rho_c >= 0.85 & C_b > 0.9 & round(rho_c_perm_pval, 3) < 0.1 & !is.na(rho_c), "good", "bad")) |>
  mutate(fit = ifelse(q_CI_h >= 1 | c_CI_l <= 0, "bad", fit))

# Plot how many good/bad species we have per sample
pdf("proportion_good_bad_species_100readsplus.pdf", width = 13, height = 8)
dat100_filt |>
  group_by(label, fit) |>
  count() |>
  ungroup() |>
  inner_join(cdata) |>
  mutate(
    internal_name_by_dom = fct_reorder(internal_name_by_dom, -age),
    label_fig = fct_reorder(label_fig, -age)
  ) |>
  ggplot(aes(internal_name_by_dom, n, fill = fit)) +
  geom_col() +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    legend.position = "top",
    legend.title = element_blank(),
    text = element_text(size = 16),
  ) +
  ylab("Number of species") +
  xlab("Sample") +
  guides(fill = guide_legend(nrow = 1)) +
  facet_grid(~label_fig, scales = "free_x", space = "free")
dev.off()

# Get good fits for the 50 reads
dat50_filt <- dat50 |>
  ungroup() |>
  inner_join(holi_data_sp_euk_50) |>
  mutate(fit = ifelse(rho_c >= 0.85 & C_b > 0.9 & round(rho_c_perm_pval, 3) < 0.1 & !is.na(rho_c), "good", "bad")) |>
  mutate(fit = ifelse(q_CI_h >= 1 | c_CI_l <= 0, "bad", fit))

# Plot how many good/bad species we have per sample
pdf("proportion_good_bad_species_less100reads.pdf", width = 13, height = 8)
dat50_filt |>
  group_by(label, fit) |>
  count() |>
  ungroup() |>
  inner_join(cdata) |>
  mutate(
    internal_name_by_dom = fct_reorder(internal_name_by_dom, -age),
    label_fig = fct_reorder(label_fig, -age)
  ) |>
  ggplot(aes(internal_name_by_dom, n, fill = fit)) +
  geom_col() +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    legend.position = "top",
    legend.title = element_blank(),
    text = element_text(size = 16),
  ) +
  ylab("Number of species") +
  xlab("Sample") +
  guides(fill = guide_legend(nrow = 1)) +
  facet_grid(~label_fig, scales = "free_x", space = "free")
dev.off()



# Let's plots some random taxa for good and bad fits

# 100 reads, good fit
tax <- dat100_filt |>
  ungroup() |>
  filter(fit == "good") |>
  group_by(internal_name_by_dom) |>
  slice_sample(n = 100) |>
  ungroup()


samples <- tax$internal_name_by_dom |> unique()

plots100_good <- purrr::map(.x = samples, dat = tax, .f = function(x, dat, orient = "fwd", pos = 25, p_breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6)) {
  data <- dat |> filter(internal_name_by_dom == x)
  grid_size <- calculate_plot_grid(length(data$tax_name))
  l <- lapply(data$name, function(X) {
    df1 <- data |>
      filter(name == X)
    p <- get_dmg_decay_fit(df1, orient = orient, pos = pos, p_breaks = p_breaks)
    p <- p + ggtitle(X)
    return(p)
  })
  plot <- ggpubr::ggarrange(plotlist = l, ncol = grid_size$cols, nrow = grid_size$rows, align = "hv")
  tit <- paste0(x, " -- ", data$label_fig |> unique())
  ggpubr::annotate_figure(plot, top = ggpubr::text_grob(tit,
                                                        color = "black", face = "bold", size = 12
  ))
}, .progress = TRUE)

names(plots100_good) <- samples

pdf("plots100good.pdf", width = 20, height = 20)
print(plots100_good)
dev.off()

# 50 reads, good fit
tax <- dat50_filt |>
  ungroup() |>
  filter(fit == "good") |>
  group_by(internal_name_by_dom) |>
  slice_sample(n = 100) |>
  ungroup()

samples <- tax$internal_name_by_dom |> unique()

plots50_good <- purrr::map(.x = samples, dat = tax, .f = function(x, dat, orient = "fwd", pos = 25, p_breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6)) {
  data <- dat |> filter(internal_name_by_dom == x)
  grid_size <- calculate_plot_grid(length(data$tax_name))
  l <- lapply(data$name, function(X) {
    df1 <- data |>
      filter(name == X)
    p <- get_dmg_decay_fit(df1, orient = orient, pos = pos, p_breaks = p_breaks)
    p <- p + ggtitle(X)
    return(p)
  })
  plot <- ggpubr::ggarrange(plotlist = l, ncol = grid_size$cols, nrow = grid_size$rows, align = "hv")
  tit <- paste0(x, " -- ", data$label_fig |> unique())
  ggpubr::annotate_figure(plot, top = ggpubr::text_grob(tit,
                                                        color = "black", face = "bold", size = 12
  ))
}, .progress = TRUE)

names(plots50_good) <- samples

pdf("plot50good.pdf", width = 20, height = 20)
print(plots50_good)
dev.off()

# 100 reads, bad fit
tax <- dat100_filt |>
  ungroup() |>
  filter(fit == "bad") |>
  group_by(internal_name_by_dom) |>
  slice_sample(n = 10) |>
  ungroup()


samples <- tax$internal_name_by_dom |> unique()

plots100_bad <- purrr::map(.x = samples, dat = tax, .f = function(x, dat, orient = "fwd", pos = 25, p_breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6)) {
  data <- dat |> filter(internal_name_by_dom == x)
  grid_size <- calculate_plot_grid(length(data$tax_name))
  l <- lapply(data$name, function(X) {
    df1 <- data |>
      filter(name == X)
    p <- get_dmg_decay_fit(df1, orient = orient, pos = pos, p_breaks = p_breaks)
    p <- p + ggtitle(X)
    return(p)
  })
  plot <- ggpubr::ggarrange(plotlist = l, ncol = grid_size$cols, nrow = grid_size$rows, align = "hv")
  tit <- paste0(x, " -- ", data$label_fig |> unique())
  ggpubr::annotate_figure(plot, top = ggpubr::text_grob(tit,
                                                        color = "black", face = "bold", size = 12
  ))
}, .progress = TRUE)

names(plots100_bad) <- samples

pdf("plots100_bad.pdf", width = 20, height = 20)
print(plots100_bad)
dev.off()

# 50 reads, bad fit
tax <- dat50_filt |>
  ungroup() |>
  filter(fit == "bad") |>
  group_by(internal_name_by_dom) |>
  slice_sample(n = 10) |>
  ungroup() |>
  inner_join(cdata |> select(label, label_fig))

samples <- tax$internal_name_by_dom |> unique()

plots50_bad <- purrr::map(.x = samples, dat = tax, .f = function(x, dat, orient = "fwd", pos = 25, p_breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6)) {
  data <- dat |> filter(internal_name_by_dom == x)
  grid_size <- calculate_plot_grid(length(data$tax_name))
  l <- lapply(data$name, function(X) {
    df1 <- data |>
      filter(name == X)
    p <- get_dmg_decay_fit(df1, orient = orient, pos = pos, p_breaks = p_breaks)
    p <- p + ggtitle(X)
    return(p)
  })
  plot <- ggpubr::ggarrange(plotlist = l, ncol = grid_size$cols, nrow = grid_size$rows, align = "hv")
  tit <- paste0(x, " -- ", data$label_fig |> unique())
  ggpubr::annotate_figure(plot, top = ggpubr::text_grob(tit,
                                                        color = "black", face = "bold", size = 12
  ))
}, .progress = TRUE)

names(plots50_bad) <- samples

pdf("plots50_bad.pdf", width = 20, height = 20)
print(plots50_bad)
dev.off()

# Use filterBAM data to identify genera that might be spurious
library(data.table)
library(tidyverse)
library(taxonomizr)
#library(unixtools)
library(dtplyr)


# acc2taxid <- fread("~/Desktop/ncbi_taxonomy_01Oct2022/combined_accession2taxid_20221112.gz",
#                    tmpdir = "~/Downloads/",
#                    nThread = 12,
#                    showProgress = TRUE
# )
names <- "~/Desktop/ncbi_taxonomy_01Oct2022/names.dmp"
nodes <- "~/Desktop/ncbi_taxonomy_01Oct2022/nodes.dmp"
# acc2taxid_file <- "/projects/lundbeck/scratch/for_antonio/ncbi_taxonomy_01Oct2022/combined_accession2taxid_20221112.gz"


read.names.sql(names, sqlFile = "nameNode.sqlite", overwrite=TRUE)
read.nodes.sql(nodes, sqlFile = "nameNode.sqlite", overwrite=TRUE)

fb_data <- read_tsv("bam-filtered.summary.tsv.gz") 


accs <- fb_data |>
  select(accession.version = reference) |>
  distinct() |>
  pull(accession.version)

tax_ids <- acc2taxid %>%
  filter(accession.version %in% accs)

tax_ids <- read_csv("tax_ids.txt")

# write_tsv(tax_ids, "./results/taxonomy/holi-fb-taxids.tsv.gz")

tax_ids_lst <- tax_ids |>
  select(taxid) |>
  distinct() |>
  pull(taxid)

holi_data <- read_csv("~/Library/CloudStorage/GoogleDrive-mikkel.w.pedersen@gmail.com/My Drive/MediterreneanProject/mediterranean_R/final_euk_data.csv", ) |>
  select(-1)

tax_data <- getTaxonomy(tax_ids_lst, "nameNode.sqlite") |>
  as_tibble() |>
  mutate(taxid = tax_ids_lst)


# species level
agg_species_stats_less10k <- tax_ids |>
  inner_join(tax_data) |>
  inner_join(fb_data |> rename(accession.version = reference)) |>
  mutate(flt = paste(label, species, sep = "--")) |>
  #filter(flt %in% (dat100_filt |> select(label, name) |> mutate(flt = paste(label, name, sep = "--")) |> pull(flt))) |>
  inner_join(dat100_filt |> 
               select(label, name, A_b) |> 
               mutate(flt = paste(label, name, sep = "--")) |> 
               select(-label, -name)) |>
  filter(read_ani_median >= 95, n_reads >= 100, A_b > 0.05) |>
  as_tibble() |>
  select(-flt) |>
  group_by(label, species) |>
  mutate(
    num_alns = sum(n_alns),
    weight = n_alns / num_alns,
    scaled_breadth_exp_ratio = breadth_exp_ratio * weight
  ) |>
  ungroup() |>
  group_by(label, species) |>
  summarise(
    n = n(),
    median_A_b = median(A_b),
    mean_read_ani_median = mean(read_ani_median),
    mean_read_ani_std = mean(read_ani_std),
    median_reference_length = median(reference_length),
    sum_reference_length = sum(reference_length),
    mean_breadth_exp_ratio = mean(breadth_exp_ratio),
    median_breadth_exp_ratio = median(breadth_exp_ratio),
    penalized_weighted_median_entropy = penalized_weighted_median(norm_entropy, n_reads, norm_entropy),
    penalized_weighted_median_gini = penalized_weighted_median(norm_gini, n_reads, norm_gini),
    penalized_weighted_median_breadth_exp_ratio = penalized_weighted_median(breadth_exp_ratio, n_reads, breadth_exp_ratio),
    median_entropy = median(norm_entropy),
    mean_entropy = mean(norm_entropy),
    median_gini = median(norm_gini),
    mean_gini = mean(norm_gini),
    median_n_reads = median(n_reads),
    mean_n_reads = mean(n_reads),
    n_reads = sum(n_reads)
  ) |>
  ungroup() |>
  mutate(rm = ifelse(penalized_weighted_median_breadth_exp_ratio > 0.8 | (penalized_weighted_median_gini < 0.6 & penalized_weighted_median_entropy > 0.75), "keep", "remove")) |> 
  filter(sum_reference_length < 1000) |>
  inner_join(cdata) |>
  inner_join(tax_data) |>
  filter(rm == 'keep')


agg_species_stats_more10k <- tax_ids |>
  inner_join(tax_data) |>
  inner_join(fb_data |> rename(accession.version = reference)) |>
  mutate(flt = paste(label, species, sep = "--")) |>
  #filter(flt %in% (dat100_filt |> select(label, name) |> mutate(flt = paste(label, name, sep = "--")) |> pull(flt))) |>
  inner_join(dat100_filt |> 
               #select(label, name, A_b) |> 
               mutate(flt = paste(label, name, sep = "--")) |> 
               select(-label, -name)) |>
  filter(read_ani_median >= 95, n_reads >= 100, A_b > 0.05) |>
  as_tibble() |>
  select(-flt) |>
  group_by(label, species) |>
  mutate(
    num_alns = sum(n_alns),
    weight = n_alns / num_alns,
    scaled_breadth_exp_ratio = breadth_exp_ratio * weight
  ) |>
  ungroup() |>
  group_by(label, species) |>
  summarise(
    n = n(),
    median_A_b = median(A_b),
    mean_read_ani_median = mean(read_ani_median),
    mean_read_ani_std = mean(read_ani_std),
    median_reference_length = median(reference_length),
    sum_reference_length = sum(reference_length),
    mean_breadth_exp_ratio = mean(breadth_exp_ratio),
    median_breadth_exp_ratio = median(breadth_exp_ratio),
    penalized_weighted_median_entropy = penalized_weighted_median(norm_entropy, n_reads, norm_entropy),
    penalized_weighted_median_gini = penalized_weighted_median(norm_gini, n_reads, norm_gini),
    penalized_weighted_median_breadth_exp_ratio = penalized_weighted_median(breadth_exp_ratio, n_reads, breadth_exp_ratio),
    median_entropy = median(norm_entropy),
    mean_entropy = mean(norm_entropy),
    median_gini = median(norm_gini),
    mean_gini = mean(norm_gini),
    median_n_reads = median(n_reads),
    mean_n_reads = mean(n_reads),
    n_reads = sum(n_reads)
  ) |>
  ungroup() |>
  mutate(rm = ifelse(penalized_weighted_median_breadth_exp_ratio > 0.8 | (penalized_weighted_median_gini < 0.6 & penalized_weighted_median_entropy > 0.75), "keep", "remove")) |> 
  filter(sum_reference_length > 1000) |>
  inner_join(cdata) |>
  inner_join(tax_data) |>
  filter(rm == 'keep')

agg_species_stats_more10k |>
  group_by(rm) |>
  count() |> 
  knitr::kable()

agg_species_stats_less10k |>
  group_by(rm) |>
  count() |> 
  knitr::kable()

unique(agg_species_stats_less10k$species)
unique(agg_species_stats_more10k$species)

agg_species_stats_more10k %>%
  filter(rm == "keep") %>%
  summarise(unique_species_count = n_distinct(species))

agg_species_stats_less10k %>%
  filter(rm == "keep") %>%
  summarise(unique_species_count = n_distinct(species))


# Print unique genus values after filtering in the first dataset
species_filtered_more10k <- agg_species_stats_more10k %>%
  filter(rm == "keep") %>%
  distinct(species)

print(species_filtered_more10k)

# Print unique genus values after filtering in the second dataset
species_filtered_less10k <- agg_species_stats_less10k %>%
  filter(rm == "keep") %>%
  distinct(species)

print(species_filtered_less10k)

# Find unique values present in the first print but not in the second
unique_in_first <- setdiff(species_filtered_more10k$species, species_filtered_less10k$species)
print(unique_in_first)

# Find unique values present in the second print but not in the first
unique_in_second <- setdiff(species_filtered_less10k$species, species_filtered_more10k$species)
print(unique_in_second)

# genus level
agg_genus_stats_less10k <- tax_ids |>
  inner_join(tax_data) |>
  inner_join(fb_data |> rename(accession.version = reference)) |>
  mutate(flt = paste(label, species, sep = "--")) |>
  #filter(flt %in% (dat100_filt |> select(label, name) |> mutate(flt = paste(label, name, sep = "--")) |> pull(flt))) |>
  inner_join(dat100_filt |> 
               select(label, name, A_b, taxa_path) |> 
               mutate(flt = paste(label, name, sep = "--")) |> 
               select(-label, -name)) |>
  filter(read_ani_median >= 95, n_reads >= 100, A_b > 0.05) |>
  as_tibble() |>
  select(-flt) |>
  group_by(label, genus) |>
  mutate(
    num_alns = sum(n_alns),
    weight = n_alns / num_alns,
    scaled_breadth_exp_ratio = breadth_exp_ratio * weight
  ) |>
  ungroup() |>
  group_by(label, genus) |>
  summarise(
    n = n(),
    median_A_b = median(A_b),
    mean_read_ani_median = mean(read_ani_median),
    mean_read_ani_std = mean(read_ani_std),
    median_reference_length = median(reference_length),
    sum_reference_length = sum(reference_length),
    mean_breadth_exp_ratio = mean(breadth_exp_ratio),
    median_breadth_exp_ratio = median(breadth_exp_ratio),
    penalized_weighted_median_entropy = penalized_weighted_median(norm_entropy, n_reads, norm_entropy),
    penalized_weighted_median_gini = penalized_weighted_median(norm_gini, n_reads, norm_gini),
    penalized_weighted_median_breadth_exp_ratio = penalized_weighted_median(breadth_exp_ratio, n_reads, breadth_exp_ratio),
    median_entropy = median(norm_entropy),
    mean_entropy = mean(norm_entropy),
    median_gini = median(norm_gini),
    mean_gini = mean(norm_gini),
    median_n_reads = median(n_reads),
    mean_n_reads = mean(n_reads),
    n_reads = sum(n_reads)
  ) |>
  ungroup() |>
  mutate(rm = ifelse(penalized_weighted_median_breadth_exp_ratio > 0.8 | (penalized_weighted_median_gini < 0.6 & penalized_weighted_median_entropy > 0.75), "keep", "remove")) |> 
  filter(sum_reference_length < 1000) |>
  inner_join(cdata) |>
  inner_join(tax_data) |>
  filter(rm == 'keep')

agg_genus_stats_more10k <- tax_ids |>
  inner_join(tax_data) |>
  inner_join(fb_data |> rename(accession.version = reference)) |>
  mutate(flt = paste(label, species, sep = "--")) |>
  #filter(flt %in% (dat100_filt |> select(label, name) |> mutate(flt = paste(label, name, sep = "--")) |> pull(flt))) |>
  inner_join(dat100_filt |> 
               select(label, name, A_b) |> 
               mutate(flt = paste(label, name, sep = "--")) |> 
               select(-label, -name)) |>
  filter(read_ani_median >= 95, n_reads >= 100, A_b > 0.05) |>
  as_tibble() |>
  select(-flt) |>
  group_by(label, genus) |>
  mutate(
    num_alns = sum(n_alns),
    weight = n_alns / num_alns,
    #median_A_b = median(A_b),
    scaled_breadth_exp_ratio = breadth_exp_ratio * weight
  ) |>
  ungroup() |>
  group_by(label, genus) |>
  summarise(
    n = n(),
    median_A_b = median(A_b),
    mean_read_ani_median = mean(read_ani_median),
    mean_read_ani_std = mean(read_ani_std),
    median_reference_length = median(reference_length),
    sum_reference_length = sum(reference_length),
    mean_breadth_exp_ratio = mean(breadth_exp_ratio),
    median_breadth_exp_ratio = median(breadth_exp_ratio),
    penalized_weighted_median_entropy = penalized_weighted_median(norm_entropy, n_reads, norm_entropy),
    penalized_weighted_median_gini = penalized_weighted_median(norm_gini, n_reads, norm_gini),
    penalized_weighted_median_breadth_exp_ratio = penalized_weighted_median(breadth_exp_ratio, n_reads, breadth_exp_ratio),
    median_entropy = median(norm_entropy),
    mean_entropy = mean(norm_entropy),
    median_gini = median(norm_gini),
    mean_gini = mean(norm_gini),
    median_n_reads = median(n_reads),
    mean_n_reads = mean(n_reads),
    n_reads = sum(n_reads)
  ) |>
  ungroup() |>
  mutate(rm = ifelse(penalized_weighted_median_breadth_exp_ratio > 0.8 | (penalized_weighted_median_gini < 0.6 & penalized_weighted_median_entropy > 0.75), "keep", "remove")) |> 
  filter(sum_reference_length > 1000) |>
  inner_join(cdata) |>
  inner_join(tax_data) |>
  filter(rm == 'keep')


agg_genus_stats_more10k |>
  group_by(rm) |>
  count() |> 
  knitr::kable()

agg_genus_stats_less10k |>
  group_by(rm) |>
  count() |> 
  knitr::kable()


unique(agg_genus_stats_less10k$genus)
unique(agg_genus_stats_more10k$genus)


agg_genus_stats_more10k %>%
  filter(rm == "keep") %>%
  summarise(unique_genus_count = n_distinct(genus))

agg_genus_stats_less10k %>%
  filter(rm == "keep") %>%
  summarise(unique_genus_count = n_distinct(genus))


# Print unique genus values after filtering in the first dataset
genus_filtered_more10k <- agg_genus_stats_more10k %>%
  filter(rm == "keep") %>%
  distinct(genus)

print(genus_filtered_more10k)

# Print unique genus values after filtering in the second dataset
genus_filtered_less10k <- agg_genus_stats_less10k %>%
  filter(rm == "keep") %>%
  distinct(genus)

print(genus_filtered_less10k)

# Find unique values present in the first print but not in the second
unique_in_first <- setdiff(genus_filtered_more10k$genus, genus_filtered_less10k$genus)
print(unique_in_first)

# Find unique values present in the second print but not in the first
unique_in_second <- setdiff(genus_filtered_less10k$genus, genus_filtered_more10k$genus)
print(unique_in_second)

pdf(file = "reference_lengths_sum_median.pdf")
agg_genus_stats_more10k %>%
  ggplot(aes(x = penalized_weighted_median_breadth_exp_ratio, y = penalized_weighted_median_gini, size = median_reference_length, color = rm)) +
  geom_point()

agg_genus_stats_more10k %>%
  ggplot(aes(x = penalized_weighted_median_breadth_exp_ratio, y = penalized_weighted_median_gini, size = sum_reference_length, color = rm)) +
  geom_point()
dev.off()







