#### Metagenome data prep ####
#### Oct 2026 ####
# 0. Setup ---------------------------------------------------------------------

## 0.1 Packages ----
suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
  library(vegan)
  library(RColorBrewer)
})

## 0.2 Paths ----
PATH_RAW  <- "rawdata/metagenomes_postpipe"                 # folder of *.agg.stat.gz
PATH_DMGR <- "other/dmg.R"                   # your dmg_fwd_CCC()
PATH_META <- "metadata.xlsx"
PATH_SITE <- "metadata_site.csv"
## come from metadata.xlsx columns named here (NULL = not available).
META_COL_RAW      <- "raw reads"
META_COL_FILTERED <- "fastP filtered reads"
META_COL_MAPPED   <- NULL
PATH_CTRL         <- "cleandata/controls"
REPORT_CTRL       <- c("extraction_blank", "library_blank")   # positives are ignored

## 0.3 Damage params ----
MIN_READS_FIT <- 100

# settings for dmg.R
CCC_RHO <- 0.75
CCC_CB  <- 0.90
CCC_P   <- 0.10
NPERM   <- 1000      # permutations inside dmg_fwd_CCC. 100 -> p in steps of 0.01
NPROC   <- 9       # workers for dmg_fwd_CCC
DMG_SEED <- 1234    

## "Damaged" / authenticated-ancient threshold on A_b: excess C->T at position 0.
DAMAGE_THRESHOLD <- 0.10
REQUIRE_GOOD_FIT <- TRUE

## 0.4 taxonomy 
RANKS <- c("domain", "kingdom", "phylum", "class", "order", "family",
           "genus", "species")
SUBSETS <- c(metazoa = 33208L, eukaryota = 2759L)
MATRIX_RANK <- "genus"
MODERN_AGE_YRBP <- -73

## 0.5 -- decontamination-
DECONTAM       <- TRUE
NEG_CONTROLS   <- NULL       # NULL = auto-detect library_blank + extraction_blank
DECONTAM_COL   <- "n_reads"  # read column compared against the blanks
BLANK_STAT     <- "max"      # "max" = beat the worst blank; "mean" = the average
BLANK_MULTIPLIER <- 1        # 1 implements exactly "fewer reads than the control"

## 0.6 core data 
CORES_BOTTOM_TO_TOP <- c("GC017", "PC012", "PC014")
CORE_ORDER_PLOT     <- rev(CORES_BOTTOM_TO_TOP)
CORE_COLS <- c(GC017 = "#D55E00", PC012 = "#0072B2", PC014 = "#009E73")

## 0.7 helping functions

say <- function(...) cat("\n>>> ", sprintf(...), "\n", sep = "")

add.alpha <- function(col, alpha = 1){
  apply(sapply(col, col2rgb) / 255, 2,
        function(x) rgb(x[1], x[2], x[3], alpha = alpha))
}

## point size on a log scale, so 100 and 100,000 reads are visibly different
scale_cex <- function(x, to = c(0.35, 3.0), floor_at = 1){
  lx  <- log10(pmax(x, floor_at))
  rng <- range(lx, finite = TRUE)
  if (!is.finite(diff(rng)) || diff(rng) == 0) return(rep(mean(to), length(lx)))
  to[1] + (lx - rng[1]) / diff(rng) * diff(to)
}

## TRUE only where the value is really TRUE; NA counts as FALSE. Used wherever a
## logical column is filtered on, because `filter(x)` silently drops NA rows and
## `sum(x)` returns NA.
is_true <- function(x) !is.na(x) & x

stopifnot(dir.exists(PATH_RAW), file.exists(PATH_META), file.exists(PATH_SITE))

# 1. Data in ---------------------------------------------------------------------

## 1.1 Read file names

files <- list.files(PATH_RAW, "\\.agg\\.stat(\\.gz)?$", full.names = TRUE)
stopifnot(length(files) > 0)
labels_from_file <- sub("\\.agg\\.stat(\\.gz)?$", "", basename(files))

say("Found %d files in %s", length(files), PATH_RAW)
print(tibble(file = basename(files), label = labels_from_file), n = Inf)

## 1.2 Check headers

headers <- vapply(files, function(f) readLines(f, n = 1), character(1))
if (n_distinct(headers) != 1)
  stop("Files do not share a header -- do not concatenate them blind.")
cn   <- strsplit(headers[[1]], "\t", fixed = TRUE)[[1]]
NPOS <- sum(grepl("^fwK[0-9]+$", cn))  # damage positions: 15
n_id <- match("A", cn) - 1L                 # width of the identity block
say("%d columns, %d damage positions (fwK0..fwK%d)", length(cn), NPOS, NPOS - 1)

## 1.3 read files in 

read_one <- function(f, lab){
  x  <- readLines(f)[-1]
  nf <- str_count(x, "\t") + 1L
  x[nf == n_id] <- str_c(x[nf == n_id], strrep("\t", length(cn) - n_id))
  bad <- which(str_count(x, "\t") + 1L != length(cn))
  if (length(bad))
    stop(sprintf("%s: line %d has %d fields, expected %d or %d",
                 basename(f), bad[1] + 1L, nf[bad[1]], n_id, length(cn)))
  read_tsv(I(x), col_names = cn, quote = "", na = c("", "NA"),
           show_col_types = FALSE, progress = FALSE) |>
    mutate(label = lab, .before = 1)
}
raw <- map2(files, labels_from_file, read_one) |> list_rbind()

## metaDMG wraps its text fields in literal quote characters: the rank column
## arrives as "species", WITH the quotes, not as species. Strip them.
## Note taxa_path keeps its inner quotes whatever readr's `quote` is set to,
## because the field starts with a digit -- so this cannot be left to the reader.
raw <- raw |> mutate(across(c(name, rank, lca, taxa_path), ~ str_remove_all(.x, '"')))

### Some checks

## Verified, not assumed. If quotes survive, `rank == "species"` matches nothing
## and every rank column downstream comes out NA -- 600 lines later, silently.
if (any(str_detect(c(raw$rank, raw$name, raw$taxa_path), '"'), na.rm = TRUE))
  stop("Quote characters survived the strip -- rank values look like '\"species\"'. ",
       "Nothing downstream will match. Check the mutate above.")
say("Quote strip verified. rank values look like: %s",
    str_c(head(unique(raw$rank), 4), collapse = ", "))

say("Read %s rows from %d samples", format(nrow(raw), big.mark = ","),
    n_distinct(raw$label))
say("Rows with no damage fit attempted (nalign = -1): %d", sum(raw$nalign == -1))

## GRAIN. Everything downstream assumes one row per sample per taxon.
stopifnot(!any(duplicated(paste(raw$label, raw$taxid))))
say("Grain check passed: (label, taxid) is unique.")

## Three arithmetic checks on real rows. If the columns have shifted by one --
## the classic failure mode with a 187-column file -- these go wrong at once and
## everything downstream is confident nonsense.
chk <- raw |> filter(nreads >= MIN_READS_FIT) |> slice_head(n = 500)
fwf_cols  <- str_c("fwf",  0:(NPOS - 1))
fwdx_cols <- str_c("fwdx", 0:(NPOS - 1))
fwK_cols  <- str_c("fwK",  0:(NPOS - 1))
fwN_cols  <- str_c("fwN",  0:(NPOS - 1))

err_model <- max(abs(outer(chk$A, 0:(NPOS - 1), \(a, x) a) *
                       outer(1 - chk$q, 0:(NPOS - 1), `^`) + chk$c -
                       as.matrix(chk[, fwdx_cols])), na.rm = TRUE)
err_freq  <- max(abs(as.matrix(chk[, fwK_cols]) / as.matrix(chk[, fwN_cols]) -
                       as.matrix(chk[, fwf_cols])), na.rm = TRUE)
err_zfit  <- max(abs(chk$A / chk$sigmaD_new - chk$Zfit_new), na.rm = TRUE)

say("CHECK  D(x) = A(1-q)^x + c  vs fwdx : %.2e", err_model)
say("CHECK  fwK/fwN              vs fwf  : %.2e", err_freq)
say("CHECK  A/sigmaD_new         vs Zfit : %.2e", err_zfit)
if (max(err_model, err_freq) > 1e-4)
  stop("Model identities do not hold -- the columns are misaligned. Stop here.")


# 2. Damage ---------------------------------------------------------------------

## 2.1 check functions 
stopifnot(file.exists(PATH_DMGR))
source(PATH_DMGR) 


## 2.2 make subset with 100 reads
fit_in <- raw |>
  filter(nreads >= MIN_READS_FIT) |>
  select(tax_name = taxid, label, all_of(c(fwf_cols, fwdx_cols)))

say("Fitting damage quality for %s of %s rows (nreads >= %d)",
    format(nrow(fit_in), big.mark = ","), format(nrow(raw), big.mark = ","),
    MIN_READS_FIT)

## 2.3 estimating quality of damage fit 
set.seed(DMG_SEED)
ccc <- dmg_fwd_CCC(fit_in, unique(fit_in$label),
                   ci = "asymptotic", nperm = NPERM, nproc = NPROC) |>
  ungroup() |>
  rename(taxid = tax_name)

write_tsv(ccc, "cleandata/MTG_euk_ccc.tsv.gz")

## 2.4 join fit data to the table 
dmg <- raw |>
  left_join(ccc, by = c("label", "taxid")) |>
  mutate(
    fit = case_when(
      is.na(rho_c)                            ~ NA_character_,
      q_CI_h >= 1 | c_CI_l <= 0               ~ "bad",   # boundary veto
      rho_c >= CCC_RHO & C_b > CCC_CB &
        round(rho_c_perm_pval, 3) < CCC_P     ~ "good",
      TRUE                                    ~ "bad"
    ),
    ## A_b is the bootstrap amplitude. We keep it NA when the fit failed, so
    ## "no damage" and "cannot tell" stay distinguishable -- the bacterial
    ## script wrote 0 for both, which merged them.
    damage = if_else(fit == "good", A_b, NA_real_),
    Zfit   = Zfit_new    # metaDMG now ships _old and _new; _new is primary
  )


## 2.5 Report outputs

say("FIT OUTCOME (NA = fewer than %d reads, never tested)", MIN_READS_FIT)
print(table(dmg$fit, useNA = "ifany"))
## Counted among the rows that were actually TESTED -- metaDMG reports these CIs
## for every taxon, including the thousands never fitted, so counting the whole
## table would report a veto rate of almost 100%.
tested <- dmg |> filter(!is.na(rho_c))
say("Of %d tested rows, vetoed by q_CI_h >= 1: %d | by c_CI_l <= 0: %d",
    nrow(tested), sum(is_true(tested$q_CI_h >= 1)), sum(is_true(tested$c_CI_l <= 0)))

## Is "ancient" just "abundant"? If the pass rate climbs steeply with read count
say("P(fit == 'good') by read-count bin -- watch for a steep climb:")
dmg |>
  filter(!is.na(fit)) |>
  mutate(bin = cut(nreads, c(100, 300, 1e3, 3e3, 1e4, Inf),
                   include.lowest = TRUE, dig.lab = 6)) |>
  count(bin, fit) |>
  pivot_wider(names_from = fit, values_from = n, values_fill = 0) |>
  mutate(p_good = good / (good + bad)) |>
  print()


# 3. Taxonomy ---------------------------------------------------------------------

## One row per (taxid, lineage element). Parse the distinct paths only -- a path
## is a property of the taxid, not of the sample.
tax_long <- dmg |>
  distinct(taxid, taxa_path) |>
  ## quotes stripped again here
  mutate(element = str_split(str_remove_all(taxa_path, '"'), ";")) |>
  unnest_longer(element, indices_to = "depth") |>      # depth 1 = the taxon itself
  separate_wider_regex(
    element,
    c(el_taxid = "[0-9]+", ":", el_name = ".*", ":", el_rank = "[^:]*")
  )

say("Parsed %s taxa into %s lineage elements",
    format(n_distinct(tax_long$taxid), big.mark = ","),
    format(nrow(tax_long), big.mark = ","))

say("Ranks appearing more than once in a single path (kept vs dropped):")
tax_long |>
  filter(el_rank %in% RANKS) |>
  arrange(taxid, el_rank, desc(depth)) |>
  summarise(kept = first(el_name), dropped = str_c(el_name[-1], collapse = " ; "),
            .by = c(taxid, el_rank)) |>
  filter(dropped != "") |>
  distinct(el_rank, kept, dropped) |>
  slice_head(n = 3, by = el_rank) |>
  print(n = Inf)

## Deepest element = nearest the root = the canonical one.
## Start from EVERY taxid, then attach the ranks. Starting from the pivot would
## silently drop root, "cellular organisms" and the viral "acellular root" -- they
## carry none of RANKS -- and their subset flags would come out NA rather than FALSE.
tax_wide <- tax_long |>
  distinct(taxid) |>
  left_join(
    tax_long |>
      filter(el_rank %in% RANKS) |>
      slice_max(depth, by = c(taxid, el_rank), n = 1, with_ties = FALSE) |>
      select(taxid, el_rank, el_name) |>
      pivot_wider(names_from = el_rank, values_from = el_name),
    by = "taxid"
  )

## Any rank absent from every path still needs its column.
for (r in setdiff(RANKS, names(tax_wide))) tax_wide[[r]] <- NA_character_
tax_wide <- tax_wide |> select(taxid, all_of(RANKS))

## Subset membership, by taxid on the path rather than by name.
for (nm in names(SUBSETS)) {
  members <- unique(tax_long$taxid[tax_long$el_taxid == as.character(SUBSETS[[nm]])])
  tax_wide[[str_c("in_", nm)]] <- tax_wide$taxid %in% members
}

## The ladder being entirely NA is the failure mode this section is prone to,
## and it is invisible until a barplot comes out empty. Catch it here.
filled <- tax_wide |> summarise(across(all_of(RANKS), ~ sum(!is.na(.x))))
say("Taxa with each rank filled in:")
print(filled)
if (all(unlist(filled) == 0)) {
  say("Ranks actually seen in taxa_path were:")
  print(sort(table(tax_long$el_rank), decreasing = TRUE)[1:10])
  stop("No lineage element matched RANKS. If the ranks above are wrapped in ",
       "quote characters the strip failed; otherwise RANKS needs updating.")
}

say("Subset sizes (distinct taxa):")
print(tax_wide |> summarise(across(starts_with("in_"), sum)))

say("Spot-check -- compare against the raw paths above:")
tax_wide |> filter(!is.na(species)) |> slice_head(n = 5) |>
  select(taxid, kingdom, phylum, class, order, family, genus, species) |> print()

dmg <- dmg |> left_join(tax_wide, by = "taxid")
stopifnot(nrow(dmg) == nrow(raw))            # a join must never multiply rows

say("Rows by kingdom (NA = lineages NCBI gives no kingdom, e.g. the SAR clade):")
print(dmg |> count(kingdom, sort = TRUE))

# 4. Metadata ---------------------------------------------------------------------

normalise_label <- function(x) if_else(str_starts(x, "WP2"), str_replace_all(x, "-", "."), x)

meta_raw <- read_excel(PATH_META)
say("metadata.xlsx: %d rows x %d cols (blank rows are dropped next)",
    nrow(meta_raw), ncol(meta_raw))

samples <- meta_raw |>
  filter(!is.na(MTG.ID)) |>
  transmute(
    label        = normalise_label(MTG.ID),
    lipid_id     = LipidID,
    core         = CoreID,
    batch        = factor(SeqExp),
    depth_raw    = as.character(Length),
    age_raw      = as.character(Est.Cal.YrBP),
    section      = Section,
    ## a few of these carry text in some rows; the coercion to NA is intended
    dup_rate     = suppressWarnings(as.numeric(`duplication rate`)),
    gc_lib       = suppressWarnings(as.numeric(`GC content`)),
    insert_peak  = suppressWarnings(as.numeric(`Insert Size Peak`)),
    reads_premap = suppressWarnings(as.numeric(`premapping reads`))
  ) |>
  mutate(
    is_modern = !is.na(age_raw) & str_to_lower(str_trim(age_raw)) == "modern",
    is_top    = !is.na(depth_raw) & str_to_upper(str_trim(depth_raw)) == "TOP",
    age_yrBP  = if_else(is_modern, MODERN_AGE_YRBP, suppressWarnings(as.numeric(age_raw))),
    depth_cm  = if_else(is_top, 0, suppressWarnings(as.numeric(depth_raw))),
    is_control = is.na(core),                   # controls are the rows with no core
    ctrl_type = case_when(
      !is_control             ~ NA_character_,
      str_detect(label, "LibNTC") ~ "library_blank",
      str_detect(label, "ExrNTC") ~ "extraction_blank",
      str_detect(label, "LibPTC") ~ "library_positive",
      str_detect(label, "ExrPTC") ~ "extraction_positive",
      TRUE                        ~ "other_control"
    )
  )

say("Usable metadata rows: %d (%d samples, %d controls)",
    nrow(samples), sum(!samples$is_control), sum(samples$is_control))
say("Parsed %d 'modern' ages and %d 'TOP' depths",
    sum(samples$is_modern), sum(samples$is_top))

bad_age <- samples |> filter(!is_control, is.na(age_yrBP))
if (nrow(bad_age)) {
  say("WARNING: %d non-control rows have an unparseable age:", nrow(bad_age))
  print(bad_age |> select(label, core, age_raw))
}

## Core coordinates. The file has a UTF-8 BOM and CRLF endings; read_csv strips
## both. The rename_with(str_trim) matters though -- the header is `lat ` with a
## trailing space, and without trimming `lat` simply does not exist.
sites <- read_csv(PATH_SITE, show_col_types = FALSE) |>
  rename_with(str_trim) |>
  transmute(core = CoreID, site_name = name, lat = as.numeric(lat), lon = as.numeric(lon))
say("Sites: %s", str_c(sites$core, collapse = ", "))

samples <- samples |> left_join(sites, by = "core")

## --- label reconciliation ---------------------------------------------------
lab_dat  <- unique(dmg$label)
lab_meta <- unique(samples$label)
say("Labels in data: %d | in metadata: %d | matched: %d",
    length(lab_dat), length(lab_meta), length(intersect(lab_dat, lab_meta)))

miss_in_meta <- setdiff(lab_dat, lab_meta)
miss_in_dat  <- setdiff(lab_meta, lab_dat)

if (length(miss_in_dat)) {
  say("In metadata but with no data rows (%d):", length(miss_in_dat))
  print(samples |> filter(label %in% miss_in_dat) |>
          select(label, core, ctrl_type, reads_premap) |> arrange(label), n = Inf)
}
if (length(miss_in_meta)) {
  say("!! IN DATA BUT NOT IN METADATA (%d):", length(miss_in_meta))
  print(miss_in_meta)
  stop("Unmatched labels. Fix normalise_label() before continuing.")
}

dmg <- dmg |> left_join(samples, by = "label")
stopifnot(nrow(dmg) == nrow(raw))


# 5. Master table ---------------------------------------------------------------

## 5.1 Build columns ----
## the 150 per-position columns (fwK0..bwdxConf14), left out of the .tsv.gz
pos_cols <- names(dmg) |> str_subset("^(fw|bw)(K|N|f|dx|dxConf)[0-9]+$")

master <- dmg |>
  mutate(
    ## damage is already NA unless the fit is good, so the fit clause is a
    ## safeguard in case `damage` is ever built differently
    is_ancient = if (REQUIRE_GOOD_FIT) is_true(fit == "good" & damage >= DAMAGE_THRESHOLD)
    else                  is_true(damage >= DAMAGE_THRESHOLD)
  ) |>
  select(
    ## who and where
    label, lipid_id, core, site_name, lat, lon, batch,
    depth_cm, age_yrBP, is_modern, is_control, ctrl_type, reads_premap,
    ## what
    taxid, name, rank, all_of(RANKS), starts_with("in_"),
    ## how much
    nalign, n_reads = nreads, mean_rlen, var_rlen, mean_gc, var_gc,
    mean_dust, var_dust, mean_nspec, var_nspec,
    ## damage
    A, q, c, phi, A_b, q_b, c_b, phi_b, Zfit, Zfit_old,
    A_CI_l, A_CI_h, q_CI_l, q_CI_h, c_CI_l, c_CI_h,
    rho_c, C_b, rho_c_perm_pval, l_shift, s_shift,
    fit, damage, is_ancient,
    ## provenance
    lca, taxa_path,
    everything()
  )

say("MASTER TABLE: %s rows x %d columns",
    format(nrow(master), big.mark = ","), ncol(master))
stopifnot(nrow(master) == nrow(raw))

## 5.2 Summary ----
say("Contents by subset and control status:")
print(master |> count(is_control, in_metazoa, in_eukaryota))

say("Reads per sample (all taxa, note these are CUMULATIVE up the tree):")
print(master |> summarise(taxa = n(), reads_at_root = sum(n_reads[taxid == 1]),
                          .by = c(label, core)) |> arrange(core, label), n = 108)



# 6. Decontamination ------------------------------------------------------------
## 6.2 flags sample rows with fewer reads than that taxon's blank threshold;
## 6.3 removes them, and every control, to make `clean` -- the analysis dataset.
## `dat` keeps every row for the control summaries here and the read counts in 7.
## Everything about the negative controls is written to cleandata/controls/;
## positive controls are ignored throughout.

dir.create(PATH_CTRL, showWarnings = FALSE, recursive = TRUE)

## 6.1 Negative controls ----
dat <- master

neg_labels <- if (is.null(NEG_CONTROLS)) {
  samples$label[is_true(samples$ctrl_type %in% c("library_blank", "extraction_blank"))]
} else NEG_CONTROLS
neg_labels <- intersect(neg_labels, unique(dat$label))

say("Negative controls in the data (%d):", length(neg_labels))
print(dat |> filter(label %in% neg_labels) |>
        summarise(n_taxa = n(), reads = sum(.data[[DECONTAM_COL]]),
                  biggest_taxon = max(.data[[DECONTAM_COL]]), .by = label))

## 6.2 Thresholds and flags ----
if (!DECONTAM || !length(neg_labels)) {
  if (!length(neg_labels)) warning("No negative controls found -- decontamination skipped.")
  say("DECONTAM off: nothing flagged.")
  dat <- dat |> mutate(blank_thr = 0, decontam_removed = FALSE)
} else {
  
  ## one threshold per taxon, over the blanks the taxon appears in
  ## !! FLAG: with BLANK_STAT = "mean", blanks where the taxon is absent are not
  ## !!       counted as zeros, so "mean" is higher than the true mean over all blanks
  blank_tbl <- dat |>
    filter(label %in% neg_labels) |>
    summarise(blank_thr = if (BLANK_STAT == "mean") mean(.data[[DECONTAM_COL]])
              else                      max(.data[[DECONTAM_COL]]),
              n_blanks_seen = n_distinct(label), .by = taxid) |>
    mutate(blank_thr = blank_thr * BLANK_MULTIPLIER)
  
  say("%d distinct taxa appear in the blanks (threshold = %s x %g)",
      nrow(blank_tbl), BLANK_STAT, BLANK_MULTIPLIER)
  say("Blank taxa by kingdom:")
  print(dat |> filter(taxid %in% blank_tbl$taxid) |> distinct(taxid, kingdom) |>
          count(kingdom, sort = TRUE))
  
  dat <- dat |>
    left_join(blank_tbl |> select(taxid, blank_thr), by = "taxid") |>
    mutate(blank_thr = replace_na(blank_thr, 0),
           decontam_removed = !is_control & .data[[DECONTAM_COL]] < blank_thr)
  
  ## -- what was flagged --
  removed <- dat |> filter(decontam_removed)
  n_samp_obs <- sum(!dat$is_control)
  
  say("FLAGGED %s of %s sample observations (%.2f%%)",
      format(nrow(removed), big.mark = ","), format(n_samp_obs, big.mark = ","),
      100 * nrow(removed) / n_samp_obs)
  say("FLAGGED %s of %s reads (%.2f%%)",
      format(sum(removed[[DECONTAM_COL]]), big.mark = ","),
      format(sum(dat[[DECONTAM_COL]][!dat$is_control]), big.mark = ","),
      100 * sum(removed[[DECONTAM_COL]]) / sum(dat[[DECONTAM_COL]][!dat$is_control]))
  say("Ancient observations flagged: %d of %d",
      sum(removed$is_ancient), sum(dat$is_ancient & !dat$is_control))
  
  ## -- cost per sample: a sample losing most of its taxa is mostly reagent --
  cost <- dat |>
    filter(!is_control) |>
    summarise(before = n(), lost = sum(decontam_removed), .by = label) |>
    mutate(pct = 100 * lost / before) |>
    arrange(desc(pct))
  say("Per-sample taxa flagged: median %.0f (%.2f%%), worst %d (%.2f%%)",
      median(cost$lost), median(cost$pct), max(cost$lost), max(cost$pct))
  print(head(cost, 5))
  
  ## -- fate of each blank taxon across the samples --
  fate <- dat |>
    filter(!is_control, taxid %in% blank_tbl$taxid) |>
    summarise(n_kept = sum(!decontam_removed),
              med_reads_kept = median(.data[[DECONTAM_COL]][!decontam_removed]),
              blank_thr = first(blank_thr), .by = c(taxid, name)) |>
    mutate(margin = med_reads_kept / blank_thr) |>
    arrange(desc(n_kept))
  say("Blank taxa eliminated everywhere: %d ; surviving somewhere: %d",
      sum(fate$n_kept == 0), sum(fate$n_kept > 0))
  print(fate |> filter(n_kept > 0) |> head(10))
  
  write_tsv(blank_tbl |>
              left_join(dat |> distinct(taxid, name, rank, kingdom), by = "taxid") |>
              relocate(name, rank, kingdom, .after = taxid) |>
              arrange(desc(blank_thr)),
            file.path(PATH_CTRL, "MTG_euk_decontam_blank_thresholds.tsv"))
  write_tsv(removed |> select(label, taxid, name, rank, kingdom, genus,
                              all_of(DECONTAM_COL), blank_thr),
            file.path(PATH_CTRL, "MTG_euk_decontam_flagged.tsv"))
  write_tsv(fate, file.path(PATH_CTRL, "MTG_euk_decontam_blank_fate.tsv"))
}

stopifnot(nrow(dat) == nrow(raw))

## 6.3 Remove contaminants and controls ----
## the analysis dataset: real samples only, flagged rows gone. The three columns
## dropped would be constant here. Removed rows are in controls/ (6.2, 6.4).
clean <- dat |>
  filter(!is_control, !decontam_removed) |>
  select(-is_control, -ctrl_type, -decontam_removed)

n_ctrl_rows <- sum(dat$is_control)
n_flag_rows <- sum(dat$decontam_removed)
stopifnot(nrow(clean) == nrow(dat) - n_ctrl_rows - n_flag_rows)

say("CLEAN DATASET: %s rows from %d samples (removed %s control rows and %s contaminant rows)",
    format(nrow(clean), big.mark = ","), n_distinct(clean$label),
    format(n_ctrl_rows, big.mark = ","), format(n_flag_rows, big.mark = ","))

## 6.4 Control data ----
## every negative-control row, without the per-position damage columns
ctrl_dat <- dat |> filter(is_control, ctrl_type %in% REPORT_CTRL)
write_tsv(ctrl_dat |> select(-all_of(pos_cols)),
          file.path(PATH_CTRL, "MTG_euk_controls.tsv.gz"))
say("Control data: %s rows from %d controls",
    format(nrow(ctrl_dat), big.mark = ","), n_distinct(ctrl_dat$label))


## 6.5 Taxonomic summary of the controls ----
## n_reads is cumulative, so summing rows double-counts. Instead work out each
## taxon's own reads -- n_reads minus its direct children's n_reads, i.e. reads
## whose LCA is exactly that taxon -- and group those by lineage: every read is
## then counted once. A taxon's parent is the second entry of its taxa_path.
ctrl_info <- samples |>
  filter(is_control, ctrl_type %in% REPORT_CTRL) |>
  transmute(label, ctrl_type, batch, reads_premap, in_data = label %in% dat$label)

ctrl_own <- ctrl_dat |>
  mutate(parent = as.numeric(str_match(taxa_path, "^[^;]*;([0-9]+):")[, 2]))
kids <- ctrl_own |>
  filter(!is.na(parent)) |>
  summarise(kid_reads = sum(n_reads), .by = c(label, parent))
ctrl_own <- ctrl_own |>
  left_join(kids, by = c("label", "taxid" = "parent")) |>
  mutate(own = n_reads - replace_na(kid_reads, 0)) |>
  select(label, taxid, domain, kingdom, phylum, genus, in_eukaryota, in_metazoa, own)

## check: no negative own reads, and own reads add up to each control's root count
chk <- ctrl_own |>
  summarise(own = sum(own), .by = label) |>
  left_join(ctrl_dat |> filter(taxid == 1) |> select(label, root = n_reads), by = "label")
stopifnot(all(ctrl_own$own >= 0), isTRUE(all(chk$own == chk$root)))
ctrl_own <- ctrl_own |> filter(own > 0)

## -- (a) one row per control --
## controls with no data rows were sequenced and produced nothing: 0, not NA
top_genus <- ctrl_own |>
  filter(!is.na(genus)) |>
  summarise(reads = sum(own), .by = c(label, genus)) |>
  slice_max(reads, n = 1, by = label, with_ties = FALSE) |>
  rename(top_genus = genus, top_genus_reads = reads)

ctrl_overview <- ctrl_info |>
  left_join(ctrl_own |>
              summarise(reads_assigned  = sum(own),
                        reads_eukaryota = sum(own[in_eukaryota]),
                        reads_metazoa   = sum(own[in_metazoa]),
                        n_genera        = n_distinct(genus, na.rm = TRUE),
                        .by = label),
            by = "label") |>
  mutate(across(c(reads_assigned, reads_eukaryota, reads_metazoa, n_genera),
                ~ replace_na(.x, 0))) |>
  left_join(top_genus, by = "label") |>
  mutate(pct_assigned  = 100 * reads_assigned / reads_premap,
         pct_eukaryota = if_else(reads_assigned > 0, 100 * reads_eukaryota / reads_assigned, NA_real_),
         pct_metazoa   = if_else(reads_assigned > 0, 100 * reads_metazoa / reads_assigned, NA_real_)) |>
  arrange(ctrl_type, label)

say("CONTROLS OVERVIEW (reads with an LCA, before any filtering):")
print(ctrl_overview |> select(label, ctrl_type, reads_premap, reads_assigned,
                              pct_eukaryota, pct_metazoa, top_genus, top_genus_reads), n = Inf)

## -- (b) reads by domain / kingdom / phylum, per control --
ctrl_taxa <- ctrl_own |>
  mutate(across(c(domain, kingdom, phylum), ~ replace_na(.x, "unassigned"))) |>
  summarise(reads = sum(own), .by = c(label, domain, kingdom, phylum)) |>
  mutate(pct = 100 * reads / sum(reads), .by = label) |>
  left_join(ctrl_info |> select(label, ctrl_type), by = "label") |>
  relocate(ctrl_type, .after = label) |>
  arrange(ctrl_type, label, desc(reads))

## -- (c) genera per control type --
## mean_reads averages over ALL controls of that type, including empty ones
n_type <- ctrl_overview |> count(ctrl_type, name = "n_controls")

ctrl_genera <- ctrl_own |>
  filter(!is.na(genus)) |>
  summarise(reads = sum(own), .by = c(label, kingdom, genus)) |>
  left_join(ctrl_info |> select(label, ctrl_type), by = "label") |>
  summarise(n_detected = n_distinct(label), total_reads = sum(reads),
            max_reads = max(reads), .by = c(ctrl_type, kingdom, genus)) |>
  left_join(n_type, by = "ctrl_type") |>
  mutate(mean_reads  = total_reads / n_controls,
         pct_of_type = 100 * total_reads / sum(total_reads), .by = ctrl_type) |>
  relocate(n_controls, .after = n_detected) |>
  arrange(ctrl_type, desc(total_reads))

say("TOP 10 GENERA PER CONTROL TYPE:")
print(ctrl_genera |> slice_head(n = 10, by = ctrl_type), n = Inf)

say("Kingdoms in the negative controls (%% of their assigned reads):")
print(ctrl_taxa |>
        filter(ctrl_type %in% c("library_blank", "extraction_blank")) |>
        summarise(reads = sum(reads), .by = kingdom) |>
        mutate(pct = round(100 * reads / sum(reads), 2)) |>
        arrange(desc(reads)), n = 15)

write_tsv(ctrl_overview, file.path(PATH_CTRL, "MTG_euk_controls_overview.tsv"))
write_tsv(ctrl_taxa,     file.path(PATH_CTRL, "MTG_euk_controls_taxa_phylum.tsv"))
write_tsv(ctrl_genera,   file.path(PATH_CTRL, "MTG_euk_controls_genera.tsv"))
say("Control files written to %s/: %s", PATH_CTRL,
    str_c(list.files(PATH_CTRL), collapse = ", "))


# 7. Read summary ---------------------------------------------------------------
## one row per library in metadata.xlsx (samples and negative controls), then group
## means and SDs for the supplement. Assigned = reads with an LCA (species level
## or finer in this run), counted before decontamination.

## 7.1 Per library ----
## fastP writes counts as text with a unit ("680.336862 M", "5.191000 K");
## plain numbers pass through unchanged
parse_count <- function(x){
  x    <- str_trim(as.character(x))
  unit <- str_extract(x, "[KMG]$")
  num  <- suppressWarnings(as.numeric(str_remove(x, "\\s*[KMG]$")))
  round(num * ifelse(is.na(unit), 1, c(K = 1e3, M = 1e6, G = 1e9)[unit]))
}

meta_reads <- function(col){
  if (is.null(col)) return(rep(NA_real_, nrow(samples)))
  if (!col %in% names(meta_raw))
    stop("Column '", col, "' is not in metadata.xlsx. Check the META_COL_* settings.")
  m <- meta_raw |> filter(!is.na(MTG.ID))
  x <- as.character(m[[col]])
  x[str_trim(x) %in% c("", "NA", "NaN")] <- NA      # the sheet has some literal "NA" text
  v <- parse_count(x)
  bad <- !is.na(x) & is.na(v)
  if (any(bad)) stop("Could not read ", sum(bad), " values in '", col, "', e.g. '",
                     x[bad][1], "'")
  if (any(is.na(v))) say("NOTE: '%s' is missing for: %s", col,
                         str_c(normalise_label(m$MTG.ID)[is.na(v)], collapse = ", "))
  v[match(samples$label, normalise_label(m$MTG.ID))]
}
node_reads <- function(id){
  x <- dat |> filter(taxid == id)
  x$n_reads[match(samples$label, x$label)]
}

read_sum <- samples |>
  transmute(label, core, is_control, ctrl_type, batch, depth_cm, age_yrBP,
            group = if_else(is_control, ctrl_type, core),
            in_data = label %in% dat$label) |>
  mutate(reads_raw       = meta_reads(META_COL_RAW),
         reads_filtered  = meta_reads(META_COL_FILTERED),
         reads_premap    = samples$reads_premap,
         reads_mapped    = meta_reads(META_COL_MAPPED),
         reads_assigned  = node_reads(1),
         reads_eukaryota = node_reads(2759),
         reads_metazoa   = node_reads(33208)) |>
  ## a library in the data with no Metazoa row has 0 metazoan reads; a library
  ## that was sequenced (has pre-mapping reads) but has no data rows assigned 0
  mutate(across(c(reads_assigned, reads_eukaryota, reads_metazoa),
                ~ if_else(is.na(.x) & (in_data | !is.na(reads_premap)), 0, .x)),
         ## percentages are NA, not NaN, when the denominator is 0
         pct_premap    = if_else(reads_raw > 0,    100 * reads_premap / reads_raw, NA_real_),       # of raw
         pct_mapped    = if_else(reads_premap > 0, 100 * reads_mapped / reads_premap, NA_real_),    # of pre-mapping
         pct_assigned  = if_else(reads_premap > 0, 100 * reads_assigned / reads_premap, NA_real_),  # of pre-mapping
         pct_eukaryota = if_else(reads_assigned > 0, 100 * reads_eukaryota / reads_assigned, NA_real_),
         pct_metazoa   = if_else(reads_assigned > 0, 100 * reads_metazoa / reads_assigned, NA_real_)) |>
  arrange(is_control, group, age_yrBP, label) |>
  filter(!is_control | ctrl_type %in% REPORT_CTRL)               # drop positives

for (v in c("META_COL_RAW", "META_COL_FILTERED", "META_COL_MAPPED"))
  if (is.null(get(v))) say("NOTE: %s not set -- that column is NA and left out of the supplement.", v)
say("Libraries summarised: %d (%d samples, %d controls; %d controls with no data rows)",
    nrow(read_sum), sum(!read_sum$is_control), sum(read_sum$is_control),
    sum(read_sum$is_control & !read_sum$in_data))

write_csv(read_sum, "cleandata/MTG_euk_read_summary_libraries.csv")

## 7.2 Group means and SDs ----
read_cols   <- c("reads_raw", "reads_filtered", "reads_premap", "reads_mapped",
                 "reads_assigned", "reads_eukaryota", "reads_metazoa",
                 "pct_premap", "pct_mapped", "pct_assigned", "pct_eukaryota", "pct_metazoa")
GROUP_ORDER <- c("All samples", CORE_ORDER_PLOT, "extraction_blank", "library_blank")

grp <- bind_rows(read_sum |> filter(!is_control) |> mutate(group = "All samples"),
                 read_sum) |>
  filter(!is.na(group)) |>
  summarise(n = n(),
            across(all_of(read_cols),
                   list(mean = ~ mean(.x, na.rm = TRUE), sd = ~ sd(.x, na.rm = TRUE))),
            .by = group) |>
  mutate(across(where(is.numeric), ~ replace(.x, is.nan(.x), NA)),
         group = factor(group, levels = union(GROUP_ORDER, group))) |>
  arrange(group)

write_csv(grp, "cleandata/MTG_euk_read_summary_groups.csv")

## 7.3 Supplement table ----
## "mean ± SD" per group; SD is blank for groups of one, columns that are all
## NA (e.g. raw reads not set) are dropped. Written with a BOM so Excel shows ±.
SUP_COLS <- c(reads_raw = "Raw reads", reads_filtered = "fastp-filtered reads",
              reads_premap = "Pre-mapping reads",
              reads_mapped = "Mapped reads", reads_assigned = "Assigned reads",
              reads_eukaryota = "Eukaryotic reads", reads_metazoa = "Metazoan reads",
              pct_premap = "% raw reads to mapping",
              pct_mapped = "% mapped", pct_assigned = "% assigned",
              pct_eukaryota = "% eukaryotic", pct_metazoa = "% metazoan")
GROUP_LABELS <- c(extraction_blank = "Extraction blanks", library_blank = "Library blanks")

fmt_msd <- function(m, s, digits){
  f <- function(x) formatC(x, format = "f", digits = digits, big.mark = ",")
  ifelse(is.na(m), "-", ifelse(is.na(s), f(m), paste0(f(m), " ± ", f(s))))
}

grp_chr <- as.character(grp$group)
sup <- tibble(Group = unname(ifelse(grp_chr %in% names(GROUP_LABELS),
                                    GROUP_LABELS[grp_chr], grp_chr)),
              n = grp$n)
for (v in names(SUP_COLS)) {
  m <- grp[[str_c(v, "_mean")]]
  if (all(is.na(m))) next
  sup[[SUP_COLS[[v]]]] <- fmt_msd(m, grp[[str_c(v, "_sd")]],
                                  digits = if (str_starts(v, "pct")) 2 else 0)
}

say("SUPPLEMENTARY READ SUMMARY (mean ± SD per library):")
print(as.data.frame(sup), row.names = FALSE)
say("Footnote: %% raw reads to mapping = pre-mapping / raw reads;")
say("%% mapped and %% assigned are of pre-mapping reads; %% eukaryotic and")
say("%% metazoan are of assigned reads. Assigned = reads given an LCA at species level")
say("or finer, before decontamination.")

write_excel_csv(sup, "cleandata/MTG_euk_read_summary_supplement.csv")


# 8. Write out ------------------------------------------------------------------
## the cleaned analysis dataset from 6.3: real samples only, contaminant rows
## removed. The analysis script starts from MTG_euk_master.rds and needs no
## further decontamination filtering. Controls are in cleandata/controls/.
## The .tsv.gz drops the per-position damage columns; the .rds keeps everything.

stopifnot(!any(clean$label %in% samples$label[samples$is_control]))

saveRDS(clean, "cleandata/MTG_euk_master.rds")
saveRDS(samples |> filter(!is_control), "cleandata/MTG_euk_samples.rds")

samp_cols <- intersect(names(samples), names(clean))
tax_cols  <- c("taxid", "name", "rank", RANKS, str_c("in_", names(SUBSETS)), "lca", "taxa_path", "blank_thr")
gh_samples <- clean |> select(all_of(samp_cols)) |> distinct()
gh_taxa    <- clean |> select(all_of(tax_cols)) |> distinct()
stopifnot(!anyDuplicated(gh_samples$label), !anyDuplicated(gh_taxa$taxid))   # truly one row each

write_tsv(gh_samples, "cleandata/MTG_euk_master_samples.tsv.gz")
write_tsv(gh_taxa,    "cleandata/MTG_euk_master_taxa.tsv.gz")
write_tsv(clean |> select(label, taxid, !any_of(c(samp_cols, tax_cols, pos_cols))),
          "cleandata/MTG_euk_master_obs.tsv.gz")

say("DONE. Analysis dataset: %s rows x %d cols, %d samples. GitHub files (MB): %s",
    format(nrow(clean), big.mark = ","), ncol(clean), n_distinct(clean$label),
    str_c(sprintf("%s %.1f", c("obs", "samples", "taxa"),
                  file.size(str_c("cleandata/MTG_euk_master_", c("obs", "samples", "taxa"), ".tsv.gz")) / 1e6),
          collapse = ", "))




