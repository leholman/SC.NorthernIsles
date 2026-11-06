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
DECONTAM_COL   <- "n_reads"  # the bacterial script used n_reads_ref, which does
BLANK_STAT     <- "mean"      # "max" = beat the worst blank; "mean" = the average
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


## 2.1 make subset with 100 reads
fit_in <- raw |>
  filter(nreads >= MIN_READS_FIT) |>
  select(tax_name = taxid, label, all_of(c(fwf_cols, fwdx_cols)))

say("Fitting damage quality for %s of %s rows (nreads >= %d)",
    format(nrow(fit_in), big.mark = ","), format(nrow(raw), big.mark = ","),
    MIN_READS_FIT)

## 2.2 estimating quality of damage fit 
set.seed(DMG_SEED)
ccc <- dmg_fwd_CCC(fit_in, unique(fit_in$label),
                   ci = "asymptotic", nperm = NPERM, nproc = NPROC) |>
  ungroup() |>
  rename(taxid = tax_name)

write_tsv(ccc, "cleandata/MTG_euk_ccc.tsv.gz")

## 2.3 join fit data to the table 
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


## 2.4 Report outputs

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



