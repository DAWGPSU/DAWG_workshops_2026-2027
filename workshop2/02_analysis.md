# Part 2: Analysing bacterial and resistance-gene profiles in R

**Alpha diversity, beta diversity and differential abundance of gut bacteria (MetaPhlAn 4) and antimicrobial resistance genes (RGI bwt), in a longitudinal study**

In [Part 1](01_run_profiling.md) raw reads were turned into one table per tool and per sample. Here we load those tables into R and ask, for the **bacterial community** and for the **resistance genes**:

1. Is the data quality adequate? (read QC)
2. How many species or resistance genes does each sample contain? (**alpha diversity / richness**)
3. How different are samples from each other, and how far does each person drift from their own baseline? (**beta diversity**, bacteria only)
4. Which individual species, resistance genes or drug classes change over time? (**differential abundance**)

> [!NOTE]
> **We are not running raw data through R in this workshop.** The raw reads were processed with the Part 1 pipeline beforehand, and we give you the resulting tables inside a **workspace image** (an `.RData` file). You load the image, and everything you need is already in memory. The code that read those tables is still printed below (marked **Read only**) so you can see how the data got here, but you will not run it.

## How to read this document

| Marker | Meaning |
|---|---|
| **IMPORTANT: Run this** (blue-purple box) | Part of the hands-on session. Run it. |
| **NOTE: Read only** (blue box) | Shown for transparency. It is skipped automatically when you load the workspace image. |
| **TIP** (green box) | Something to look at in the output. |
| **WARNING** (orange box) | A common mistake. |

---

## Contents

1. [The study](#1-the-study)
2. [Set-up and the workspace image](#2-set-up-and-the-workspace-image)
3. [Statistical design at a glance](#3-statistical-design-at-a-glance)
4. [Read QC](#4-read-qc)
5. [Bacteria](#5-bacteria)
6. [Antimicrobial resistance genes](#6-antimicrobial-resistance-genes)
7. [Save the workspace image](#7-save-the-workspace-image)
8. [Pitfalls](#8-pitfalls)
9. [References](#9-references)

---

## 1. The study

We use the dataset from **Palleja et al. 2018, *Nature Microbiology* 3, 1255-1265** ([doi:10.1038/s41564-018-0257-9](https://doi.org/10.1038/s41564-018-0257-9)), "Recovery of gut microbiota of healthy adults following antibiotic exposure". The reads are public (European Nucleotide Archive project PRJEB20800).

- **Who:** 12 healthy young men.
- **What was done:** a 4-day course of three broad-spectrum antibiotics (meropenem, gentamicin and vancomycin), chosen to mimic the preventive treatment sometimes given in intensive care.
- **What was measured:** stool samples at baseline (before treatment) and at four later timepoints over about six months, sequenced by shotgun metagenomics.
- **What the authors found:** bacterial diversity fell sharply, then recovered to close to baseline within about 1.5 months. Nine common species were still missing in most people at six months. Bacteria carrying beta-lactam resistance genes were selected for during and after treatment.

This is a **longitudinal** design: the *same* 12 people are sampled repeatedly. That is why the statistics below always account for "which person is this sample from". It also means there is no untreated control group, so we describe change **relative to each person's own baseline** rather than comparing treatment arms.

> [!NOTE]
> Only the paired-end reads of the project were processed with the Part 1 pipeline.

---

## 2. Set-up and the workspace image

### 2.1 What the metadata looks like

One row per sample, three columns:

| Sample | Participant | Timepoint |
|---|---|---|
| S01 | P01 | Baseline |
| S02 | P01 | Day4 |
| S03 | P01 | Day8 |
| S04 | P01 | Day42 |
| S05 | P01 | Day180 |

`Sample` must match the `SampleID` from the Part 1 manifest. Timepoints are ordered, with **Baseline first**: every later timepoint is compared with the first one.

### 2.2 R packages

Tested with R 4.3 or newer.

```r
install.packages(c(
  "magrittr", "dplyr", "tidyr", "tibble", "readr", "stringr", "purrr",
  "ggplot2", "scales", "jsonlite", "vegan", "permute", "lme4",
  "lmerTest", "broom.mixed", "BiocManager"
))
BiocManager::install("biobakery/maaslin3")   # route given in the MaAsLin3 documentation
```

### 2.3 Working directory and workspace image

> [!IMPORTANT]
> **Run this.** Change `working_directory` to a folder on your computer, and put the workspace image we gave you in its `rdata/` sub-folder.

The script keeps everything in one working directory with two sub-folders:

```
your_working_directory/
├── figures/    every plot the script saves
└── rdata/      the workspace image, MaAsLin3 output and result tables
```

```r
# ---- Packages --------------------------------------------------------------------
required_packages <- c(
  "magrittr", "dplyr", "tidyr", "tibble", "readr", "stringr", "purrr",
  "ggplot2", "scales", "jsonlite", "vegan", "permute", "lme4",
  "lmerTest", "broom.mixed", "maaslin3"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop("Install these packages first: ", paste(missing_packages, collapse = ", "))
}
library(magrittr)   # provides %>%; every other function is written as package::function

# ---- Working directory with figures/ and rdata/ ------------------------------------
working_directory <- "~/dawg_workshop_part2"      # <-- CHANGE THIS
dir.create(working_directory, showWarnings = FALSE, recursive = TRUE)
setwd(working_directory)
dir.create("figures", showWarnings = FALSE)
dir.create("rdata",   showWarnings = FALSE)

# ---- Workspace image: load it if it exists, otherwise it is built and saved later ---
workspace_image_path <- file.path("rdata", "part2_workspace.RData")
loaded_from_image    <- FALSE
if (file.exists(workspace_image_path)) {
  load(workspace_image_path)      # brings back every object the presenter had in memory
  loaded_from_image <- TRUE       # set AFTER load(), which would otherwise overwrite it
  message("Loaded workspace image: ", workspace_image_path)
} else {
  message("No workspace image found: the tables will be read from the Part 1 output.")
}

# ---- Analysis settings -----------------------------------------------------------------
set.seed(1)
minimum_post_host_reads  <- 100000   # samples with fewer reads (R1 + R2) are dropped
min_prevalence           <- 0.10     # a feature must be detected in >= 10% of samples
permutations             <- 999      # number of shuffles in the permutation tests
fdr_threshold            <- 0.10     # false discovery rate cut-off for differential abundance
minimum_percent_coverage <- 80       # resistance genes: >= 80% of the gene must be covered
minimum_mapped_reads     <- 5        # resistance genes: >= 5 reads must map
```

> [!TIP]
> If you see `Loaded workspace image`, you are set. The next section is skipped for you.

### 2.4 Reading the Part 1 tables

> [!NOTE]
> **Read only: skipped when the workspace image is loaded.** The workshop presenters ran this block once to build the image. It reads the files written by the Part 1 pipeline and puts them into tidy tables. It is not run in the hands-on session, and attendees do not have the raw data or the `results/` folder.

```r
if (!loaded_from_image) {

  # ---- Where the Part 1 output and the metadata live (presenter's computer) --------
  results_directory <- "/path/to/workshop/results"   # output folder of Part 1
  metadata_path     <- "/path/to/metadata.tsv"       # see section 2.1
  # Timepoint order, baseline FIRST. Edit to match the labels in your metadata.
  timepoint_levels  <- c("Baseline", "Day4", "Day8", "Day42", "Day180")

  # ---- Metadata ----------------------------------------------------------------------
  metadata_raw <- readr::read_tsv(metadata_path, show_col_types = FALSE) %>%
    dplyr::select(Sample, Participant, Timepoint) %>%
    dplyr::mutate(
      Participant = as.character(Participant),
      Timepoint   = factor(Timepoint, levels = timepoint_levels)
    )
  stopifnot(
    "Sample IDs must be unique"            = !anyDuplicated(metadata_raw$Sample),
    "Some Timepoint values are not in timepoint_levels" = !anyNA(metadata_raw$Timepoint)
  )

  # ---- Read QC: fastp reports plus read counts after host removal ---------------------
  fastp_qc <- list.files(
    file.path(results_directory, "fastp"),
    pattern = "\\.fastp\\.json$", full.names = TRUE
  ) %>%
    purrr::map(function(fastp_file) {
      fastp_json <- jsonlite::fromJSON(fastp_file)
      return(tibble::tibble(
        Sample           = stringr::str_remove(basename(fastp_file), "\\.fastp\\.json$"),
        reads_raw        = fastp_json$summary$before_filtering$total_reads,
        reads_post_fastp = fastp_json$summary$after_filtering$total_reads
      ))
    }) %>%
    purrr::list_rbind()

  host_removal_counts <- list.files(
    file.path(results_directory, "qc"),
    pattern = "\\.read_counts\\.tsv$", full.names = TRUE
  ) %>%
    purrr::map(readr::read_tsv, show_col_types = FALSE) %>%
    purrr::list_rbind() %>%
    dplyr::rename(Sample = SampleID)

  read_qc <- fastp_qc %>%
    dplyr::inner_join(host_removal_counts, by = "Sample") %>%
    dplyr::filter(Sample %in% metadata_raw$Sample)

  # ---- MetaPhlAn profiles (all taxonomic ranks, one row per clade per sample) ---------
  # Profiles were made with -t rel_ab_w_read_stats, so column 5 is the estimated number
  # of reads from each clade. `coverage` is "-" on the UNCLASSIFIED row, so column
  # types are declared rather than guessed. Lines starting with # are headers.
  metaphlan_long <- list.files(
    file.path(results_directory, "metaphlan"),
    pattern = "\\.metaphlan\\.txt$", full.names = TRUE
  ) %>%
    purrr::map(function(profile_file) {
      profile_rows <- readr::read_tsv(
        profile_file,
        comment   = "#",
        col_names = c("clade_name", "clade_taxid", "relative_abundance",
                      "coverage", "estimated_reads"),
        col_types = readr::cols(
          clade_name         = readr::col_character(),
          clade_taxid        = readr::col_character(),
          relative_abundance = readr::col_double(),
          coverage           = readr::col_character(),
          estimated_reads    = readr::col_double()
        )
      ) %>%
        dplyr::mutate(Sample = stringr::str_remove(basename(profile_file),
                                                   "\\.metaphlan\\.txt$"))
      return(profile_rows)
    }) %>%
    purrr::list_rbind() %>%
    dplyr::filter(Sample %in% metadata_raw$Sample)

  # ---- RGI bwt gene-level tables ----------------------------------------------------------
  # Read as text first: a sample with no hits writes a header-only file, and an empty
  # column would otherwise be guessed as logical and clash when the tables are joined.
  rgi_genes <- list.files(
    file.path(results_directory, "rgi"),
    pattern = "\\.gene_mapping_data\\.txt$", full.names = TRUE
  ) %>%
    purrr::map(function(rgi_file) {
      rgi_rows <- readr::read_tsv(
        rgi_file, col_types = readr::cols(.default = readr::col_character())
      ) %>%
        dplyr::mutate(Sample = stringr::str_remove(basename(rgi_file),
                                                   "\\.gene_mapping_data\\.txt$"))
      return(rgi_rows)
    }) %>%
    purrr::list_rbind() %>%
    dplyr::filter(Sample %in% metadata_raw$Sample) %>%
    dplyr::transmute(
      Sample,
      aro_term            = `ARO Term`,
      drug_class          = `Drug Class`,
      mapped_reads        = as.numeric(`All Mapped Reads`),
      percent_coverage    = as.numeric(`Average Percent Coverage`),
      reference_length_bp = as.numeric(`Reference Length`)
    )

  # ---- MicrobeCensus: genome equivalents per sample -------------------------------------
  genome_equivalents <- list.files(
    file.path(results_directory, "microbecensus"),
    pattern = "\\.microbecensus\\.txt$", full.names = TRUE
  ) %>%
    purrr::map(function(census_file) {
      census_text <- paste(readLines(census_file, warn = FALSE), collapse = "\n")
      return(tibble::tibble(
        Sample = stringr::str_remove(basename(census_file), "\\.microbecensus\\.txt$"),
        genome_equivalents = as.numeric(stringr::str_match(
          census_text, "genome_equivalents[:\\s]+([0-9.eE+-]+)")[, 2])
      ))
    }) %>%
    purrr::list_rbind() %>%
    dplyr::filter(Sample %in% metadata_raw$Sample)
}
```

Objects created by that block, and available to you after loading the image: `metadata_raw`, `read_qc`, `metaphlan_long`, `rgi_genes`, `genome_equivalents`.

---

## 3. Statistical design at a glance

**Terms you will meet**

- **Mixed model:** a regression that knows some samples come from the same person. It gives each person their own starting level (a *random intercept*), so repeated samples are not treated as independent.
- **Bray-Curtis dissimilarity:** compares two samples by how much of each species' abundance they share. 0 means identical, 1 means nothing shared.
- **PCoA:** a map in which samples that are similar sit close together.
- **PERMANOVA:** asks whether the groups on that map differ more than expected by chance, by shuffling the group labels many times.
- **Q-value (FDR):** a p-value corrected for testing many features at once. We call a result significant when q is below the cut-off.

| Question | What we do | Main R functions |
|---|---|---|
| Is the data quality adequate? | Count reads left after each cleaning step and drop samples with too few | `dplyr` |
| How many species or resistance genes does a sample have (richness)? | Count what is detected. Compare timepoints with a mixed model, with sequencing depth as an extra term because deeper samples detect more | `lme4::glmer.nb`, `lmerTest::lmer` |
| How diverse and even is a bacterial sample (Shannon index)? | Mixed model comparing timepoints | `vegan::diversity`, `lmerTest::lmer` |
| How different are samples overall? | Bray-Curtis dissimilarity, drawn as a PCoA map | `vegan::vegdist`, `stats::cmdscale` |
| Does composition differ between timepoints? | PERMANOVA, shuffling timepoint labels only *within* each person | `vegan::adonis2`, `permute::how` |
| Are samples equally spread out at each timepoint? | PERMDISP (a companion check to PERMANOVA) | `vegan::betadisper` |
| How far does each person drift from their own baseline? | Bray-Curtis distance to their own baseline sample, then a mixed model | `lmerTest::lmer` |
| Which species change? | MaAsLin3: one model for how much of a species is present, one for whether it is detected at all | `maaslin3::maaslin3` |
| Which resistance genes (ARO terms) and drug classes change? | Scale reads to a per-genome value (RPKG), then MaAsLin3 | `maaslin3::maaslin3` |

> [!NOTE]
> Some choices have no field-wide standard: the sequencing-depth term in the richness models, and the coverage and read-count cut-offs used to filter resistance genes. They are analyst decisions, shown as settings at the top, and should be reported in your methods.

---

## 4. Read QC

> [!IMPORTANT]
> **Run this.** Samples with too few reads are removed here, and the log of the read count becomes a **depth** term used in later models.

fastp reports reads before and after trimming, and Part 1 also counted the reads left after **host removal** (human reads removed). Low depth makes detection-based measures unreliable.

```r
# ---- Flag samples that pass the depth cut-off -------------------------------------------
read_qc_checked <- read_qc %>%
  dplyr::mutate(
    host_percent        = 100 * (reads_post_fastp - reads_post_host) / reads_post_fastp,
    passes_depth_filter = reads_post_host >= minimum_post_host_reads
  )
summary(read_qc_checked[, c("reads_raw", "reads_post_host", "host_percent")])

# ---- Plot: reads kept at each stage (one line per sample) -------------------------------
read_stage_plot <- read_qc_checked %>%
  tidyr::pivot_longer(
    c(reads_raw, reads_post_fastp, reads_post_host),
    names_to = "stage", values_to = "reads"
  ) %>%
  dplyr::mutate(stage = factor(
    stage,
    levels = c("reads_raw", "reads_post_fastp", "reads_post_host"),
    labels = c("Raw", "After fastp", "After host removal")
  )) %>%
  ggplot2::ggplot(ggplot2::aes(x = stage, y = reads, group = Sample)) +
  ggplot2::geom_line(colour = "grey60", alpha = 0.6) +
  ggplot2::geom_point(ggplot2::aes(colour = passes_depth_filter), size = 1.5) +
  ggplot2::scale_y_log10(labels = scales::label_comma()) +
  ggplot2::scale_colour_manual(values = c("TRUE" = "black", "FALSE" = "red")) +
  ggplot2::labs(x = NULL, y = "Reads (R1 + R2)", colour = "Passes depth filter") +
  ggplot2::theme_bw() +
  ggplot2::theme(legend.position = "bottom")
print(read_stage_plot)
ggplot2::ggsave(file.path("figures", "read_qc_stages.pdf"),
                read_stage_plot, width = 120, height = 90, units = "mm")

# ---- Keep passing samples and add the depth terms ----------------------------------------
# `metadata` is the table used for everything below.
metadata <- metadata_raw %>%
  dplyr::inner_join(
    dplyr::select(read_qc_checked, Sample, reads_post_host, passes_depth_filter),
    by = "Sample"
  ) %>%
  dplyr::filter(passes_depth_filter) %>%
  dplyr::mutate(
    depth_log10        = log10(reads_post_host),
    depth_log10_scaled = as.numeric(scale(depth_log10))
  )
message(nrow(metadata), " samples from ", dplyr::n_distinct(metadata$Participant),
        " participants kept; ", nrow(metadata_raw) - nrow(metadata),
        " removed for low depth or missing QC.")
baseline_level <- levels(metadata$Timepoint)[1]
```

> [!TIP]
> Look for a steep drop at the host-removal stage (a human-rich or mislabelled sample) and for samples sitting close to the depth cut-off.

---

## 5. Bacteria

### 5.1 The species table

> [!IMPORTANT]
> **Run this.** It turns the MetaPhlAn profiles into a samples-by-species table of read counts.

MetaPhlAn writes a hierarchy: a read assigned to a species is *also* counted at its genus, family and so on. Mixing ranks in one table would count the same read many times, so we keep only **species-level rows** (rows containing `|s__` but no `|t__`).

```r
species_long <- metaphlan_long %>%
  dplyr::filter(
    stringr::str_detect(clade_name, stringr::fixed("|s__")),
    !stringr::str_detect(clade_name, stringr::fixed("|t__")),
    stringr::str_detect(clade_name, "^k__(Bacteria|Archaea)"),
    Sample %in% metadata$Sample
  ) %>%
  dplyr::transmute(
    Sample,
    species         = stringr::str_match(clade_name, "\\|s__([^|]+)")[, 2],
    estimated_reads = round(estimated_reads)
  )

# Samples x species matrix of read counts, rows in the same order as `metadata`
species_counts <- species_long %>%
  tidyr::pivot_wider(names_from = species, values_from = estimated_reads,
                     values_fn = sum, values_fill = 0) %>%
  tibble::column_to_rownames("Sample") %>%
  as.matrix()
stopifnot("Some samples have no MetaPhlAn species" =
            all(metadata$Sample %in% rownames(species_counts)))
species_counts <- species_counts[metadata$Sample, , drop = FALSE]

# A second table keeping species found in >= 10% of samples, used for beta diversity.
# Alpha diversity uses the UNfiltered table, because filtering would change richness.
species_counts_filtered <- species_counts[
  , colMeans(species_counts > 0) >= min_prevalence, drop = FALSE
]
message(ncol(species_counts), " species detected; ", ncol(species_counts_filtered),
        " found in at least ", 100 * min_prevalence, "% of samples.")
```

### 5.2 Alpha diversity

> [!IMPORTANT]
> **Run this.**

- **Richness:** number of species detected. It depends on sequencing depth, so the model includes depth as an extra term.
- **Shannon index:** combines how many species there are and how evenly abundant they are.

Each person is sampled repeatedly, so every model includes a **random intercept per participant**. Richness is a count, so it uses a negative-binomial model (a count model that allows extra spread). The model output compares **each timepoint with baseline**.

```r
stopifnot(identical(rownames(species_counts), metadata$Sample))

bacteria_alpha <- metadata %>%
  dplyr::mutate(
    Richness = rowSums(species_counts > 0),
    Shannon  = vegan::diversity(species_counts, index = "shannon")
  )

# Richness depends on depth: how strong is the relationship?
stats::cor(bacteria_alpha$Richness, bacteria_alpha$depth_log10, method = "spearman")

bacteria_alpha_plot <- bacteria_alpha %>%
  tidyr::pivot_longer(c(Richness, Shannon), names_to = "metric", values_to = "value") %>%
  ggplot2::ggplot(ggplot2::aes(x = Timepoint, y = value)) +
  ggplot2::geom_line(ggplot2::aes(group = Participant), colour = "grey60", alpha = 0.5) +
  ggplot2::geom_boxplot(ggplot2::aes(fill = Timepoint), outlier.shape = NA, alpha = 0.7) +
  ggplot2::facet_wrap(~ metric, scales = "free_y") +
  ggplot2::scale_fill_viridis_d() +
  ggplot2::labs(x = NULL, y = NULL) +
  ggplot2::theme_bw() +
  ggplot2::theme(legend.position = "none")
print(bacteria_alpha_plot)
ggplot2::ggsave(file.path("figures", "bacteria_alpha.pdf"),
                bacteria_alpha_plot, width = 150, height = 70, units = "mm")

# ---- Models: each timepoint versus baseline ---------------------------------------------
shannon_model <- lmerTest::lmer(Shannon ~ Timepoint + (1 | Participant),
                                data = bacteria_alpha)
richness_model <- lme4::glmer.nb(
  Richness ~ Timepoint + depth_log10_scaled + (1 | Participant),
  data    = bacteria_alpha,
  control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
)

bacteria_alpha_results <- dplyr::bind_rows(
  dplyr::mutate(broom.mixed::tidy(shannon_model, effects = "fixed"),
                response = "Shannon"),
  dplyr::mutate(broom.mixed::tidy(richness_model, effects = "fixed"),
                response = "Richness")
) %>%
  dplyr::filter(term != "(Intercept)") %>%
  dplyr::select(response, term, estimate, std.error, statistic, p.value)
readr::write_tsv(bacteria_alpha_results, file.path("rdata", "bacteria_alpha_models.tsv"))
bacteria_alpha_results
```

> [!TIP]
> Each `Timepoint` row compares that timepoint with baseline. A negative Shannon estimate means lower diversity than at baseline. For richness the estimate is on the log scale: `exp(estimate)` is the fold-change in the number of species.

> [!WARNING]
> A message like `boundary (singular) fit` means the person-to-person variation was estimated as zero, so the model behaves like an ordinary regression. It is a warning to be aware of, not an error.

### 5.3 Beta diversity

> [!IMPORTANT]
> **Run this.** Beta diversity measures how different samples are *from each other*. We compare each person only with **themselves over time**.

```r
# ---- Bray-Curtis dissimilarity on relative abundance --------------------------------------
bacteria_relative <- vegan::decostand(
  species_counts_filtered[rowSums(species_counts_filtered) > 0, , drop = FALSE],
  method = "total"
)
bray_bacteria <- vegan::vegdist(bacteria_relative, method = "bray")

# Metadata in exactly the same order as the distances (vegan requires this)
metadata_bacteria <- metadata[match(labels(bray_bacteria), metadata$Sample), ]

# ---- PCoA map: lines join each person's samples in time order ------------------------------
pcoa_bacteria <- stats::cmdscale(bray_bacteria, k = 2, eig = TRUE)
pcoa_percent  <- 100 * pcoa_bacteria$eig[1:2] / sum(pcoa_bacteria$eig[pcoa_bacteria$eig > 0])

bacteria_pcoa_plot <- tibble::tibble(
  Sample = rownames(pcoa_bacteria$points),
  PC1    = pcoa_bacteria$points[, 1],
  PC2    = pcoa_bacteria$points[, 2]
) %>%
  dplyr::inner_join(metadata, by = "Sample") %>%
  dplyr::arrange(Participant, Timepoint) %>%
  ggplot2::ggplot(ggplot2::aes(x = PC1, y = PC2)) +
  ggplot2::geom_path(
    ggplot2::aes(group = Participant), colour = "grey60", alpha = 0.6,
    arrow = grid::arrow(length = grid::unit(1.5, "mm"), type = "closed")
  ) +
  ggplot2::geom_point(ggplot2::aes(fill = Timepoint),
                      shape = 21, size = 2.5, colour = "black", stroke = 0.3) +
  ggplot2::scale_fill_viridis_d() +
  ggplot2::coord_equal() +
  ggplot2::labs(x = paste0("PC1 (", round(pcoa_percent[1], 1), "%)"),
                y = paste0("PC2 (", round(pcoa_percent[2], 1), "%)")) +
  ggplot2::theme_bw() +
  ggplot2::theme(legend.position = "bottom")
print(bacteria_pcoa_plot)
ggplot2::ggsave(file.path("figures", "bacteria_pcoa.pdf"),
                bacteria_pcoa_plot, width = 120, height = 110, units = "mm")
```

**Do timepoints differ in overall composition? (PERMANOVA)** Samples from one person are not interchangeable, so timepoint labels are shuffled **only within each person**.

```r
within_participant_design <- permute::how(
  nperm = permutations, blocks = metadata_bacteria$Participant
)
bacteria_permanova <- vegan::adonis2(
  bray_bacteria ~ Timepoint,
  data         = metadata_bacteria,
  permutations = within_participant_design
) %>%
  tibble::as_tibble(rownames = "term")
readr::write_tsv(bacteria_permanova, file.path("rdata", "bacteria_permanova.tsv"))
bacteria_permanova
```

> [!TIP]
> Report the **R2** value as well as the p-value: R2 is the fraction of the community differences that timepoint explains.

**Are samples equally spread out at each timepoint? (PERMDISP)** PERMANOVA can be significant simply because one timepoint is more variable than another, so we check the spread as well.

```r
bacteria_dispersion <- vegan::betadisper(bray_bacteria,
                                         group = metadata_bacteria$Timepoint,
                                         type  = "centroid")
vegan::permutest(bacteria_dispersion, permutations = within_participant_design)
```

**How far does each person drift from their own baseline?** For every later sample, take the Bray-Curtis distance to the same person's baseline sample.

```r
bray_matrix_bacteria <- as.matrix(bray_bacteria)

baseline_samples <- metadata_bacteria %>%
  dplyr::filter(Timepoint == baseline_level) %>%
  dplyr::select(Participant, baseline_sample = Sample)
stopifnot("Need exactly one baseline sample per participant" =
            !anyDuplicated(baseline_samples$Participant))

bacteria_within <- metadata_bacteria %>%
  dplyr::filter(Timepoint != baseline_level) %>%
  dplyr::inner_join(baseline_samples, by = "Participant") %>%
  dplyr::mutate(
    bray_to_baseline = bray_matrix_bacteria[cbind(Sample, baseline_sample)],
    Timepoint        = droplevels(Timepoint)
  )

bacteria_within_plot <- ggplot2::ggplot(
  bacteria_within, ggplot2::aes(x = Timepoint, y = bray_to_baseline)
) +
  ggplot2::geom_line(ggplot2::aes(group = Participant), colour = "grey60", alpha = 0.5) +
  ggplot2::geom_boxplot(ggplot2::aes(fill = Timepoint), outlier.shape = NA, alpha = 0.7) +
  ggplot2::scale_fill_viridis_d() +
  ggplot2::labs(x = NULL, y = paste0("Bray-Curtis distance to ", baseline_level)) +
  ggplot2::theme_bw() +
  ggplot2::theme(legend.position = "none")
print(bacteria_within_plot)
ggplot2::ggsave(file.path("figures", "bacteria_distance_to_baseline.pdf"),
                bacteria_within_plot, width = 90, height = 80, units = "mm")

# The first post-baseline timepoint is the reference here, because baseline rows are absent.
bacteria_within_results <- lmerTest::lmer(
  bray_to_baseline ~ Timepoint + (1 | Participant), data = bacteria_within
) %>%
  broom.mixed::tidy(effects = "fixed")
readr::write_tsv(bacteria_within_results,
                 file.path("rdata", "bacteria_distance_to_baseline_models.tsv"))
bacteria_within_results
```

> [!TIP]
> If the antibiotics disturbed the community and it then recovered, distance to baseline should peak early and fall at later timepoints.

### 5.4 Differential abundance with MaAsLin3

> [!IMPORTANT]
> **Run this.** This step answers: *which individual species change over time?*

**MaAsLin3** fits two models for every species:

- an **abundance** model: when the species is present, is there more or less of it?
- a **prevalence** model: is the species more or less likely to be detected at all?

Both include the participant random intercept. Two settings matter:

- Abundance is *relative* (each sample sums to one), so if one species rises, others appear to fall. MaAsLin3 guards against this by comparing each species with the *median* change across species (`median_comparison_abundance = TRUE`) instead of with zero.
- Sequencing depth is included as a term, because deeper samples detect more species and would otherwise distort the prevalence model.

```r
# Rows = samples, with the Sample ID as row names (they must match the species table)
maaslin_metadata <- metadata %>%
  dplyr::mutate(Participant = as.factor(Participant)) %>%
  dplyr::select(Sample, Participant, Timepoint, depth_log10) %>%
  tibble::column_to_rownames("Sample")

differential_formula <- "~ Timepoint + depth_log10 + (1|Participant)"

bacteria_da_directory <- file.path("rdata", "maaslin3_bacteria_species")
set.seed(1)
maaslin3::maaslin3(
  input_data                   = species_counts,   # raw counts: MaAsLin3 scales them itself
  input_metadata               = maaslin_metadata,
  output                       = bacteria_da_directory,
  formula                      = differential_formula,
  normalization                = "TSS",
  transform                    = "LOG",
  standardize                  = TRUE,
  min_prevalence               = min_prevalence,
  max_significance             = fdr_threshold,
  median_comparison_abundance  = TRUE,
  median_comparison_prevalence = FALSE,
  cores                        = 1
)

bacteria_da_significant <- readr::read_tsv(
  file.path(bacteria_da_directory, "all_results.tsv"), show_col_types = FALSE
) %>%
  dplyr::filter(metadata != "depth_log10", !is.na(qval_individual),
                qval_individual < fdr_threshold) %>%
  dplyr::arrange(qval_individual) %>%
  dplyr::select(dplyr::any_of(c("feature", "value", "model", "coef", "stderr",
                                "N", "N_not_zero", "qval_individual")))
readr::write_tsv(bacteria_da_significant,
                 file.path("rdata", "bacteria_species_maaslin3_significant.tsv"))
bacteria_da_significant
```

**Reading the table**

- `value` is the timepoint compared with baseline. `model` is `abundance` or `prevalence`.
- In an abundance model, a `coef` of *c* means a 2^*c*-fold change in relative abundance. In a prevalence model it is a change in the log-odds of the species being detected.
- `qval_individual` is the corrected p-value for that one association.
- The output folder also holds `summary_plot.pdf` and one plot per association. Open them before trusting a hit.

> [!TIP]
> Antibiotics should reduce many commensals at Day4 and Day8 and let others rise. Compare your top hits with the species named in the study (section 1).

---

## 6. Antimicrobial resistance genes

### 6.1 From read counts to RPKG

> [!IMPORTANT]
> **Run this.**

RGI bwt reports how many reads mapped to each resistance gene (identified by its **ARO term**, a CARD database label). Raw counts cannot be compared directly, because (a) deeper samples have more reads on everything, (b) longer genes catch more reads, and (c) communities differ in how much DNA each genome contributes. **RPKG** fixes all three: *reads per kilobase of gene per genome equivalent*. It divides by gene length in kb, then by the number of genome equivalents (MicrobeCensus) in the sample.

Two filters remove spurious low-level hits: at least 80% of the gene covered, and at least 5 reads mapped. These are workshop defaults, **not** an official standard. RGI leaves this choice to the analyst, so report what you use.

Genes are summarised two ways: by individual **ARO term** and by **drug class**. A gene can belong to several drug classes and is counted in each, so drug-class totals overlap and must not be added together.

```r
rgi_rpkg <- rgi_genes %>%
  dplyr::filter(Sample %in% metadata$Sample,
                percent_coverage >= minimum_percent_coverage,
                mapped_reads     >= minimum_mapped_reads) %>%
  dplyr::inner_join(genome_equivalents, by = "Sample") %>%
  dplyr::mutate(
    rpkg = (mapped_reads / (reference_length_bp / 1000)) / genome_equivalents
  )
stopifnot("No resistance genes pass the filters; relax the thresholds" =
            nrow(rgi_rpkg) > 0)

# Samples x feature matrices of RPKG. Samples with no passing gene get zeros (not NA).
# Drug classes are separated by "; " when a gene belongs to several.
amr_matrices <- list()
for (feature_column in c("aro_term", "drug_class")) {
  feature_wide <- rgi_rpkg %>%
    tidyr::separate_longer_delim(dplyr::all_of(feature_column), delim = "; ") %>%
    dplyr::group_by(Sample, feature = .data[[feature_column]]) %>%
    dplyr::summarise(rpkg = sum(rpkg), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = feature, values_from = rpkg, values_fill = 0)

  amr_matrices[[feature_column]] <- metadata %>%
    dplyr::select(Sample) %>%
    dplyr::left_join(feature_wide, by = "Sample") %>%
    dplyr::mutate(dplyr::across(-Sample, function(rpkg_values) {
      return(tidyr::replace_na(rpkg_values, 0))
    })) %>%
    tibble::column_to_rownames("Sample") %>%
    as.matrix()
}
aro_matrix <- amr_matrices$aro_term
message(ncol(aro_matrix), " resistance genes (ARO terms) and ",
        ncol(amr_matrices$drug_class), " drug classes pass the filters.")
```

### 6.2 Richness of resistance genes

> [!IMPORTANT]
> **Run this.** Resistance-gene **richness** is the number of different ARO terms detected in a sample (after the filters above).

It is a detection count, so it depends on sequencing depth. The model includes depth, exactly like bacterial richness.

```r
stopifnot(identical(rownames(aro_matrix), metadata$Sample))

amr_richness <- metadata %>%
  dplyr::mutate(Richness = rowSums(aro_matrix > 0))

amr_richness_plot <- ggplot2::ggplot(
  amr_richness, ggplot2::aes(x = Timepoint, y = Richness)
) +
  ggplot2::geom_line(ggplot2::aes(group = Participant), colour = "grey60", alpha = 0.5) +
  ggplot2::geom_boxplot(ggplot2::aes(fill = Timepoint), outlier.shape = NA, alpha = 0.7) +
  ggplot2::scale_fill_viridis_d() +
  ggplot2::labs(x = NULL, y = "Resistance genes detected (ARO terms)") +
  ggplot2::theme_bw() +
  ggplot2::theme(legend.position = "none")
print(amr_richness_plot)
ggplot2::ggsave(file.path("figures", "amr_richness.pdf"),
                amr_richness_plot, width = 90, height = 80, units = "mm")

amr_richness_model <- lme4::glmer.nb(
  Richness ~ Timepoint + depth_log10_scaled + (1 | Participant),
  data    = amr_richness,
  control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
)
amr_richness_results <- broom.mixed::tidy(amr_richness_model, effects = "fixed") %>%
  dplyr::filter(term != "(Intercept)") %>%
  dplyr::select(term, estimate, std.error, statistic, p.value)
readr::write_tsv(amr_richness_results, file.path("rdata", "amr_richness_model.tsv"))
amr_richness_results
```

> [!TIP]
> As before, each `Timepoint` row compares with baseline and `exp(estimate)` is the fold-change in the number of resistance genes detected.

### 6.3 Differential abundance of resistance genes

> [!IMPORTANT]
> **Run this.** We test which **individual resistance genes (ARO terms)** and which **drug classes** change over time.

The compositionality worry from section 5.4 is different here. RPKG is scaled **per genome equivalent**, not as a share of all resistance reads, so a rise in one gene does not mechanically push the others down. For data already on an absolute-type scale like this, the MaAsLin3 manual's setting is `normalization = "NONE"` with `median_comparison_abundance = FALSE`, so each gene is tested against zero change. Coefficients read as changes in gene copies per genome (on the log scale).

> [!WARNING]
> Choosing "NONE" here is an analyst judgement. Report it in your methods. As a sensitivity check you can repeat the analysis with `normalization = "TSS"` and `median_comparison_abundance = TRUE`, as for bacteria.

> [!NOTE]
> Individual ARO terms are sparse, and reads from near-identical alleles can map to several references, so a single ARO-term hit is less certain than a drug-class result. Treat ARO-term results as hypotheses to follow up.

```r
# Feature names contain spaces and brackets. They are converted to R-safe names for
# the model, then mapped back to the original names in the results.
for (feature_level in c("aro_term", "drug_class")) {
  feature_matrix <- amr_matrices[[feature_level]]
  feature_lookup <- tibble::tibble(
    feature_original = colnames(feature_matrix),
    feature          = make.names(colnames(feature_matrix), unique = TRUE)
  )
  colnames(feature_matrix) <- feature_lookup$feature

  feature_directory <- file.path("rdata", paste0("maaslin3_amr_", feature_level))
  set.seed(1)
  maaslin3::maaslin3(
    input_data                   = feature_matrix,
    input_metadata               = maaslin_metadata,
    output                       = feature_directory,
    formula                      = differential_formula,
    normalization                = "NONE",
    transform                    = "LOG",
    standardize                  = TRUE,
    min_prevalence               = min_prevalence,
    max_significance             = fdr_threshold,
    median_comparison_abundance  = FALSE,
    median_comparison_prevalence = FALSE,
    cores                        = 1
  )

  feature_significant <- readr::read_tsv(
    file.path(feature_directory, "all_results.tsv"), show_col_types = FALSE
  ) %>%
    dplyr::filter(metadata != "depth_log10", !is.na(qval_individual),
                  qval_individual < fdr_threshold) %>%
    dplyr::left_join(feature_lookup, by = "feature") %>%
    dplyr::arrange(qval_individual) %>%
    dplyr::select(feature_original,
                  dplyr::any_of(c("value", "model", "coef", "stderr",
                                  "N", "N_not_zero", "qval_individual")))
  readr::write_tsv(
    feature_significant,
    file.path("rdata", paste0("amr_", feature_level, "_maaslin3_significant.tsv"))
  )
  print(feature_significant)
}
```

> [!TIP]
> Compare the drug-class results with the study: antibiotic exposure included a beta-lactam (meropenem), an aminoglycoside (gentamicin) and a glycopeptide (vancomycin). The authors did not report a simple enrichment of resistance genes for all three classes, so an unexpected or absent signal is not necessarily a mistake.

---

## 7. Save the workspace image

> [!NOTE]
> **Presenters only.** This runs only when the script built everything from the Part 1 tables (no image was loaded). It saves **every object in memory**, so attendees can load it and get all intermediate tables without the raw data. Attendees who loaded the image skip it automatically.

```r
if (!loaded_from_image) {
  save.image(file = workspace_image_path)
  message("Saved workspace image: ", workspace_image_path)
}
```

The image holds data only. Attendees still need the R packages from section 2.2, and they should use a similar R version to the one that created it.

---

## 8. Pitfalls

- **Depth.** Richness and prevalence follow sequencing depth. Keep the depth term in the models and report it.
- **Repeated samples.** Treating samples from the same person as independent makes results look more significant than they are. Keep the participant random intercept, and shuffle labels only within people for PERMANOVA.
- **Small studies.** With 12 participants, models can be singular and p-values coarse. Report effect sizes, not just p-values.
- **Many tests.** The q-value corrects across all tests in the model, including the depth term.
- **Overlapping drug classes.** A gene is counted once in every class it belongs to, so classes cannot be added up.
- **Thresholds are choices.** The coverage, read-count and prevalence cut-offs change the answer. Try alternatives and see whether your conclusions hold.
- **Database versions.** Report the MetaPhlAn index and the CARD version (both recorded by Part 1).

---

## 9. References

- Palleja A, et al. Recovery of gut microbiota of healthy adults following antibiotic exposure. *Nat Microbiol* 3, 1255-1265 (2018). doi:10.1038/s41564-018-0257-9
- Nickols WA, et al. MaAsLin 3: refining and extending generalized multivariable linear models for meta-omic association discovery. *bioRxiv* (2024). doi:10.1101/2024.12.13.628459
- Anderson MJ. A new method for non-parametric multivariate analysis of variance. *Austral Ecology* 26, 32-46 (2001). (PERMANOVA)
- Anderson MJ. Distance-based tests for homogeneity of multivariate dispersions. *Biometrics* 62, 245-253 (2006). (PERMDISP)
- Oksanen J, et al. *vegan: Community Ecology Package* (R package); Simpson GL. *permute* (R package).
- Bates D, Maechler M, Bolker B, Walker S. Fitting linear mixed-effects models using lme4. *J Stat Softw* 67(1), 1-48 (2015).
- Kuznetsova A, Brockhoff PB, Christensen RHB. lmerTest package: tests in linear mixed effects models. *J Stat Softw* 82(13), 1-26 (2017).
- Nayfach S, Pollard KS. Average genome size estimation improves comparative metagenomics and sheds light on the functional ecology of the human microbiome. *Genome Biol* 16, 51 (2015). (MicrobeCensus, RPKG)
- Blanco-Miguez A, et al. *Nat Biotechnol* 41, 1633-1644 (2023). (MetaPhlAn 4)
- Alcock BP, et al. *Nucleic Acids Res* 51, D690-D699 (2023). (CARD)

**Back to:** [Part 1: Profiling reads](01_run_profiling.md)
