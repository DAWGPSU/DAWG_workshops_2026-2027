# DAWG Workshop #1 — 16S Amplicon Sequencing (DADA2 + ALDEx3)

## Overview

This workshop walks trainees through a complete 16S rRNA amplicon
sequencing analysis — from raw paired-end FASTQ reads to a differential
abundance result — entirely in R, on their own laptop, with no HPC
account and no QIIME 2 install.

Trainees compare **vegetated** and **bare** soil sites from the Atacama
Desert (Neilson et al. 2017) using **DADA2** for denoising,
**phyloseq/vegan** for diversity statistics, and **ALDEx3** for
differential abundance.

This folder is for code version control only — trainee materials are
distributed separately via OneDrive.

**Trainee data & scripts (zip):** [[OneDrive link](https://pennstateoffice365-my.sharepoint.com/:f:/g/personal/ako5205_psu_edu/IgBjXeYu747DSJxtg5BBC7MsAe86MtynooxtQXxxZ1mF6j0?e=opRf4A)]

**R version:** 4.4 or newer, on your own laptop (no HPC account needed).
Please download, unzip, run `setup_packages.R` once.

---

## The Dataset

**Atacama Desert soil microbiome** (Neilson et al. 2017, *mSystems*,
[doi:10.1128/mSystems.00195-16](https://doi.org/10.1128/mSystems.00195-16))

* Also the official QIIME 2 tutorial dataset — well-documented, widely
  used for teaching
* Sequenced with the same V4 (515F/806R) primers and Illumina platform
  the slide deck already covers

---


## Slides

The lecture half of the workshop is a separate slide deck covering amplicon sequencing background, the DADA2
algorithm, compositional-data statistics, and the ALDEx3 method — the
conceptual counterpart to the hands-on scripts and walkthrough doc here.

---


## What You'll Need

- This document for instructions and code. Preferably as a separate tab.
- R (4.4+) and RStudio, installed on your own laptop.
- **`setup_packages.R` already run** — see this folder's `README.md` for
  the one-time setup. 

## Assumptions

This workflow assumes that your sequencing data:

- Is for 16S rRNA, not ITS (fungal).
- Has already been demultiplexed (split into per-sample FASTQ files) —
  done for you, already sitting in `data/` in this folder.
- Is paired-end, with forward and reverse FASTQ files for each sample.

------------------------------------------------------------------------

## Getting Set Up

1. Unzip the workshop folder you downloaded somewhere on your computer.
2. Open R or RStudio with **this folder** (the one containing `data/` and
   `scripts/`) as your working directory. In RStudio: `Session > Set
   Working Directory > Choose Directory`, and pick it. (If a `.Rproj` file
   is included, just double-click it instead — that sets the working
   directory for you.)
3. Check you're in the right place:

``` r
list.files()  # should show "data", "scripts", this .md file, etc.
```

------------------------------------------------------------------------

## Background

**Raw data:** Real soil samples along two elevation/aridity transects in
the Atacama Desert (vegetated and unvegetated sites), 16S V4 region,
515F/806R primers, Illumina MiSeq (2×150 bp). Already demultiplexed and
sitting in `data/` in this folder — see `data/data-info.md`.

**Question we'll actually answer today:** does soil bacterial community
composition and diversity differ between vegetated and bare (unvegetated)
sites? That's the `vegetation` column in the metadata, and it's what we'll
run through ALDEx3 at the end.

------------------------------------------------------------------------

## Workflow Steps

0.  Set Your Working Directory
1.  Load Packages
2.  Import Data
3.  Inspect Read Quality
4.  Denoise Reads with DADA2
5.  Track Reads Through the Pipeline
6.  Build a Phylogenetic Tree & Assign Taxonomy (provided — you load, not build)
7.  Alpha and Beta Diversity
8.  Taxonomic Classification (see Step 6)
9.  Differential Abundance with ALDEx3

------------------------------------------------------------------------

## Step 0: Set Your Working Directory

Everything in this doc assumes R's working directory is this folder (the
one with `data/` and `scripts/` in it) — covered in "Getting Set Up"
above. Quick way to confirm you're set up right:

``` r
list.files("data")  # should include manifest.tsv and sample-metadata.tsv
```

If that errors with "no such file or directory," your working directory
isn't set to the right place yet — fix that before continuing.

------------------------------------------------------------------------

## Step 1: Load Packages

``` r
library(dada2)
library(phyloseq)
library(DECIPHER)
library(phangorn)
library(vegan)
library(tidyverse)
library(ALDEx3)
```

> **Note:** if any `library()` call fails with "there is no package called
> ...", go back and run `setup_packages.R` (see this folder's `README.md`)
> — it installs everything this doc needs, once, ahead of time.

------------------------------------------------------------------------

## Step 2: Import Data

Two supporting files, just like `data-manifest.txt` and `metadata.tsv` in
the QIIME 2 workshop:

- **`manifest.tsv`**: sample name → path to that sample's R1/R2 FASTQ files.
- **`sample-metadata.tsv`**: `vegetation`, `transect-name`, `elevation`,
  and other columns describing each sample.

``` r
manifest <- read_tsv("data/manifest.tsv")
metadata <- read_tsv("data/sample-metadata.tsv")

# Named vectors of file paths -- this is the format DADA2's functions expect
fnFs <- setNames(manifest$forward, manifest$sample)
fnRs <- setNames(manifest$reverse, manifest$sample)

length(fnFs)   # how many samples did we get?
head(metadata) # sanity-check the metadata loaded correctly
```

🔍 **What to check:** does `length(fnFs)` match the number of rows in
`sample-metadata.tsv`? If not, stop here — a mismatch means something went
wrong when the data was prepared, before it ever reached you.

------------------------------------------------------------------------

## Step 3: Inspect Read Quality

This is the R equivalent of `demux summarize` — before deciding where to
trim, look at where quality actually drops off.

``` r
plotQualityProfile(fnFs[1:4])  # forward reads, first 4 samples
plotQualityProfile(fnRs[1:4])  # reverse reads, first 4 samples
```

👀 Each plot shows a heat map of quality score by position, with the
median (green line) and quartiles (orange). Reverse reads almost always
degrade faster than forward reads — that's normal, not a problem with your
data.

🔍 **What to look for:** the position where the green median line starts
dropping below ~Q30 (see slide 12's Phred score refresher). That position
is roughly where you'll want to truncate in Step 4. Do this for a few more
samples, not just the first one — a truncation length chosen from a single
sample can be a bad fit for the rest.

------------------------------------------------------------------------

## Step 4: Denoise Reads with DADA2

This step bundles what were four separate commands in the QIIME 2 doc
(filter, learn errors, denoise, merge) — DADA2's R interface keeps them as
separate function calls, but conceptually it's one job: turn raw reads into
an ASV table.

``` r
filtFs <- file.path("filtered", paste0(names(fnFs), "_F_filt.fastq.gz"))
filtRs <- file.path("filtered", paste0(names(fnRs), "_R_filt.fastq.gz"))
names(filtFs) <- names(fnFs)
names(filtRs) <- names(fnRs)

# --- Filter and trim ---------------------------------------------------------
# truncLen: cut each read to this length, based on what you saw in Step 3.
# Starting point below assumes quality holds up reasonably well to ~140bp on
# both reads (2x150bp MiSeq run) -- adjust to what YOUR quality plots show.
# The V4 amplicon is ~253bp, so forward + reverse need to add up to at least
# ~273bp (253 + ~20bp overlap) after truncation, or mergePairs in the next
# chunk will fail to merge most reads.
out <- filterAndTrim(
  fnFs, filtFs, fnRs, filtRs,
  truncLen = c(140, 140),
  maxN = 0, maxEE = c(2, 2), truncQ = 2, rm.phix = TRUE,
  compress = TRUE, multithread = TRUE
)
head(out)  # reads.in vs reads.out per sample -- lost a lot? loosen maxEE or truncLen.

# --- Learn the error model ----------------------------------------------------
# This is the DADA2 step from the slides: learn what real sequencing error
# looks like in THIS run, from the data itself, before deciding whether a
# rare sequence is a true biological variant or a misread.
errF <- learnErrors(filtFs, multithread = TRUE)
errR <- learnErrors(filtRs, multithread = TRUE)
plotErrors(errF, nominalQ = TRUE)  # sanity check: points should track the black line

# --- Denoise -------------------------------------------------------------------
dadaFs <- dada(filtFs, err = errF, multithread = TRUE)
dadaRs <- dada(filtRs, err = errR, multithread = TRUE)

# --- Merge forward and reverse into full-length ASVs --------------------------
# maxMismatch = 1: with a ~27bp forward/reverse overlap at these truncLen
# values, the default (0 mismatches allowed) throws out a pair if even one
# letter in the overlap was misread. Tolerating 1 real error recovers a
# meaningful chunk of otherwise-good reads.
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, maxMismatch = 1, verbose = TRUE)

# --- Build the ASV table ---------------------------------------------------
seqtab <- makeSequenceTable(mergers)

# Sanity check before we go further: our V4 amplicon should be ~253bp, so
# merged sequences should cluster tightly around that. A lot of sequences
# well outside that range usually means non-specific priming somewhere
# upstream -- worth noticing now, not after taxonomy assignment.
table(nchar(getSequences(seqtab)))

# --- Remove chimeras ---------------------------------------------------------
seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread = TRUE, verbose = TRUE)

dim(seqtab.nochim)                                   # samples x ASVs
sum(seqtab.nochim) / sum(seqtab)                      # fraction of reads kept after chimera removal
```

> ⏱️ `learnErrors()` and `dada()` are the slow steps — a few minutes each,
> depending on how many samples we're running. The terminal will look
> stuck; that's normal.

------------------------------------------------------------------------

## Step 5: Track Reads Through the Pipeline

The R equivalent of `denoising-stats.qzv` — build the same read-tracking
table by hand, since each DADA2 step returns its own counts.

``` r
getN <- function(x) sum(getUniques(x))
track <- tibble(
  sample     = names(fnFs),
  input      = out[, "reads.in"],
  filtered   = out[, "reads.out"],
  denoisedF  = sapply(dadaFs, getN),
  denoisedR  = sapply(dadaRs, getN),
  merged     = sapply(mergers, getN),
  nonchim    = rowSums(seqtab.nochim)
)
track
```

🔍 **What to look for**, same as the slide deck's step 5 warning: how many
reads survive each stage. `merged` is usually where the biggest drop
happens — if it's much lower than `denoisedF`/`denoisedR` *for a sample
that had real depth to begin with*, your forward and reverse reads probably
aren't overlapping enough, which usually means Step 4's `truncLen` was too
aggressive. A sample with very few reads in the first place will always
show a poor merge rate regardless — that's not a truncLen problem, it's
just not enough data, and it's exactly what Step 7's read-count filter
below is for.

------------------------------------------------------------------------

## Step 6: Build a Phylogenetic Tree & Assign Taxonomy (provided)

This step is **provided for you**, not run live — here's why, and what to
do instead.

Building a tree (needed for UniFrac in Step 7) means aligning every ASV
and working out how they relate to each other — for a dataset this size
(thousands of ASVs), that alignment alone takes a few minutes even on a
fast laptop, and a fuller tree-refinement step some tutorials use can run
for a very long time without ever finishing. Multiply that by everyone in
the room running it at once, and it eats the workshop. So instead, this
was already run once (by the instructor, ahead of time) and the result —
a phylogenetic tree, **and** taxonomy assigned against the SILVA
database (the next step, normally) — is sitting ready for you to load:

``` r
ps <- readRDS("data/checkpoints/phyloseq_object.rds")
ps
```

That single object already contains your ASV table, the sample metadata,
a phylogenetic tree, taxonomy assignments, and has had very-low-read
samples (sequencing noise, not real biology) filtered out — everything
Steps 6-8 would otherwise have built. You're picking back up as if you'd
just run all three steps yourself.

🔍 **For your own understanding** (not something to run now): the tree
was built by aligning ASVs with `DECIPHER::AlignSeqs()` and building a
neighbor-joining tree from the resulting distances (`phangorn::NJ()`) —
a fast, standard approach at this scale, skipping the much slower
full maximum-likelihood refinement some tutorials use, which isn't
practical outside small example datasets. Taxonomy was assigned with
`assignTaxonomy()` against the SILVA reference — same idea as Step 8
below. If you want to see or re-run this exact process later on your own
data, `instructor-only/README.md` and `scripts/02_diversity_stats.R` have
the full code.

------------------------------------------------------------------------

## Step 7: Alpha and Beta Diversity

Picking up with the `ps` object you just loaded — already has the ASV
table, metadata, tree, and taxonomy all bundled together.

### Alpha diversity (within-sample)

``` r
alpha <- estimate_richness(ps, measures = c("Observed", "Shannon", "Simpson"))
alpha$sample <- rownames(alpha)
alpha <- left_join(alpha, metadata, by = c("sample" = "sample-id"))

# Boxplot (the standard way to show a distribution by group), with the raw
# points jittered on top so nobody mistakes the box for hiding a small n.
alpha_long <- alpha %>%
  pivot_longer(cols = c(Observed, Shannon, Simpson), names_to = "measure", values_to = "value")

ggplot(alpha_long, aes(x = vegetation, y = value, fill = vegetation)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.15, alpha = 0.5, size = 1) +
  facet_wrap(~ measure, scales = "free_y") +
  labs(title = "Alpha diversity by vegetation", x = NULL, y = NULL) +
  theme_bw() +
  theme(legend.position = "none")

# A plot is not a test (slide 21/22) -- back it up with a rank-based test:
wilcox.test(Shannon ~ vegetation, data = alpha)
```

🔍 Same message as the slides: don't stop at the boxplot. `vegetation` has
two groups, so Wilcoxon is the right test here (Kruskal–Wallis if you're
comparing `transect-name`, which has more than two levels).

### Beta diversity (between-sample)

``` r
# Rarefy first for a fair comparison across samples of very different depth
ps.rare <- rarefy_even_depth(ps, rngseed = 1, replace = FALSE)

bray <- phyloseq::distance(ps.rare, method = "bray")
ord.bray <- ordinate(ps.rare, method = "PCoA", distance = bray)

# stat_ellipse() draws a "typical spread" region per group -- a visual aid
# for whether the two point clouds separate, not a statistical test on its
# own; PERMANOVA below is the actual test.
plot_ordination(ps.rare, ord.bray, color = "vegetation") +
  stat_ellipse(aes(group = vegetation), type = "t", linewidth = 0.8) +
  labs(title = "Bray-Curtis PCoA")

# The plot alone doesn't tell you if groups actually differ -- PERMANOVA does:
adonis2(bray ~ vegetation, data = as(sample_data(ps.rare), "data.frame"))
```

🔍 Report the R² from `adonis2()`, not just the p-value (slide 31, point 5)
— it tells you how much of the community variation `vegetation` actually
explains, which a p-value alone doesn't.

You can swap `method = "bray"` for `"unifrac"` or `"wunifrac"` (now that we
have a tree) to see whether the answer changes with a phylogeny-aware
metric — a good sensitivity check, same spirit as the gamma sensitivity
check in Step 9.

### Genus-level composition: what's actually there

Step 9 (next) asks a strict question — which ASVs are *statistically*
different between groups — and with thousands of ASVs tested at once,
often only a handful (or one) will clear that bar. Here's a simpler,
purely descriptive question instead: of everything living in each sample,
which genera make up most of it, and does that overall mix look different
between bare and vegetated soil? Picture each sample as a pie chart
sliced by genus, then unrolled into a single vertical bar so many samples
can stand side by side.

``` r
# Collapse every ASV into its genus, then convert counts to relative
# abundance (% of that sample's total reads) so samples with different
# sequencing depth are comparable.
ps_genus <- tax_glom(ps, taxrank = "Genus", NArm = FALSE)
ps_genus_rel <- transform_sample_counts(ps_genus, function(x) x / sum(x))

genus_df <- psmelt(ps_genus_rel) %>%
  mutate(Genus = if_else(is.na(Genus), "Unclassified", Genus))

# Keep the plot readable: only the most abundant genera overall get their
# own color (capped at 8), everything else is lumped into "Other."
TOP_N_GENERA <- 8
top_genera <- genus_df %>%
  group_by(Genus) %>%
  summarise(total = sum(Abundance), .groups = "drop") %>%
  arrange(desc(total)) %>%
  slice_head(n = TOP_N_GENERA) %>%
  pull(Genus)

genus_df <- genus_df %>%
  mutate(GenusLabel = if_else(Genus %in% top_genera, Genus, "Other"))

# A colorblind-safe categorical palette, assigned in a fixed order (most
# abundant genus gets slot 1, etc.). "Other" gets a neutral gray on
# purpose -- it's a catch-all, not its own identity.
genus_colors <- c(
  setNames(
    c("#2a78d6", "#eb6834", "#1baf7a", "#eda100",
      "#e87ba4", "#008300", "#4a3aa7", "#e34948")[seq_along(top_genera)],
    top_genera
  ),
  "Other" = "#898781"
)

ggplot(genus_df, aes(x = Sample, y = Abundance, fill = GenusLabel)) +
  geom_col() +
  scale_fill_manual(values = genus_colors) +
  facet_wrap(~ vegetation, scales = "free_x") +
  labs(title = "Genus-level composition by sample",
       x = NULL, y = "Relative abundance", fill = "Genus") +
  theme_bw() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
```

🔍 This is a *description*, not a test — two groups can look visually
similar here and still differ significantly in the alpha/beta diversity
tests above, or in Step 9's ALDEx3 results (or vice versa). It's a
companion to those results, not a replacement for them.

------------------------------------------------------------------------

## Step 8: Taxonomic Classification

Already done — this was included in the `ps` object you loaded back in
Step 6 (taxonomy was assigned there together with the tree, both provided
ahead of time for speed). Nothing to run here.


------------------------------------------------------------------------

## Step 9: Differential Abundance with ALDEx3

Everything up to here described *how similar samples are* (diversity) or
*what's present* (taxonomy). This step answers a different question:
**which specific ASVs differ between vegetated and unvegetated soil** —
and it's built to handle the scale-uncertainty problem from Part 4 of the
slides, rather than pretending the normalization is neutral.

``` r
# ALDEx3 wants taxa x samples -- the transpose of what phyloseq's otu_table
# gives you. Get this backwards and you'll get nonsense, not an error
# (same warning as on the slide).
Y <- t(as(otu_table(ps), "matrix"))
meta_df <- as(sample_data(ps), "data.frame")

fit <- aldex(
  Y,
  ~ vegetation,
  data  = meta_df,
  scale = tss.sm,   # scale model -- see slide 27 for what this assumes
  gamma = 0.5       # admit scale uncertainty instead of pretending we know it
)

res <- summary(fit)
# estimate, std.error, p.val.adj (BH), log2 scale -- same columns as the slide

# The recommended sensitivity check from the slides: compare gamma = 0 (old,
# overconfident behavior) against gamma = 0.5 (realistic). Taxa that are
# "significant" at gamma=0 but vanish at gamma=0.5 were only significant
# because of an unstated assumption about total microbial load.
fit0 <- aldex(Y, ~ vegetation, data = meta_df, scale = tss.sm, gamma = 0)
res0 <- summary(fit0)

sig_gamma0  <- res0 %>% filter(p.val.adj < 0.05)
sig_gamma05 <- res  %>% filter(p.val.adj < 0.05)
cat("Significant at gamma=0:  ", nrow(sig_gamma0),  "\n")
cat("Significant at gamma=0.5:", nrow(sig_gamma05), "\n")

# Attach taxonomy so results read as genera, not ASV hashes. IMPORTANT: join
# on the "entity" column, not on rownames(res) -- ALDEx3's summary() numbers
# each row 1, 2, 3... (just a row counter), while the ASV's actual identity
# (its full DNA sequence, matching tax_table(ps)'s row names) lives in the
# "entity" column. Joining on the row number instead silently matches
# nothing, and every taxonomy column comes back NA -- not because those
# ASVs are unclassified, but because the join key was wrong.
tax_df <- as.data.frame(as(tax_table(ps), "matrix")) %>% rownames_to_column("entity")
res_annotated <- left_join(as.data.frame(res), tax_df, by = "entity") %>%
  rename(ASV = entity)
```

🔍 **What to look for:** which ASVs survive both gamma values (those are
the safer calls), and which only showed up at gamma=0. `res_annotated`
above gives you the genus for each significant ASV — that's the figure
this whole workshop has been building toward.

``` r
# Bar chart: the ASVs ALDEx3 found most different, colored by direction,
# solid = significant. Only the top 20 by adjusted p-value are shown --
# with thousands of ASVs, plotting all of them is unreadable.
plot_df <- res_annotated %>%
  mutate(
    label       = if_else(is.na(Genus), paste0("Unclassified (", substr(ASV, 1, 8), "...)"), Genus),
    label       = paste0(label, " #", row_number()),
    direction   = if_else(estimate > 0, "Enriched: vegetated", "Enriched: bare"),
    significant = p.val.adj < 0.05
  ) %>%
  arrange(p.val.adj) %>%
  slice_head(n = 20)

ggplot(plot_df, aes(x = reorder(label, estimate), y = estimate,
                     fill = direction, alpha = significant)) +
  geom_col() +
  scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.35), guide = "none") +
  coord_flip() +
  labs(
    title    = "Top 20 ASVs by ALDEx3 effect size (gamma = 0.5)",
    subtitle = "Solid bars = significant (BH-adjusted p < 0.05); faded = not significant",
    x = NULL, y = "Effect size (CLR difference, vegetated vs bare)", fill = NULL
  ) +
  theme_bw()
```

------------------------------------------------------------------------

## 🎉 Congratulations!

You've completed the full pipeline — from raw paired-end FASTQ files to a
differential abundance result you can defend. Along the way you inspected
read quality, denoised with DADA2, tracked reads at every step, built a
phylogenetic tree, tested alpha and beta diversity properly (not just
plotted them), assigned taxonomy, and used ALDEx3 to find which taxa
actually differ between vegetated and bare soil.

If everything worked, `scripts/` has the same code as standalone `.R`
files (`01_dada2_pipeline.R`, `02_diversity_stats.R`,
`03_aldex3_differential_abundance.R`) — useful for running the whole thing
start to finish on your own data later.
