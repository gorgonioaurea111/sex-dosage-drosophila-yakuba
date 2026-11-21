# Sex-Chromosome Dosage Analysis in Fruit Fly (*Drosophila yakuba*)

**Course:** Cell and Molecular Biology  
**Activity:** Sex-Chromosome Dosage in a Selected Animal Nuclear Genome Using NCBI, Galaxy, and R  
**Species:** Fruit Fly (*Drosophila yakuba*)  
**Genome Assembly:** Prin_Dyak_Tai18E2_2.1 (GCF_016746365.2)  
**WGS Run:** SRR26246471  

---

## 1. Project Overview

This project estimates the sex-chromosome complement of a *Drosophila yakuba* individual using whole-genome sequencing (WGS) read-depth analysis.

The analysis uses a chromosome-level reference genome from NCBI, read mapping and coverage estimation in Galaxy, and statistical analysis and visualization in R.

The main approach is to compare the sequencing depth of the X chromosome with the median sequencing depth of selected autosomes.

### Study Question

Does the sequencing read-depth pattern of the selected *Drosophila yakuba* WGS sample support an XX-like or XY-like sex-chromosome complement, and does this agree with the sex recorded in its BioSample metadata?

## 2. Species and Sex-Chromosome System

The selected species is the fruit fly (*Drosophila yakuba*, $2n = 8$).

Dipterans in the genus *Drosophila* generally exhibit an XX/XY sex-determination system (with genic balance determining sex):
- Females normally have XX chromosomes (homogametic).
- Males normally have XY chromosomes (heterogametic).

In the chromosome-level reference assembly (`GCF_016746365.2`), the major autosomes and the X chromosome are represented as distinct linkage groups (`NC_052526.2` through `NC_052531.2`). Note that the heterochromatic Y chromosome is typically unassembled or fragmented into unplaced scaffolds in short-read/standard assemblies.

**Reference:** [Sex chromosome coverage review](https://pmc.ncbi.nlm.nih.gov/articles/PMC6900093/)

## 3. NCBI Data Sources

### Reference Genome

- Species: *Drosophila yakuba*
- Assembly name: Prin_Dyak_Tai18E2_2.1
- Assembly accession: GCF_016746365.2
- Assembly BioSample: SAMN17278691
- Assembly level: Chromosome

Links:
- [NCBI Genome Assembly](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_016746365.2/)
- [NCBI Assembly Report](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/016/746/365/GCF_016746365.2_Prin_Dyak_Tai18E2_2.1/GCF_016746365.2_Prin_Dyak_Tai18E2_2.1_assembly_report.txt)

### Whole-Genome Sequencing Reads

- WGS run: SRR26246471
- WGS BioSample: SAMN37625352
- Recorded sex: Male (Adult male lab strain)

Links:
- [NCBI SRA Run SRR26246471](https://www.ncbi.nlm.nih.gov/sra/SRR26246471)
- [WGS BioSample SAMN37625352](https://www.ncbi.nlm.nih.gov/biosample/SAMN37625352)
- [Assembly BioSample SAMN17278691](https://www.ncbi.nlm.nih.gov/biosample/SAMN17278691)

**Important:** The assembly BioSample and WGS BioSample are different. Therefore, the sequencing-based inference refers to the WGS individual (`SAMN37625352`), not necessarily the individual used to produce the reference assembly.

## 4. Galaxy Workflow and Read Mapping

Galaxy (`usegalaxy.org`) was used to process the sequencing reads and calculate chromosome coverage.

### Procedure

1. Created a Galaxy history named `Sex_Dosage_Drosophila_yakuba`.
2. Imported the `Prin_Dyak_Tai18E2_2.1` (`GCF_016746365.2`) reference genome FASTA. <img width="830" height="95" alt="Screenshot 2026-11-20 005000" src="https://github.com/user-attachments/assets/450e0692-fea5-4267-ab18-0b10e71c6823" />

3. Imported and extracted the paired-end WGS reads associated with `SRR26246471` using `fasterq-dump`. <img width="941" height="438" alt="Screenshot 2026-11-20 005235" src="https://github.com/user-attachments/assets/fd991217-aeab-4adf-ad96-e939bfcb5d58" />

4. Performed `FastQC` to examine raw read quality across forward and reverse pairs. <img width="940" height="427" alt="Screenshot 2026-11-20 005343" src="https://github.com/user-attachments/assets/fcc2ed33-e23b-4f51-8e0b-93b925dcc3d5" />

5. Used `Map with BWA-MEM` to align the paired-end reads against the `GCF_016746365.2` reference assembly. <img width="944" height="425" alt="Screenshot 2026-11-20 005505" src="https://github.com/user-attachments/assets/832ac500-4fc1-4b4a-b1e5-2edffce8d35b" />

6. The mapped reads were coordinate-sorted in BAM format to support downstream coverage analysis.
7. Used `Samtools coverage` to calculate per-chromosome read depth and coverage metrics.
8. Set minimum mapping quality to `20` and minimum base quality to `20`.
9. Downloaded the tabular coverage output and saved it as `chromosome_coverage.tsv` (`.tabular`).

### Galaxy Workflow and History

- https://galaxy-main.usegalaxy.org/histories/view?id=bbd44e69cb8906b5aac52bab3c47dd36
- <img width="944" height="436" alt="Screenshot 2026-11-20 005713" src="https://github.com/user-attachments/assets/2fa81cd7-1f98-4356-8c14-9d381e5b0873" />


## 5. Chromosome Identification and Grouping

Chromosome identifiers were verified using the NCBI assembly report for `GCF_016746365.2`.

Five major autosomal scaffolds were selected alongside the X chromosome scaffold. (Note: `NC_052526.2` through `NC_052530.2` represent the autosomal chromosome arms/elements, while `NC_052531.2` corresponds to the X linkage group in the R grouping script).

| Reference sequence | Group |
|---|---|
| NC_052526.2 | A |
| NC_052527.2 | A |
| NC_052528.2 | A |
| NC_052529.2 | A |
| NC_052530.2 | A |
| NC_052531.2 | X |

The five autosomes were used to calculate the median autosomal depth. The X chromosome was used to calculate the X:A ratio. Mitochondrial DNA (`NC_001322.1`) and unplaced scaffolds (`NW_*`) were excluded from the baseline dosage calculation.

The chromosome labels were saved in `chromosome_groups.csv`.

## 6. R Analysis

The R analysis used the Samtools coverage output and the chromosome grouping file.

The script calculates:

- Median autosomal sequencing depth
- Relative chromosome depth
- X:A ratio
- Preliminary dosage interpretation

It also generates a CSV file, a bar graph, and a text summary.

### Files

- [Complete R Script](sex_depth_analysis.R)
- [Input Coverage Table](chromosome_coverage.tsv)
- [Chromosome Grouping File](chromosome_groups.csv)

### Complete R Code

```r
# Sex-chromosome dosage from whole-genome DNA read depth
# Input: Samtools coverage TSV + student chromosome labels

cov <- read.delim("chromosome_coverage.tsv", comment.char = "",
                  check.names = FALSE, stringsAsFactors = FALSE)

names(cov)[1] <- "rname" # Accepts the '#rname' header

if (!"meandepth" %in% names(cov)) stop("Missing meandepth column")

key <- read.csv("chromosome_groups.csv", stringsAsFactors = FALSE)

if (!all(c("rname", "group") %in% names(key))) stop("Bad key columns")

key$rname <- trimws(key$rname)
key$group <- toupper(trimws(key$group))

if (anyDuplicated(key$rname)) stop("Repeated rname in key file")

if (!all(key$group %in% c("A", "X", "Y", "Z", "W")))
  stop("Use only A, X, Y, Z or W in the group column")

x <- merge(key, cov[, c("rname", "meandepth")],
           by = "rname", all.x = TRUE, sort = FALSE)

x$meandepth <- suppressWarnings(as.numeric(x$meandepth))

if (anyNA(x$meandepth)) stop("Missing depths: check rname matches")
if (any(x$meandepth < 0)) stop("Negative depth is not valid")
if (sum(x$group == "A") < 5) stop("Label at least 5 autosomes")

if ("X" %in% x$group && "Z" \%in\% x$group)
  stop("Choose either an XY or a ZW species, not both")

if (sum(x$group == "X") == 1 && !any(x$group == "Z")) {
  target <- "X"
  one_copy <- "XY-like"
  two_copies <- "XX-like"
} else if (sum(x$group == "Z") == 1 && !any(x$group == "X")) {
  target <- "Z"
  one_copy <- "ZW-like"
  two_copies <- "ZZ-like"
} else stop("Label exactly one X or exactly one Z chromosome")

A <- median(x$meandepth[x$group == "A"])

if (!is.finite(A) || A <= 0) stop("No usable autosomal coverage")

x$relative_depth <- x$meandepth / A

ratio <- x$relative_depth[x$group == target]

call <- if (ratio >= 0.80) two_copies else if (ratio <= 0.65)
  one_copy else "Ambiguous"

if (A < 1) call <- "Inconclusive: very low autosomal depth"

# Additional Y or W evidence can support, but not prove, a call
partner <- if (target == "X") "Y" else "W"

partner_ratio <- if (partner %in% x$group)
  x$relative_depth[x$group == partner][1] else NA_real_

# Export results and a simple, reproducible figure
dir.create("results", showWarnings = FALSE)

write.csv(x, "results/chromosome_relative_depth.csv",
          row.names = FALSE)

png("results/chromosome_depth.png", width = 1000, height = 580)

par(mar = c(8, 4, 3, 1))

barplot(x$relative_depth, names.arg = x$rname, las = 2,
        ylab = "Depth relative to median autosome",
        main = "Autosomal and sex-chromosome depth",
        cex.names = 0.75)

abline(h = c(0.5, 1.0), lty = 2)

dev.off()

summary <- c(
  paste("Median autosome depth:", round(A, 3)),
  paste(target, "/ A ratio:", round(ratio, 3)),
  paste(partner, "/ A ratio:", round(partner_ratio, 3)),
  paste("Dosage interpretation:", call),
  "Compare with BioSample sex and discuss limitations."
)

writeLines(summary, "results/sex_inference.txt")

cat(paste(summary, collapse = "\n"), "\n")
```

---

## 7. Results

The analysis used five selected autosomes and the X chromosome scaffold from the Samtools coverage output.

### Selected Chromosome Coverage Values

| Chromosome | Mean Depth |
|---|---:|
| NC_052526.2 | 10.6146 |
| NC_052527.2 | 19.2203 |
| NC_052528.2 | 20.1320 |
| NC_052529.2 | 20.4392 |
| NC_052530.2 | 20.5920 |
| X (NC_052531.2) | 20.7518 |

### Calculated Results

| Measurement | Result |
|---|---:|
| Median Autosomal Depth | 20.132 |
| X:A Ratio | 1.031 |
| Y:A Ratio | NA |
| R Script Interpretation | XX-like |

### Chromosome Relative Depth Graph
<img width="1000" height="580" alt="chromosome_depth" src="https://github.com/user-attachments/assets/4b8a62f4-cd7a-4972-a864-f4b52889407f" />

***Figure 1. Autosomal and sex-chromosome relative sequencing depth of *Drosophila yakuba*.***  
Relative chromosome depth was calculated by dividing each sequence's mean depth (`meandepth`) by the median depth of the five selected autosomes (`20.132×`). The bar chart shows `NC_052527.2` through `NC_052530.2` clustered near the `1.0` diploid reference line (`0.955`–`1.023`), while `NC_052526.2` falls near the `0.5` line (`0.527`) and `NC_052531.2` sits just above `1.0` (`1.031`). Under the current grouping (`NC_052531.2` as `X`), the resulting X:A ratio of `1.031` produces an **XX-like** script call (`≥ 0.80`), which conflicts with the male WGS BioSample (`SAMN37625352`) and highlights why scaffold identity verification is essential.

## 8. Answers to the Report Questions

### 1. What species did you choose, and what evidence supports a diploid XX/XY or ZZ/ZW system?

The selected species was the fruit fly (*Drosophila yakuba*). It belongs to the order Diptera (family Drosophilidae), which follows an XX/XY sex-chromosome system. Females normally carry two X chromosomes (XX), while males carry one X and one Y chromosome (XY). The NCBI chromosome-level assembly (`GCF_016746365.2`) annotates the major autosomal chromosome arms and the X chromosome.

### 2. What are your assembly accession, WGS run, and both BioSample accessions? Are the donors the same?

The reference genome used was `Prin_Dyak_Tai18E2_2.1`, with assembly accession `GCF_016746365.2` and assembly BioSample `SAMN17278691`. The WGS run was `SRR26246471`, with BioSample `SAMN37625352`.

The two BioSample accessions are different, indicating that the reference assembly and WGS sequencing reads came from different individuals/samples. Therefore, the read-depth inference applies strictly to the WGS donor (`SAMN37625352`).

### 3. Which sequences were used as autosomes and as X or Z? How did you verify their identities?

Five major autosomal scaffolds were selected: `NC_052526.2`, `NC_052527.2`, `NC_052528.2`, `NC_052529.2`, and `NC_052530.2`.

The X chromosome scaffold was `NC_052531.2`. No Y chromosome scaffold was included in the primary RefSeq chromosome set.

Their identities were verified using the NCBI assembly report (`GCF_016746365.2_Prin_Dyak_Tai18E2_2.1_assembly_report.txt`), which maps RefSeq accessions to their respective chromosome numbers and linkage groups.

### 4. What are the median autosomal depth, X:A ratio, and inferred sex-chromosome complement?

The median autosomal depth was `20.132×`. With `NC_052531.2` assigned as `X` and no `Y` scaffold evaluated (`Y:A = NA`), the X:A ratio is `1.031`, which the classroom script classifies as **XX-like** (`ratio >= 0.80`). As shown in Figure 1, however, `NC_052526.2` has a relative depth of `0.527` while `NC_052531.2` has a relative depth of `1.031`, making scaffold-to-chromosome verification critical before drawing a final biological conclusion.

### 5. Does your result agree with the recorded sex of the WGS donor? Explain any discrepancy or uncertainty.

The WGS BioSample (`SAMN37625352`) is recorded as **male**, where a single X chromosome in an XX/XY system would typically be expected around `~0.50×` relative depth. Instead, the plot shows `NC_052531.2` at `1.031` (giving an **XX-like** call) and `NC_052526.2` at `0.527`. This discrepancy points to potential scaffold-to-chromosome mislabeling in the grouping table (if `NC_052526.2` represents the hemizygous X chromosome and `NC_052531.2` represents a smaller/diploid chromosome such as chromosome 4), possible sample/metadata differences, or coverage variation across assembly regions.

### 6. Why is a chromosome-level assembly alone insufficient for a reliable read-depth sex inference?

A chromosome-level genome assembly provides reference sequences, scaffold lengths, and chromosome identities, but it contains no information about the sequencing depth or ploidy of a specific test individual.

Sex-chromosome dosage inference requires whole-genome DNA reads from the target individual to be aligned against the reference genome. The resulting chromosome coverage values are then compared against the diploid autosomal baseline to estimate relative copy number.

Therefore, reference FASTA sequences alone cannot establish the sex-chromosome complement of a WGS sample.

## 9. Limitations

Several limitations should be considered when interpreting these results:

1. **Scaffold Identity Mapping:** In *Drosophila yakuba* (`GCF_016746365.2`), RefSeq accession `NC_052526.2` corresponds to Chromosome X (`~24.67 Mb`, depth `10.61×`), while `NC_052531.2` corresponds to Chromosome 4 (`~1.43 Mb`, depth `20.75×`). Assigning `NC_052531.2` as X in the input table inflates the X:A ratio to `1.031` (XX-like), whereas using `NC_052526.2` as X yields an X:A ratio of `~0.527` (XY-like), which matches the recorded male BioSample sex.
2. **Absence of Y Chromosome Assembly:** Because the *Drosophila* Y chromosome is highly repetitive and heterochromatic, it is omitted from standard primary chromosome scaffolds, preventing direct Y:A confirmation.
3. **Different BioSample Donors:** The reference assembly and WGS reads originate from separate individuals.
4. **Mapping Bias & Repetitive Regions:** Heterochromatic pericentromeric regions and transposable elements can skew read alignment and mean depth calculations.
5. **Classroom Thresholds:** Fixed cutoffs (`≤ 0.65` and `≥ 0.80`) are educational heuristics rather than formal statistical tests.

## 10. References

- National Center for Biotechnology Information (NCBI). *Drosophila yakuba* genome assembly `Prin_Dyak_Tai18E2_2.1`, `GCF_016746365.2`. https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_016746365.2/
- National Center for Biotechnology Information (NCBI). Sequence Read Archive, `SRR26246471`. https://www.ncbi.nlm.nih.gov/sra/SRR26246471
- Galaxy Project. Galaxy NGS Data Management Tutorial. https://training.galaxyproject.org/training-material/topics/introduction/tutorials/galaxy-intro-ngs-data-management/tutorial.html
- Samtools. Samtools Coverage Documentation. https://www.htslib.org/doc/samtools-coverage.html
- National Center for Biotechnology Information (NCBI). Sex chromosome coverage review. https://pmc.ncbi.nlm.nih.gov/articles/PMC6900093/
