
# matchTag

## Overview

**matchTag** is an R package that resolves key challenges in marker trait association (MTA) analysis across multiple DArTSeq datasets. It provides tools to map loci tags to reference genomes, identify homologous markers across datasets, and compute allele frequencies for marker development and validation.

## Main Features

### 1. **Genome Mapping**
Map DArTSeq marker tags of SNPs and SIL markers to diferent reference genomes. Usefull to deal with reference bias on association studies. Tag sequences are used to generate a fasta file and BBmap to map them to the input reference genome. Using the tag orientation and CIGAR string is inferred the locus position. 

### 2. **Homologous Loci Discovery**
Find homologous loci across datasets when MTAs are not detected in the original dataset. Using tag sequence similarity, you can identify candidate loci that may represent the same genomic region, even when they have different marker IDs.

### 3. **Allele Frequency Assessment**
Check favorable allelic frequencies for specific loci in different datasets. Identify whether beneficial alleles segregate in other populations and assess their potential for marker assay design.

## Prerequisites

Before using matchTag, ensure you have:

- **Java**: Required to run BBmap for sequence alignment
- **BBtools**: Including the BBmap module for mapping marker tags to reference genomes
  - Download from: https://sourceforge.net/projects/bbmap/

## Installation

You can install the development version of matchTag from GitHub:

``` r
# install.packages("devtools")
devtools::install_github("serbier/matchTag")
```

## Quick Start

### DArTSet Class

The `DArTSet` class represents a DArTSeq dataset (SNP or SIL markers) and provides methods for mapping, homology detection, and allele frequency calculation.

#### Creating a DArTSet Object

```r
library(matchTag)
# load the test DartR gl object
gl_object <- readRDS("tests/testthat/data/gl_snp.RDS")
# Create a DArTSet object from a genlight object
dart_set <- DArTSet$new(
  dataset = gl_object,        # genlight object from dartR
  dataset_type = "SNP",       # or "SIL" for presence/absence
  ref = "tests/testthat/pv21.fa" # use uncompressed fasta file
)
```

#### DArTSet Public Methods

##### `map_tags2ref(bbmap_dir, memory = "4g", use_index = TRUE, threads = 4)`
Maps marker tags to a reference genome using BBmap. Generates mapping alignments and infers locus positions in the reference.

**Parameters:**
- `bbmap_dir`: Directory containing BBmap binaries
- `memory`: Memory allocation for Java (e.g., "4g", "8g")
- `use_index`: Whether to use existing reference index (faster for repeated mappings)
- `threads`: Number of CPU threads to use

**Returns:** Data frame stored in `self$aligments` with mapping metadata

```r
dart_set$map_tags2ref(
  bbmap_dir = "/path/to/bbtools",
  memory = "8g",
  threads = 8
)
```

##### `find_homolog(mta, cutoff = 5)`
Finds homologous loci in the dataset for a given MTA using sequence similarity.

**Parameters:**
- `mta`: An MTA object representing the query marker
- `cutoff`: Maximum Hamming distance to call a homolog (default: 5 bases)

**Returns:** Named numeric vector of distances for matches, or `NA` if no homolog found

```r
# Assuming mta is an MTA object
homologs <- dart_set$find_homolog(mta, cutoff = 3)
print(homologs)  # Shows distances for matching loci
```

##### `get_favourable_freq(alleleID, favorable_allele)`
Retrieves the frequency of a favorable allele for a specific locus in the dataset.

**Parameters:**
- `alleleID`: Marker identifier (must be in `locNames`)
- `favorable_allele`: The favorable allele state ("0", "1" for SIL; allele code for SNP)

**Returns:** Numeric scalar with the favorable allele frequency (0-1)

```r
freq <- dart_set$get_favourable_freq(
  alleleID = "1-1000-A/G", # This locus id should exist in the DArTSet instance
  favorable_allele = "G"
)
print(freq)  # e.g., 0.75
```
### MTA Class

The `MTA` class represents metadata for a marker trait association, including the marker identifier, allele sequence, and favorable allele state.

#### Creating an MTA Object

```r
# SNP MTA
snp_mta <- MTA$new(
    alleleID = "11111111-23-T/A", # Always snp position in tag 23 (0-based)
    alleleSequence = "TGCAGCTTTCAACTCACCAAACATGCGACCAAACTCTTATGTCCAGAAGGGAGCACGTGTTTGGAAGCA",
    favorableAllele = "A"
  )
# SIL (presence/absence) MTA
sil_mta <- MTA$new(
    alleleID = "11111111",
    alleleSequence = "TGCAGGATGGCAGAATGACACCCACCTCTGAACAGCTTGTACTCCAATTAGTTAATACTATTTCTTGTG",
    favorableAllele = 1 # 0 absence | 1 presence
  )
```

#### MTA Public Fields

- `alleleID`: Marker identifier (format: "chr-pos-allele1/allele2" for SNP; "chr-pos" for SIL)
- `alleleSequence`: Full sequence context for the marker
- `favorableAllele`: The beneficial allele state
- `mtaType`: "SNP" or "SIL"
- `alleles`: Vector of possible allele states
- `loc_tag_pos`: Position of the polymorphism in the sequence (SNP only)

## Workflow Example

```r
library(matchTag)
library(dartR)

# Load your genlight object
gl_object <- readRDS("tests/testthat/data/gl_snp.RDS")
# Create DArTSet
dart_set <- DArTSet$new(
  dataset = gl_object,
  dataset_type = "SNP",
  ref = "tests/testthat/pv21.fa"
)

# Map tags to reference
dart_set$map_tags2ref(
  bbmap_dir = "~/software/bbmap",
  memory = "8g",
  threads = 4
)

# Create an MTA for a marker of interest
snp_mta <- MTA$new(
    alleleID = "11111111-23-T/A", # Always snp position in tag 23 (0-based)
    alleleSequence = "TGCAGCTTTCAACTCACCAAACATGCGACCAAACTCTTATGTCCAGAAGGGAGCACGTGTTTGGAAGCA",
    favorableAllele = "A"
  )
# Find homologous loci
homologs <- dart_set$find_homolog(mta, cutoff = 5)
# returns "109760619-23-T/A"

# Check allele frequency in this dataset
freq <- dart_set$get_favourable_freq(
  alleleID = "109760619-23-T/A",
  favorable_allele = "T"
)

cat("Favorable allele frequency:", freq, "\n")
```