# BRB-Seq Pipeline #

A comprehensive pipeline for RNA-Seq analysis of 3' RNA-Seq data (BRB-Seq) supporting multi-species analysis.

## Features ##

- **Elegant and maintainable code structure** with modular functions
- **Robust error handling** with strict mode (`set -euo pipefail`)
- **Timestamped logging** for better debugging and tracking
- **Centralized configuration** for easy customization
- **Consistent code style** with proper quoting and validation

## Requirements ##

### Software Dependencies ###
* [Trimmomatic](https://github.com/timflutre/trimmomatic) - Read trimming and quality control
* [STAR](https://github.com/alexdobin/STAR) - RNA-Seq alignment
* [featureCounts](https://anaconda.org/bioconda/subread) - Read counting
* [RSeQC](https://rseqc.sourceforge.net) - RNA-Seq quality control
* [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) - Quality control reports
* [SAMtools](http://www.htslib.org/) - BAM file manipulation
* [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) - Short read aligner (for rRNA mapping)
* [GNU Parallel](https://www.gnu.org/software/parallel/) - Parallel processing
* [AGAT](https://github.com/NBISweden/AGAT) - GFF/GTF file manipulation

### Installation ###
```bash
git clone git@github.com:zhaijj/BRB_Seq_Pipeline.git
cd BRB_Seq_Pipeline
```

## Usage ##

```bash
bash BRB_Seq_Pipeline.sh <command>

Commands:
    index             Generate genome indices
    QC                Quality control and trimming
    mapping           Map reads to genome
    stat              Generate mapping statistics
    featureCounts     Count features
    readSaturation    Analyze read saturation
    gff2bed           Convert GFF3 to BED format
    map2rRNA          Map reads to rRNA
```


## Pipeline Workflow ##

### 1. Build genome index with STAR ###

Prepare a key file (`keyFile.txt`) for each species, which contains the species name and the path to the reference genome and annotation file. The example key file is [**here**](keyFile.txt).

**Format:** `species_name <tab> genome.fasta <tab> annotation.gtf`

```bash
bash BRB_Seq_Pipeline.sh index
```

**Output:** Genome indices will be created in `genomeIndex/<species_name>/` directory.

### 2. Quality control and trimming with Trimmomatic ###

Before this step, prepare a metadata file (`metadata.txt`). The example metadata file is [**here**](metadata.txt).

**Metadata format:**
- Column 1: Species name (no spaces, e.g., `Zea_mays`)
- Column 2: Plate ID
- Column 3: Plate position (e.g., `A01`, `A02`)

**Prerequisites:**
- Raw reads must be in the `demultiplexed/` directory
- Adapter file `TruSeq3-SE.fa` must be in the working directory

```bash
bash BRB_Seq_Pipeline.sh QC
```

**Output:** Trimmed reads and FastQC reports in `01_trimmed_reads/` directory.

### 3. Map reads with STAR ###

```bash
bash BRB_Seq_Pipeline.sh mapping
```

**Output:** Alignment files (sorted BAM) in `02_mapping/` directory.

**Features:**
- Two-pass mapping for improved accuracy
- Sorted BAM output
- Automatic skipping of EMPTY species entries

### 4. Generate mapping statistics ###

```bash
bash BRB_Seq_Pipeline.sh stat
```

**Output:** `summary_statistics.txt` containing:
- Raw read counts
- Clean read counts (post-trimming)
- Mapped read counts
- Uniquely mapped read counts

### 5. Count features with featureCounts ###

#### 5.1 Prepare annotation files ####

Create symbolic links to your annotation files in the `genome_annotation/` directory:

```bash
mkdir -p genome_annotation
cd genome_annotation

# Create symlinks for each species annotation
while read line; do
    species_name=$(echo "${line}" | awk '{print $1}')
    gff=$(echo "${line}" | awk '{print $3}')
    ln -sf "${gff}" "${species_name}.gff3"
done < ../keyFile.txt

cd ..
```

#### 5.2 Count features ####

```bash
bash BRB_Seq_Pipeline.sh featureCounts
```

**Output:** Feature counts and logs in `03_featureCounts/` directory.

**Features:**
- Counts primary alignments only
- Uses gene-level features with ID attribute
- Generates detailed log files for each sample

### 6. Read saturation analysis ###

#### 6.1 Install dependencies and convert GFF3 to BED ####

```bash
# Install AGAT for GFF3 to BED conversion
conda install -c bioconda agat --experimental-solver=libmamba

# Verify installation
agat_convert_sp_gff2bed.pl --help

# Install RSeQC for saturation analysis
pip3 install RSeQC

# Verify installation
RPKM_saturation.py -h

# Convert GFF3 files to BED12 format
bash BRB_Seq_Pipeline.sh gff2bed
```

**Output:** BED files in `genome_annotation/` directory (one per species).

#### 6.2 Run read saturation analysis ####

```bash
bash BRB_Seq_Pipeline.sh readSaturation
```

**Output:** Saturation analysis results in `04_readSaturation/` directory.

**Features:**
- Uses GNU Parallel for efficient processing of multiple samples
- Generates RPKM saturation curves
- Example visualization code available: [readSaturationPlot.R](readSaturationPlot.R)

### 7. rRNA mapping (optional quality control) ###

Check the percentage of reads mapping to ribosomal RNA to assess sample quality.

```bash
bash BRB_Seq_Pipeline.sh map2rRNA
```

**Output:** `map2rRNA_percentage.txt` with alignment rates for each sample.

**Features:**
- Automatically downloads Arabidopsis thaliana ncRNA reference
- Builds Bowtie2 index (cached for subsequent runs)
- Maps reads to identify rRNA contamination

## Configuration ##

The pipeline uses centralized configuration at the top of `BRB_Seq_Pipeline.sh`. You can adjust:

- **Thread counts** for different tools (STAR, Trimmomatic, featureCounts, Parallel, Bowtie2)
- **Directory paths** for input/output locations
- **Adapter sequences** for trimming

## Troubleshooting ##

- **Error handling:** The pipeline uses strict mode (`set -euo pipefail`) and will exit on errors
- **Logging:** All operations are timestamped for easy tracking and debugging
- **Empty species:** Samples marked as "EMPTY" in the metadata file are automatically skipped

## Citation ##

If you use this pipeline, please cite the relevant tools:
- STAR: Dobin et al. (2013) Bioinformatics
- Trimmomatic: Bolger et al. (2014) Bioinformatics
- featureCounts: Liao et al. (2014) Bioinformatics
- RSeQC: Wang et al. (2012) Bioinformatics
