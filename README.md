# Poaceae RNA Seq #

How to run the BRB-Seq pipeline for multi-species.

### 0. Requirements ###
```
git@github.com:zhaijj/BRB_Seq_Pipeline.git
```

* [Trimmomatic](https://github.com/timflutre/trimmomatic)  
* [STAR](https://github.com/alexdobin/STAR)
* [featureCounts](https://anaconda.org/bioconda/subread)
* [RseQC](https://rseqc.sourceforge.net)


### 1. Build genome index with STAR ###

Prepare a key file for each species, which contains the species name and the path to the reference genome and annotation file. The example key file is [**here**](keyFile.txt).

```
bash BRB_Seq_Pipeline.sh index
```

Then the genome index will be built for each species under current directory.

### 2. Trim reads with Trimmomatic ###
Before this step, make sure you have the metadata file ready. The example metadata file is [**here**](metadata.txt). The metadata file should be a tab-delimited file with the first column as the species name, recommended no space in the name, e.g., Zea_mays instead of Zea mays. The 2nd column should be plate ID, 3rd column should be plate positions, e.g., A01, A02, etc.

Note:
- **raw reads are supposed to be under `demultiplexed` directory where the shell is running**
- **`TruSeq3-SE.fa` is also expected to be under the same directory where the shell is running**

```
bash BRB_Seq_Pipeline.sh QC # using trimmomatic under /programs/trimmomatic/trimmomatic-0.36.jar
```

Then the trimmed reads will be generated for each samples under current `01_trimmed_reads` directory.

### 3. Align reads with STAR ###

```
bash BRB_Seq_Pipeline.sh mapping
```

Then the alignment files will be generated for each samples under current `02_mapping` directory.

### 4. Statistics of alignment ###

```shell
bash BRB_Seq_Pipeline.sh stat
```

Then output will be `summary_statistics.txt`.

### 5. Count reads with featureCounts ###

#### 5.1 Prepare annotation file ####
```bash
if [ ! -d genome_annotation ]; then
    mkdir -p genome_annotation
fi
cd genome_annotation
cat keyFile.txt | while read line
do
    speciesName=$(echo ${line} | awk '{print $1}')
    gff=$(echo ${line} | awk '{print $3}')
    ln -sf ${gff} ${speciesName}.gff3
done
```

#### 5.2 Count reads ####
```shell
bash BRB_Seq_Pipeline.sh featureCounts
```

Then the read counts will be generated for each samples under current `03_readCounts` directory.

### 6. Read saturation analysis ### 

#### 6.1 Prepare BED12 file based on GFF3 ####
```bash
# install agat and RPKM_saturation  first
conda install -c bioconda agat --experimental-solver=libmamba
# check if agat is installed successfully
agat_convert_sp_gff2bed.pl --help
# install RSeQC
pip3 install RSeQC
# check if RPKM_saturation.py is installed successfully
RPKM_saturation.py -h 
# generate BED12 file
bash BRB_Seq_Pipeline.sh gff2bed
```

#### 6.2 Read saturation analysis ####
```bash
bash BRB_Seq_Pipeline.sh readSaturation
```
Then the read saturation analysis will be generated for each samples under current `04_readSaturation` directory. And the example visualization code is [here](readSaturationPlot.R)


### 7. rRNA percentage check (optional) ###

**Make sure bowtie2 is installed**
```bash
bash BRB_Seq_Pipeline.sh map2rRNA
```

Then the rRNA percentage will be recorded in `map2rRNA_percentage.txt` under current directory.
