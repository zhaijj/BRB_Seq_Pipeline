#!/bin/bash
set -euo pipefail

# ============================================================================
# BRB-Seq Pipeline
# ============================================================================
# Author: Jingjing Zhai
# Date: 2023-10-05
# Last Modified: 2024-Jan-03
# Version: 0.4
# Description: Pipeline for RNA-Seq analysis of 3' RNA-Seq data (BRB-Seq)
# Contact: Jingjing Zhai, jz963@cornell.edu, zhaijingjing603@gmail.com
# ============================================================================

# ============================================================================
# Configuration
# ============================================================================
readonly THREADS=32
readonly THREADS_TRIMMOMATIC=16
readonly THREADS_FEATURECOUNTS=16
readonly THREADS_PARALLEL=24
readonly THREADS_BOWTIE=32

readonly ADAPTER="TruSeq3-SE.fa"
readonly METADATA_FILE="metadata.txt"
readonly KEYFILE="keyFile.txt"

readonly DIR_DEMULTIPLEXED="demultiplexed/"
readonly DIR_TRIMMED="01_trimmed_reads/"
readonly DIR_MAPPING="02_mapping/"
readonly DIR_FEATURECOUNTS="03_featureCounts/"
readonly DIR_SATURATION="04_readSaturation/"
readonly DIR_GENOME_INDEX="genomeIndex/"
readonly DIR_GENOME_ANNOTATION="genome_annotation/"
readonly DIR_RRNA="map2rRNA/"
readonly DIR_BOWTIE_INDEX="bowtie2_index/"

# ============================================================================
# Helper Functions
# ============================================================================

# Print usage information
usage() {
    cat << EOF
Usage: bash $(basename "$0") <command>

Commands:
    index             Generate genome indices
    QC                Quality control and trimming
    mapping           Map reads to genome
    stat              Generate mapping statistics
    featureCounts     Count features
    readSaturation    Analyze read saturation
    gff2bed           Convert GFF3 to BED format
    map2rRNA          Map reads to rRNA

EOF
    exit 1
}

# Logging function with timestamp
log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"
}

# Error logging
log_error() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] ERROR: $*" >&2
}

# Create directory if it doesn't exist
ensure_dir() {
    local dir="$1"
    if [[ ! -d "${dir}" ]]; then
        mkdir -p "${dir}"
        log "Created directory: ${dir}"
    fi
}

# Check if file exists
check_file() {
    local file="$1"
    if [[ ! -f "${file}" ]]; then
        log_error "Required file not found: ${file}"
        exit 1
    fi
}

# Check if species is empty
is_empty_species() {
    local species="$1"
    [[ "${species}" == "EMPTY" ]]
}

# ============================================================================
# Input Validation
# ============================================================================

if [[ $# -ne 1 ]]; then
    usage
fi

readonly COMMAND="$1"

# ============================================================================
# Pipeline Commands
# ============================================================================

# Generate genome indices
run_index() {
    log "Starting genome index generation"
    check_file "${KEYFILE}"
    ensure_dir "${DIR_GENOME_INDEX}"

    while read -r line; do
        local genome_fa=$(echo "${line}" | awk '{print $2}')
        local annotation=$(echo "${line}" | awk '{print $3}')
        local species_name=$(echo "${line}" | awk '{print $1}')

        local species_index_dir="${DIR_GENOME_INDEX}${species_name}"
        ensure_dir "${species_index_dir}"

        log "Generating index for ${species_name}"
        STAR --runThreadN "${THREADS}" \
             --runMode genomeGenerate \
             --genomeDir "${species_index_dir}" \
             --genomeFastaFiles "${genome_fa}" \
             --sjdbGTFfile "${annotation}" \
             --genomeSAindexNbases 13 &
    done < "${KEYFILE}"

    wait
    log "Genome index generation completed"
}

# Quality control and trimming
run_qc() {
    log "Starting quality control and trimming"
    check_file "${METADATA_FILE}"
    ensure_dir "${DIR_TRIMMED}"

    cut -f3 "${METADATA_FILE}" | sed '1d' | while read -r sample; do
        log "Processing sample: ${sample}"

        java -jar /programs/trimmomatic/trimmomatic-0.36.jar SE \
            -threads "${THREADS_TRIMMOMATIC}" \
            -phred33 \
            "${DIR_DEMULTIPLEXED}${sample}.fastq.gz" \
            "${DIR_TRIMMED}${sample}_trimmed.fq" \
            ILLUMINACLIP:"${ADAPTER}":2:30:10 \
            LEADING:3 \
            TRAILING:3 \
            SLIDINGWINDOW:4:20 \
            MINLEN:36

        pigz "${DIR_TRIMMED}${sample}_trimmed.fq"
        fastqc "${DIR_TRIMMED}${sample}_trimmed.fq.gz" -O "${DIR_TRIMMED}"
    done

    log "Quality control completed"
}

# Map reads to genome
run_mapping() {
    log "Starting read mapping"
    check_file "${METADATA_FILE}"
    ensure_dir "${DIR_MAPPING}"

    sed '1d' "${METADATA_FILE}" | while read -r line; do
        local species_name=$(echo "${line}" | awk '{print $1}')
        local plate_pos=$(echo "${line}" | awk '{print $3}')

        if is_empty_species "${species_name}"; then
            log "Skipping ${plate_pos} (EMPTY species)"
            continue
        fi

        local genome_index="${DIR_GENOME_INDEX}${species_name}"
        log "Mapping ${plate_pos} to ${species_name}"

        STAR --readFilesCommand zcat \
             --outFileNamePrefix "${DIR_MAPPING}${plate_pos}" \
             --outSAMtype BAM SortedByCoordinate \
             --outSAMstrandField intronMotif \
             --genomeDir "${genome_index}" \
             --runThreadN "${THREADS}" \
             --readFilesIn "${DIR_TRIMMED}${plate_pos}_trimmed.fq.gz" \
             --twopassMode Basic \
             --limitGenomeGenerateRAM 128000000000
    done

    log "Read mapping completed"
}

# Generate mapping statistics
run_stat() {
    log "Starting statistics generation"
    check_file "${METADATA_FILE}"

    local output_file="summary_statistics.txt"
    echo -e "platePOS\tplateID\trawReads\tcleanReads\tmappedReads\tuniqMappedReads" > "${output_file}"

    sed '1d' "${METADATA_FILE}" | while read -r line; do
        local species_name=$(echo "${line}" | awk '{print $1}')
        local plate_id=$(echo "${line}" | awk '{print $2}')
        local plate_pos=$(echo "${line}" | awk '{print $3}')

        if is_empty_species "${species_name}"; then
            log "Skipping ${plate_pos} (EMPTY species)"
            continue
        fi

        log "Generating statistics for ${plate_pos}"

        local bam_file="${DIR_MAPPING}${plate_pos}Aligned.sortedByCoord.out.bam"
        local log_file="${DIR_MAPPING}${plate_pos}Log.final.out"
        local stats_file="${bam_file}.stats"

        samtools stats "${bam_file}" > "${stats_file}"

        local raw_reads=$(( $(zcat "${DIR_DEMULTIPLEXED}${plate_pos}.fastq.gz" | wc -l) / 4 ))
        local clean_reads=$(grep 'Number of input reads' "${log_file}" | awk -F '|' '{print $2}')
        local uniq_mapped=$(grep 'Uniquely mapped reads number' "${log_file}" | awk -F '|' '{print $2}')
        local mapped_reads=$(grep 'reads mapped:' "${stats_file}" | awk -F ':' '{print $2}')

        echo -e "${plate_pos}\t${plate_id}\t${raw_reads}\t${clean_reads}\t${mapped_reads}\t${uniq_mapped}" >> "${output_file}"
    done

    log "Statistics generation completed"
}

# Count features
run_feature_counts() {
    log "Starting feature counting"
    check_file "${METADATA_FILE}"
    ensure_dir "${DIR_FEATURECOUNTS}"

    sed '1d' "${METADATA_FILE}" | while read -r line; do
        local species_name=$(echo "${line}" | awk '{print $1}')
        local plate_pos=$(echo "${line}" | awk '{print $3}')

        if is_empty_species "${species_name}"; then
            log "Skipping ${plate_pos} (EMPTY species)"
            continue
        fi

        log "Counting features for ${plate_pos}"

        local annotation="${DIR_GENOME_ANNOTATION}${species_name}.gff3"
        local bam_file="${DIR_MAPPING}${plate_pos}Aligned.sortedByCoord.out.bam"
        local output="${DIR_FEATURECOUNTS}${species_name}_${plate_pos}_featureCounts.txt"
        local log_file="${DIR_FEATURECOUNTS}${species_name}_${plate_pos}_featureCounts.log"

        featureCounts \
            --primary \
            -T "${THREADS_FEATURECOUNTS}" \
            -t gene \
            -g ID \
            -a "${annotation}" \
            -o "${output}" \
            "${bam_file}" 2> "${log_file}"
    done

    log "Feature counting completed"
}

# Analyze read saturation
run_read_saturation() {
    log "Starting read saturation analysis"
    check_file "${METADATA_FILE}"
    ensure_dir "${DIR_SATURATION}"

    # Define function for parallel execution
    process_saturation() {
        local line="$1"
        local plate_pos=$(echo "${line}" | awk '{print $3}')
        local species_name=$(echo "${line}" | awk '{print $1}')

        if [[ "${species_name}" == "EMPTY" ]]; then
            return
        fi

        local bam_file="${DIR_MAPPING}${plate_pos}Aligned.sortedByCoord.out.bam"
        local bed_file="${DIR_GENOME_ANNOTATION}${species_name}.bed"
        local output="${DIR_SATURATION}${species_name}_${plate_pos}"

        RPKM_saturation.py -i "${bam_file}" -r "${bed_file}" -o "${output}"
    }

    export -f process_saturation
    export DIR_MAPPING DIR_GENOME_ANNOTATION DIR_SATURATION

    sed '1d' "${METADATA_FILE}" | parallel -j "${THREADS_PARALLEL}" process_saturation {}

    log "Read saturation analysis completed"
}

# Convert GFF3 to BED format
run_gff2bed() {
    log "Starting GFF3 to BED conversion"

    # Note: Requires agat installation
    # conda install -c bioconda agat --experimental-solver=libmamba

    local gff_files="${DIR_GENOME_ANNOTATION}*.gff3"
    if ! compgen -G "${gff_files}" > /dev/null; then
        log_error "No GFF3 files found in ${DIR_GENOME_ANNOTATION}"
        return 1
    fi

    # Convert .gff3 to .bed using parallel processing
    ls "${DIR_GENOME_ANNOTATION}"*.gff3 | \
        parallel -j 4 'agat_convert_sp_gff2bed.pl -gff {} -o {= s/\.gff3$/.bed/ =}'

    log "GFF3 to BED conversion completed"
}

# Map reads to rRNA
run_map2rrna() {
    log "Starting rRNA mapping"
    check_file "${METADATA_FILE}"
    ensure_dir "${DIR_RRNA}"
    ensure_dir "${DIR_BOWTIE_INDEX}"

    local ncrna_file="Ath_ncRNA.fa"
    local ncrna_gz="${ncrna_file}.gz"
    local ncrna_url="https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-58/fasta/arabidopsis_thaliana/ncrna/Arabidopsis_thaliana.TAIR10.ncrna.fa.gz"

    # Download ncRNA reference if not exists
    if [[ ! -f "${ncrna_file}" ]]; then
        log "Downloading ncRNA reference"
        wget -O "${ncrna_gz}" "${ncrna_url}"
        gunzip "${ncrna_gz}"
    fi

    # Build bowtie2 index for rRNA
    local index_prefix="${DIR_BOWTIE_INDEX}rRNA"
    if [[ ! -f "${index_prefix}.1.bt2" ]]; then
        log "Building bowtie2 index for rRNA"
        bowtie2-build "${ncrna_file}" "${index_prefix}"
    fi

    # Map reads to rRNA
    sed '1d' "${METADATA_FILE}" | while read -r line; do
        local plate_pos=$(echo "${line}" | awk '{print $3}')
        log "Mapping ${plate_pos} to rRNA"

        bowtie2 -p "${THREADS_BOWTIE}" \
                -x "${index_prefix}" \
                -U "${DIR_TRIMMED}${plate_pos}_trimmed.fq.gz" \
                -S "${DIR_RRNA}${plate_pos}.sam" \
                2> "${DIR_RRNA}${plate_pos}.log"
    done

    # Summarize alignment rates
    log "Summarizing rRNA alignment rates"
    grep 'overall alignment rate' "${DIR_RRNA}"*.log | \
        sed -e 's/ overall alignment rate//g' \
            -e "s|${DIR_RRNA}||g" \
            -e 's/.log:/\t/' > map2rRNA_percentage.txt

    log "rRNA mapping completed"
}

# ============================================================================
# Main Execution
# ============================================================================

case "${COMMAND}" in
    index)
        run_index
        ;;
    QC)
        run_qc
        ;;
    mapping)
        run_mapping
        ;;
    stat)
        run_stat
        ;;
    featureCounts)
        run_feature_counts
        ;;
    readSaturation)
        run_read_saturation
        ;;
    gff2bed)
        run_gff2bed
        ;;
    map2rRNA)
        run_map2rrna
        ;;
    *)
        log_error "Unknown command: ${COMMAND}"
        usage
        ;;
esac

log "Pipeline step '${COMMAND}' completed successfully"