#!/bin/bash

# Generate CPM-normalized bigWig files from existing CUT&Tag and/or RNA BAMs.
# bigWig files are written next to their source BAM files:
#   CUT&Tag: dedup.bam -> dedup.CPM.bw
#   RNA:     Aligned.sortedByCoord.out.bam -> Aligned.sortedByCoord.out.CPM.bw
#
# Usage:
#   bash raw_data_preprocessing/4_make_bigwig_from_existing_bams.sh [SEARCH_ROOT] [MODE]
#
#   SEARCH_ROOT: directory to search recursively (default: current directory)
#   MODE:        all | cuttag | rna (default: all)

set -euo pipefail

SEARCH_ROOT="${1:-$(pwd)}"
MODE="${2:-all}"

BAMCOVERAGE_BIN="/s2/chengd/miniconda3/bin/bamCoverage"
SAMTOOLS_BIN="/s2/chengd/miniconda3/envs/technical_qc_comparisons/bin/samtools"
BEDTOOLS_BIN="/s2/chengd/miniconda3/envs/technical_qc_comparisons/bin/bedtools"
BLACKLIST="/s1/chengd/blacklist_v1/mm10-blacklist.bed.gz"

THREADS="${THREADS:-18}"
MAPQ_CUTOFF="${MAPQ_CUTOFF:-10}"
BIGWIG_BINSIZE="${BIGWIG_BINSIZE:-10}"

if [ ! -x "${BAMCOVERAGE_BIN}" ]; then
    BAMCOVERAGE_BIN=$(command -v bamCoverage || true)
fi
if [ -z "${BAMCOVERAGE_BIN}" ] || [ ! -x "${BAMCOVERAGE_BIN}" ]; then
    echo "错误：找不到 bamCoverage，请检查 deepTools 环境"
    exit 1
fi

if [ ! -x "${SAMTOOLS_BIN}" ]; then
    SAMTOOLS_BIN=$(command -v samtools || true)
fi
if [ -z "${SAMTOOLS_BIN}" ] || [ ! -x "${SAMTOOLS_BIN}" ]; then
    echo "错误：找不到 samtools"
    exit 1
fi

if [ ! -x "${BEDTOOLS_BIN}" ]; then
    BEDTOOLS_BIN=$(command -v bedtools || true)
fi
if [ -z "${BEDTOOLS_BIN}" ] || [ ! -x "${BEDTOOLS_BIN}" ]; then
    echo "错误：找不到 bedtools"
    exit 1
fi

case "${MODE}" in
    all|cuttag|rna) ;;
    *)
        echo "错误：MODE 必须是 all、cuttag 或 rna"
        exit 1
        ;;
esac

prepare_bamcoverage_blacklist() {
    local out_bed="${SEARCH_ROOT}/mm10-blacklist.merged_for_bamCoverage.bed"
    if [ ! -f "${out_bed}" ]; then
        echo "准备 bamCoverage blacklist: ${out_bed}" >&2
        if [[ "${BLACKLIST}" == *.gz ]]; then
            zcat "${BLACKLIST}"
        else
            cat "${BLACKLIST}"
        fi | awk 'BEGIN{OFS="\t"} NF>=3 {print $1,$2,$3}' | \
            sort -k1,1 -k2,2n | \
            "${BEDTOOLS_BIN}" merge -i - \
            > "${out_bed}"
    fi
    echo "${out_bed}"
}

BAMCOVERAGE_BLACKLIST=$(prepare_bamcoverage_blacklist)

make_index_if_needed() {
    local bam=$1
    local bai="${bam}.bai"
    if [ ! -f "${bai}" ] && [ -f "${bam%.*}.bai" ]; then
        bai="${bam%.*}.bai"
    fi
    if [ ! -f "${bai}" ] || [ "${bam}" -nt "${bai}" ]; then
        echo "  index: ${bam}"
        "${SAMTOOLS_BIN}" index -@ "${THREADS}" "${bam}"
    fi
}

make_cuttag_bw() {
    local bam=$1
    local bam_dir log_dir out_bw
    bam_dir=$(dirname "${bam}")
    log_dir=$(dirname "${bam_dir}")/logs
    out_bw="${bam_dir}/dedup.CPM.bw"
    mkdir -p "${log_dir}"

    if [ -f "${out_bw}" ]; then
        echo "  skip existing: ${out_bw}"
        return 0
    fi

    make_index_if_needed "${bam}"
    echo "  CUT&Tag bigWig: ${out_bw}"
    "${BAMCOVERAGE_BIN}" \
        -b "${bam}" \
        -o "${out_bw}" \
        --outFileFormat bigwig \
        --normalizeUsing CPM \
        --binSize "${BIGWIG_BINSIZE}" \
        --extendReads \
        --samFlagInclude 64 \
        --samFlagExclude 2308 \
        --minMappingQuality "${MAPQ_CUTOFF}" \
        --blackListFileName "${BAMCOVERAGE_BLACKLIST}" \
        --ignoreForNormalization chrM \
        --numberOfProcessors "${THREADS}" \
        2> "${log_dir}/bamCoverage_cuttag.err"
}

make_rna_bw() {
    local bam=$1
    local bam_dir log_dir out_bw
    bam_dir=$(dirname "${bam}")
    log_dir=$(dirname "${bam_dir}")/logs
    out_bw="${bam_dir}/Aligned.sortedByCoord.out.CPM.bw"
    mkdir -p "${log_dir}"

    if [ -f "${out_bw}" ]; then
        echo "  skip existing: ${out_bw}"
        return 0
    fi

    make_index_if_needed "${bam}"
    echo "  RNA bigWig: ${out_bw}"
    "${BAMCOVERAGE_BIN}" \
        -b "${bam}" \
        -o "${out_bw}" \
        --outFileFormat bigwig \
        --normalizeUsing CPM \
        --binSize "${BIGWIG_BINSIZE}" \
        --samFlagExclude 2308 \
        --minMappingQuality "${MAPQ_CUTOFF}" \
        --ignoreForNormalization chrM \
        --numberOfProcessors "${THREADS}" \
        2> "${log_dir}/bamCoverage_rna.err"
}

echo "搜索目录: ${SEARCH_ROOT}"
echo "模式:     ${MODE}"
echo "binSize:  ${BIGWIG_BINSIZE}"
echo "threads:  ${THREADS}"
echo "============================================================================"

if [ "${MODE}" = "all" ] || [ "${MODE}" = "cuttag" ]; then
    echo "--> CUT&Tag BAM -> bigWig"
    while IFS= read -r -d '' bam; do
        make_cuttag_bw "${bam}"
    done < <(find "${SEARCH_ROOT}" -type f -path "*/outs/dedup.bam" -print0 | sort -z)
fi

if [ "${MODE}" = "all" ] || [ "${MODE}" = "rna" ]; then
    echo "--> RNA BAM -> bigWig"
    while IFS= read -r -d '' bam; do
        make_rna_bw "${bam}"
    done < <(find "${SEARCH_ROOT}" -type f -path "*/outs/Aligned.sortedByCoord.out.bam" -print0 | sort -z)
fi

echo "============================================================================"
echo "bigWig 生成完成"
