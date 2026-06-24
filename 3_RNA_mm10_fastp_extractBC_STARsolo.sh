#!/bin/bash

# RNA 数据批量处理脚本（mm10）
# 输入：已提取 barcode 的样本目录（R1=barcode+UMI，R2=cDNA）
# 流程：STARsolo 比对定量 → BAM 索引 → QC 指标
#
# 用法：
#   bash 3_RNA_mm10_fastp_extractBC_STARsolo.sh [INPUT_DIR] [OUTPUT_DIR]
#
#   INPUT_DIR:  包含样本子目录的根目录，每个子目录下有 *_R1*.fastq.gz 和 *_R2*.fastq.gz
#   OUTPUT_DIR: 结果输出目录
#
# 示例：
#   bash 3_RNA_mm10_fastp_extractBC_STARsolo.sh \
#       /nas31/ypj/figure2/Hw_mix/RNA/v17/2.extract_bc \
#       /s1/chengd/project/multimodal_cuttag_rna/part3_hippocampus/data/v17_260109DDC_hw

set -e

# ============================================================================
# 输入/输出目录
# ============================================================================

INPUT_DIR="${1:-$(pwd)}"
OUTPUT_DIR="${2:-$(pwd)}"

echo "输入目录: ${INPUT_DIR}"
echo "输出目录: ${OUTPUT_DIR}"
mkdir -p "${OUTPUT_DIR}"

# ============================================================================
# 全局配置
# ============================================================================

STAR_BIN="/s2/chengd/miniconda3/envs/technical_qc_comparisons/bin/STAR"
SAMTOOLS_BIN="/s2/chengd/miniconda3/envs/technical_qc_comparisons/bin/samtools"
BAMCOVERAGE_BIN="/s2/chengd/miniconda3/bin/bamCoverage"
SCRIPTS_DIR="/s1/chengd/project/multiome/01_preprocessing/scripts"

STAR_IDX="/s1/chengd/referece/star_2.7.11b_mm10"
WHITELIST="/s1/chengd/referece/barcode/737K-cratac-v1.txt"

THREADS=40
BIGWIG_BINSIZE=10

index_bam_if_needed() {
    local BAM=$1
    local BAI="${BAM}.bai"
    if [ ! -f "${BAI}" ] && [ -f "${BAM%.*}.bai" ]; then
        BAI="${BAM%.*}.bai"
    fi
    if [ ! -f "${BAI}" ] || [ "${BAM}" -nt "${BAI}" ]; then
        ${SAMTOOLS_BIN} index -@ ${THREADS} "${BAM}"
    fi
}

# STAR 参数（R1 结构: 16bp Barcode + 2bp GG + 10bp UMI）
SOLOTYPE="CB_UMI_Simple"
CB_START=1
CB_LEN=16
UMI_START=19
UMI_LEN=10
CLIP3P=0

# ============================================================================
# 样本检测
# ============================================================================

SAMPLES=()
for d in "${INPUT_DIR}"/*/; do
    sample=$(basename "$d")
    if ls "${d}"*_R1*.f*q.gz 2>/dev/null | head -1 | grep -q .; then
        SAMPLES+=("$sample")
    fi
done

[ ${#SAMPLES[@]} -eq 0 ] && echo "错误：在 ${INPUT_DIR} 中未找到样本" && exit 1

echo "============================================================================"
echo "检测到 ${#SAMPLES[@]} 个样本："
for s in "${SAMPLES[@]}"; do echo "  - $s"; done
echo "============================================================================"

# ============================================================================
# 单样本处理
# ============================================================================

process_sample() {
    local SAMPLE=$1
    local SAMPLE_IN="${INPUT_DIR}/${SAMPLE}"
    local SAMPLE_OUT="${OUTPUT_DIR}/${SAMPLE}_RNA"

    echo ""
    echo "============================================================================"
    echo "处理样本: ${SAMPLE}"
    echo "============================================================================"

    # 输入文件（R1=barcode+UMI，R2=cDNA）
    local R1 R2
    R1=$(ls "${SAMPLE_IN}"/*_R1*.f*q.gz 2>/dev/null | head -1)
    R2=$(ls "${SAMPLE_IN}"/*_R2*.f*q.gz 2>/dev/null | head -1)

    for f in "$R1" "$R2"; do
        [ ! -f "$f" ] && echo "错误：输入文件不存在: $f" && return 1
    done

    mkdir -p "${SAMPLE_OUT}/outs" "${SAMPLE_OUT}/logs"

    # ------------------------------------------------------------------
    # 步骤1: STARsolo 比对与定量
    # ------------------------------------------------------------------
    echo "--> 步骤1: STARsolo 比对..."
    if [ ! -f "${SAMPLE_OUT}/outs/Aligned.sortedByCoord.out.bam" ]; then
        ${STAR_BIN} \
            --genomeDir ${STAR_IDX} \
            --outFileNamePrefix ${SAMPLE_OUT}/outs/ \
            --readFilesCommand zcat \
            --readFilesIn ${R2} ${R1} \
            --soloType ${SOLOTYPE} \
            --soloCBstart ${CB_START} \
            --soloCBlen ${CB_LEN} \
            --soloUMIstart ${UMI_START} \
            --soloUMIlen ${UMI_LEN} \
            --clip3pNbases ${CLIP3P} \
            --soloCBwhitelist ${WHITELIST} \
            --soloCellFilter EmptyDrops_CR 10000 0.99 10 45000 90000 500 0.01 20000 0.01 10000 \
            --runThreadN ${THREADS} \
            --outSAMattributes CB UB \
            --outSAMtype BAM SortedByCoordinate \
            --soloFeatures Gene GeneFull \
            --soloMultiMappers EM \
            --soloBarcodeReadLength 0 \
            2> ${SAMPLE_OUT}/logs/star.err
    fi

    # ------------------------------------------------------------------
    # 步骤2: BAM 索引
    # ------------------------------------------------------------------
    echo "--> 步骤2: BAM 索引..."
    if [ -f "${SAMPLE_OUT}/outs/Aligned.sortedByCoord.out.bam" ] && \
       [ ! -f "${SAMPLE_OUT}/outs/Aligned.sortedByCoord.out.bam.bai" ]; then
        ${SAMTOOLS_BIN} index ${SAMPLE_OUT}/outs/Aligned.sortedByCoord.out.bam
    fi

    # ------------------------------------------------------------------
    # 步骤2b: RNA BAM -> bigWig（CPM normalized，保存到 BAM 所在目录）
    # ------------------------------------------------------------------
    echo "--> 步骤2b: 生成 RNA bigWig..."
    local RNA_BAM="${SAMPLE_OUT}/outs/Aligned.sortedByCoord.out.bam"
    local RNA_BW="${SAMPLE_OUT}/outs/Aligned.sortedByCoord.out.CPM.bw"
    if [ -f "${RNA_BAM}" ] && [ ! -f "${RNA_BW}" ]; then
        index_bam_if_needed "${RNA_BAM}"
        if [ ! -x "${BAMCOVERAGE_BIN}" ]; then
            BAMCOVERAGE_BIN=$(command -v bamCoverage || true)
        fi
        if [ -z "${BAMCOVERAGE_BIN}" ] || [ ! -x "${BAMCOVERAGE_BIN}" ]; then
            echo "错误：找不到 bamCoverage，请检查 deepTools 环境" && return 1
        fi
        ${BAMCOVERAGE_BIN} \
            -b ${RNA_BAM} \
            -o ${RNA_BW} \
            --outFileFormat bigwig \
            --normalizeUsing CPM \
            --binSize ${BIGWIG_BINSIZE} \
            --samFlagExclude 2308 \
            --minMappingQuality 10 \
            --ignoreForNormalization chrM \
            --numberOfProcessors ${THREADS} \
            2> ${SAMPLE_OUT}/logs/bamCoverage_rna.err
    fi
    echo "  bigWig: ${RNA_BW}"

    # ------------------------------------------------------------------
    # 步骤3: QC 指标
    # ------------------------------------------------------------------
    echo "--> 步骤3: QC 指标..."
    local MATRIX_FILE="${SAMPLE_OUT}/outs/Solo.out/GeneFull/filtered/matrix.mtx"
    local BARCODES_FILE="${SAMPLE_OUT}/outs/Solo.out/GeneFull/filtered/barcodes.tsv"
    local METRICS_FILE="${SAMPLE_OUT}/outs/Solo.out/GeneFull/filtered/metrics.csv"

    if [ ! -f "${METRICS_FILE}" ] && [ -f "${MATRIX_FILE}" ]; then
        python3 ${SCRIPTS_DIR}/output_qc_mtx.py \
            ${MATRIX_FILE} ${BARCODES_FILE} ${METRICS_FILE}
    fi

    rm -rf ${SAMPLE_OUT}/outs/_STARtmp/ 2>/dev/null || true

    echo ""
    echo "=== 完成: ${SAMPLE} ==="
    echo "  filtered matrix: ${SAMPLE_OUT}/outs/Solo.out/GeneFull/filtered/"
    echo "  BAM:             ${SAMPLE_OUT}/outs/Aligned.sortedByCoord.out.bam"
    echo "  bigWig:          ${SAMPLE_OUT}/outs/Aligned.sortedByCoord.out.CPM.bw"
}

# ============================================================================
# 主循环
# ============================================================================

FAILED=()
for SAMPLE in "${SAMPLES[@]}"; do
    if process_sample "$SAMPLE"; then
        echo "✓ ${SAMPLE} 处理成功"
    else
        echo "✗ ${SAMPLE} 处理失败"
        FAILED+=("$SAMPLE")
    fi
done

echo ""
echo "============================================================================"
echo "全部完成！成功: $((${#SAMPLES[@]} - ${#FAILED[@]})) / ${#SAMPLES[@]}"
[ ${#FAILED[@]} -gt 0 ] && echo "失败：${FAILED[*]}"
echo "============================================================================"
