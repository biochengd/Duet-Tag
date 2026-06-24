#!/bin/bash

# CUT&TAG 数据批量处理脚本（mm10）- MACS2 全宽峰版本
#
# 与 macs_style 版本的差异：
#   - 所有 mark（包括 H3K4me3/H3K27ac/H3K4me1）统一使用宽峰参数
#   - 参数与 nano-CT 论文一致：--llocal 100000 --min-length 1000 --max-gap 1000 --broad-cutoff 0.1
#   - 输出目录后缀：_scMTR（与原版相同，方便对比）
#
# 用法：
#   bash 2_CUTnTAG_mm10_bowtie2_macs_broadonly_style.sh [INPUT_DIR] [OUTPUT_DIR]

set -e

SINTO_THREADS=4   # overcommit_memory=2 下 fork 受限，用 4 线程避免 OOM

# ============================================================================
# 输入/输出目录（支持命令行参数）
# ============================================================================

INPUT_DIR="${1:-$(pwd)}"
OUTPUT_DIR="${2:-$INPUT_DIR}"

echo "输入目录: ${INPUT_DIR}"
echo "输出目录: ${OUTPUT_DIR}"
mkdir -p "${OUTPUT_DIR}"

# ============================================================================
# 全局配置
# ============================================================================

BOWTIE2_BIN="/s2/chengd/miniconda3/envs/cuttag_env/bin/bowtie2"
SAMTOOLS_BIN="/s2/chengd/miniconda3/envs/technical_qc_comparisons/bin/samtools"
PICARD_BIN="/s2/chengd/miniconda3/envs/cuttag_env/bin/picard"
SINTO_BIN="/s2/chengd/miniconda3/envs/cuttag_env/bin/sinto"
BGZIP_BIN="/s2/chengd/miniconda3/envs/technical_qc_comparisons/bin/bgzip"
TABIX_BIN="/s2/chengd/miniconda3/envs/technical_qc_comparisons/bin/tabix"
BEDTOOLS_BIN="/s2/chengd/miniconda3/envs/technical_qc_comparisons/bin/bedtools"
MACS2_BIN="/s2/chengd/miniconda3/envs/technical_qc_comparisons/bin/macs2"
STAR_BIN="/s2/chengd/miniconda3/envs/technical_qc_comparisons/bin/STAR"
PYTHON3_BIN="/s2/chengd/miniconda3/envs/cuttag_env/bin/python3"
BAMCOVERAGE_BIN="/s2/chengd/miniconda3/bin/bamCoverage"

export JAVA_HOME="/s2/chengd/miniconda3/envs/cuttag_env/lib/jvm"
export PATH="${JAVA_HOME}/bin:${PATH}"

GENOME_FA="/s2/SHARE/00_ref_genecode/refdata-cellranger-arc-mm10-2020-A-2.0.0/fasta/genome.fa"
BOWTIE2_IDX="/s1/chengd/referece/mm10_bowtie2/mm10"
BLACKLIST="/s1/chengd/blacklist_v1/mm10-blacklist.bed.gz"
SCRIPTS_DIR="/s1/chengd/project/multiome/01_preprocessing/scripts"

THREADS=8
MAPQ_CUTOFF=10
MACS2_GSIZE="mm"
BIGWIG_BINSIZE=10

# 所有 mark 统一使用的宽峰参数（nano-CT 风格）
MACS2_PARAMS="--llocal 100000 --min-length 1000 --max-gap 1000 --broad --broad-cutoff 0.1 --keep-dup all"

prepare_bamcoverage_blacklist() {
    local OUT_BED="${OUTPUT_DIR}/mm10-blacklist.merged_for_bamCoverage.bed"
    if [ ! -f "${OUT_BED}" ]; then
        echo "准备 bamCoverage blacklist: ${OUT_BED}" >&2
        if [[ "${BLACKLIST}" == *.gz ]]; then
            zcat "${BLACKLIST}"
        else
            cat "${BLACKLIST}"
        fi | awk 'BEGIN{OFS="\t"} NF>=3 {print $1,$2,$3}' | \
            sort -k1,1 -k2,2n | \
            ${BEDTOOLS_BIN} merge -i - \
            > "${OUT_BED}"
    fi
    echo "${OUT_BED}"
}

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

# ============================================================================
# 样本检测
# ============================================================================

SAMPLES=()
for d in "${INPUT_DIR}"/*/; do
    sample=$(basename "$d")
    [[ "$sample" == *_extracted ]] && continue
    if [ -d "${INPUT_DIR}/${sample}_extracted" ]; then
        SAMPLES+=("$sample")
    elif [ -f "${d}${sample}_S1_L001_R1_001.fastq.gz" ]; then
        SAMPLES+=("$sample")
    fi
done

[ ${#SAMPLES[@]} -eq 0 ] && echo "错误：在 ${INPUT_DIR} 中未找到样本" && exit 1

extract_mark() {
    local SAMPLE=$1
    for mark in "H3K27ac" "H3K4me3" "H3K4me1" "H3K27me3" "H3K9me3" "H3K36me3" "ATAC"; do
        [[ "$SAMPLE" == *"$mark"* ]] && echo "$mark" && return 0
    done
    echo ""; return 1
}

echo "============================================================================"
echo "MACS2 全宽峰版本，检测到 ${#SAMPLES[@]} 个样本："
for s in "${SAMPLES[@]}"; do
    mark=$(extract_mark "$s")
    [ -z "$mark" ] && echo "  - $s  [警告：无法识别 MARK]" || echo "  - $s  [MARK: $mark]"
done
echo "============================================================================"

# ============================================================================
# 单样本处理
# ============================================================================

process_sample() {
    local SAMPLE=$1
    local MARK
    MARK=$(extract_mark "$SAMPLE")
    [ -z "$MARK" ] && echo "跳过 ${SAMPLE}：无法识别 MARK" && return 1

    echo ""
    echo "============================================================================"
    echo "处理样本: ${SAMPLE}  MARK: ${MARK}  [MACS2 全宽峰版本]"
    echo "============================================================================"

    local EXTRACT_DIR
    if [ -d "${INPUT_DIR}/${SAMPLE}_extracted" ]; then
        EXTRACT_DIR="${INPUT_DIR}/${SAMPLE}_extracted"
    else
        EXTRACT_DIR="${INPUT_DIR}/${SAMPLE}"
    fi

    local R1_OUT="${EXTRACT_DIR}/${SAMPLE}_S1_L001_R1_001.fastq.gz"
    local R3_OUT="${EXTRACT_DIR}/${SAMPLE}_S1_L001_R3_001.fastq.gz"

    for f in "$R1_OUT" "$R3_OUT"; do
        [ ! -f "$f" ] && echo "错误：$f 不存在" && return 1
    done

    local SAMPLE_OUT="${OUTPUT_DIR}/${SAMPLE}"
    mkdir -p "${SAMPLE_OUT}/outs/raw_mtx" "${SAMPLE_OUT}/outs/filtered_mtx"
    mkdir -p "${SAMPLE_OUT}/macs2_pk" "${SAMPLE_OUT}/logs"

    # ------------------------------------------------------------------
    # 步骤1: 准备 genome.sizes
    # ------------------------------------------------------------------
    if [ ! -f "${SAMPLE_OUT}/genome.sizes" ]; then
        [ ! -f "${GENOME_FA}.fai" ] && ${SAMTOOLS_BIN} faidx ${GENOME_FA}
        cut -f1,2 ${GENOME_FA}.fai | sort -k1,1 -k2,2n > ${SAMPLE_OUT}/genome.sizes
    fi

    # ------------------------------------------------------------------
    # 步骤2: Bowtie2 比对
    # ------------------------------------------------------------------
    echo "--> 步骤2: Bowtie2 比对（MAPQ > ${MAPQ_CUTOFF}）..."
    local BAM_SORTED="${SAMPLE_OUT}/outs/sorted.bam"
    local BAM_DEDUP="${SAMPLE_OUT}/outs/dedup.bam"

    if [ ! -f "${BAM_DEDUP}" ]; then
        local R1_BC="${SAMPLE_OUT}/outs/R1_bc.fastq.gz"
        local R3_BC="${SAMPLE_OUT}/outs/R3_bc.fastq.gz"
        local R2_BC_FILE="${EXTRACT_DIR}/${SAMPLE}_S1_L001_R2_001.fastq.gz"

        if [ ! -f "${R1_BC}" ] && [ ! -f "${BAM_SORTED}" ]; then
            echo "    写入 barcode 到 read name..."
            paste <(zcat ${R2_BC_FILE}) <(zcat ${R1_OUT}) | \
                awk '{
                    if (NR%4==1) { r1_name=$2 }
                    else if (NR%4==2) { bc=$1; seq=$2 }
                    else if (NR%4==3) { }
                    else if (NR%4==0) {
                        split(r1_name,a,"@")
                        print "@"bc":"a[2]
                        print seq
                        print "+"
                        print $2
                    }
                }' | gzip > ${R1_BC} &
            paste <(zcat ${R2_BC_FILE}) <(zcat ${R3_OUT}) | \
                awk '{
                    if (NR%4==1) { r1_name=$2 }
                    else if (NR%4==2) { bc=$1; seq=$2 }
                    else if (NR%4==3) { }
                    else if (NR%4==0) {
                        split(r1_name,a,"@")
                        print "@"bc":"a[2]
                        print seq
                        print "+"
                        print $2
                    }
                }' | gzip > ${R3_BC} &
            wait
        fi

        if [ ! -f "${BAM_SORTED}" ]; then
            ${BOWTIE2_BIN} \
                -p ${THREADS} \
                --very-sensitive --end-to-end \
                --no-mixed --no-discordant --no-unal \
                -I 10 -X 700 \
                --rg-id "${SAMPLE}" \
                --rg "SM:${SAMPLE}" \
                --rg "PL:ILLUMINA" \
                -x ${BOWTIE2_IDX} \
                -1 ${R1_BC} -2 ${R3_BC} \
                2> ${SAMPLE_OUT}/logs/bowtie2.err | \
            ${SAMTOOLS_BIN} view -bS -q ${MAPQ_CUTOFF} -@ 4 - | \
            ${SAMTOOLS_BIN} sort -@ 8 -o ${BAM_SORTED}
            ${SAMTOOLS_BIN} index ${BAM_SORTED}
            rm -f ${R1_BC} ${R3_BC}
        fi

        mkdir -p ${SAMPLE_OUT}/outs/tmp
        ${PICARD_BIN} MarkDuplicates \
            -I ${BAM_SORTED} \
            -O ${BAM_DEDUP} \
            -M ${SAMPLE_OUT}/logs/markdup_metrics.txt \
            --REMOVE_DUPLICATES true \
            --TMP_DIR ${SAMPLE_OUT}/outs/tmp \
            2> ${SAMPLE_OUT}/logs/markdup.err

        ${SAMTOOLS_BIN} index -@ 8 ${BAM_DEDUP}
        rm -f ${BAM_SORTED} ${BAM_SORTED}.bai
        rm -rf ${SAMPLE_OUT}/outs/tmp
    fi
    echo "  比对完成"

    # ------------------------------------------------------------------
    # 步骤2b: BAM -> bigWig（CPM normalized，保存到 BAM 所在目录）
    # ------------------------------------------------------------------
    echo "--> 步骤2b: 生成 CUT&Tag bigWig..."
    local CUTTAG_BW="${SAMPLE_OUT}/outs/dedup.CPM.bw"
    if [ ! -f "${CUTTAG_BW}" ]; then
        index_bam_if_needed "${BAM_DEDUP}"
        if [ ! -x "${BAMCOVERAGE_BIN}" ]; then
            BAMCOVERAGE_BIN=$(command -v bamCoverage || true)
        fi
        if [ -z "${BAMCOVERAGE_BIN}" ] || [ ! -x "${BAMCOVERAGE_BIN}" ]; then
            echo "错误：找不到 bamCoverage，请检查 deepTools 环境" && return 1
        fi
        local BW_BLACKLIST
        BW_BLACKLIST=$(prepare_bamcoverage_blacklist)
        ${BAMCOVERAGE_BIN} \
            -b "${BAM_DEDUP}" \
            -o "${CUTTAG_BW}" \
            --outFileFormat bigwig \
            --normalizeUsing CPM \
            --binSize ${BIGWIG_BINSIZE} \
            --extendReads \
            --samFlagInclude 64 \
            --samFlagExclude 2308 \
            --minMappingQuality ${MAPQ_CUTOFF} \
            --blackListFileName "${BW_BLACKLIST}" \
            --ignoreForNormalization chrM \
            --numberOfProcessors ${THREADS} \
            2> ${SAMPLE_OUT}/logs/bamCoverage_cuttag.err
    fi
    echo "  bigWig: ${CUTTAG_BW}"

    # ------------------------------------------------------------------
    # 步骤3: BAM → fragments.tsv
    # ------------------------------------------------------------------
    echo "--> 步骤3: 生成 fragments..."
    if [ ! -f "${SAMPLE_OUT}/outs/fragments.tsv.gz" ]; then
        ${PYTHON3_BIN} - ${BAM_DEDUP} ${SAMPLE_OUT}/outs/dedup_cb.bam <<'PYEOF'
import pysam, sys
in_bam, out_bam = sys.argv[1], sys.argv[2]
with pysam.AlignmentFile(in_bam, "rb") as fin, \
     pysam.AlignmentFile(out_bam, "wb", header=fin.header) as fout:
    for read in fin:
        bc = read.query_name.split(":")[0]
        if len(bc) >= 8:
            read.set_tag("CB", bc)
        fout.write(read)
PYEOF
        ${SAMTOOLS_BIN} index ${SAMPLE_OUT}/outs/dedup_cb.bam

        ${SINTO_BIN} fragments \
            -b ${SAMPLE_OUT}/outs/dedup_cb.bam \
            -f ${SAMPLE_OUT}/outs/fragments.tsv \
            -t CB \
            -p ${SINTO_THREADS} \
            --min_mapq ${MAPQ_CUTOFF}

        rm -f ${SAMPLE_OUT}/outs/dedup_cb.bam ${SAMPLE_OUT}/outs/dedup_cb.bam.bai

        sort -k1,1 -k2,2n ${SAMPLE_OUT}/outs/fragments.tsv \
            > ${SAMPLE_OUT}/outs/fragments.sorted.tsv
        mv ${SAMPLE_OUT}/outs/fragments.sorted.tsv ${SAMPLE_OUT}/outs/fragments.tsv
        ${BGZIP_BIN} ${SAMPLE_OUT}/outs/fragments.tsv
        ${TABIX_BIN} -p bed ${SAMPLE_OUT}/outs/fragments.tsv.gz
    fi

    # ------------------------------------------------------------------
    # 步骤4: MACS2 Peak Calling（所有 mark 统一宽峰参数）
    # ------------------------------------------------------------------
    echo "--> 步骤4: MACS2 Peak Calling（全宽峰，${MARK}）..."
    local PEAK_FILE="${SAMPLE_OUT}/macs2_pk/aggregate_peaks.broadPeak"

    if [ ! -f "${PEAK_FILE}" ]; then
        ${MACS2_BIN} callpeak \
            -t ${BAM_DEDUP} \
            -g ${MACS2_GSIZE} -f BAMPE \
            ${MACS2_PARAMS} \
            -B --SPMR \
            --outdir ${SAMPLE_OUT}/macs2_pk \
            -n aggregate \
            2> ${SAMPLE_OUT}/logs/macs2.err
    fi

    local FORMATTED_PEAK="${SAMPLE_OUT}/macs2_pk/aggregate_peaks_formatted.bed"
    if [ ! -f "${FORMATTED_PEAK}" ]; then
        cut -f1-4 ${PEAK_FILE} | sed '/chrM/d' | \
            awk '$1 ~ /^chr([0-9]+|[XYM])$/' | \
            ${BEDTOOLS_BIN} intersect -a - -b ${BLACKLIST} -v | \
            sort -k1,1 -k2,2n > ${SAMPLE_OUT}/macs2_pk/temp_filtered.bed

        ${BEDTOOLS_BIN} merge -i ${SAMPLE_OUT}/macs2_pk/temp_filtered.bed -d 1000 | \
            awk 'BEGIN{OFS="\t"} {print $1,$2,$3,"peak_"NR}' \
            > ${FORMATTED_PEAK}

        rm ${SAMPLE_OUT}/macs2_pk/temp_filtered.bed
    fi
    echo "  MACS2 完成，peaks: $(wc -l < ${FORMATTED_PEAK})"

    PEAK_FILE="${FORMATTED_PEAK}"

    # ------------------------------------------------------------------
    # 步骤5: peak 矩阵生成
    # ------------------------------------------------------------------
    echo "--> 步骤5: 生成 peak 矩阵..."
    if [ ! -f "${SAMPLE_OUT}/outs/peak_read_ov.tsv.gz" ]; then
        zcat ${SAMPLE_OUT}/outs/fragments.tsv.gz | sed '/^#/d' | \
            awk 'BEGIN{OFS="\t"} $1 ~ /^chr([0-9]+|[XYM])$/ {print $1,$2,$3,$4}' | \
            sort -k1,1 -k2,2n | \
            ${BEDTOOLS_BIN} intersect \
                -a ${PEAK_FILE} \
                -b - \
                -wo -sorted -g ${SAMPLE_OUT}/genome.sizes | \
            sort -k8,8 | \
            ${BEDTOOLS_BIN} groupby -g 8 -c 4 -o freqdesc | \
            gzip > ${SAMPLE_OUT}/outs/peak_read_ov.tsv.gz
    fi

    if [ ! -f "${SAMPLE_OUT}/outs/raw_mtx/barcodes.tsv" ]; then
        cut -f1-3 ${PEAK_FILE} > ${SAMPLE_OUT}/outs/raw_mtx/peaks.bed
        zcat ${SAMPLE_OUT}/outs/peak_read_ov.tsv.gz | cut -f1 \
            > ${SAMPLE_OUT}/outs/raw_mtx/barcodes.tsv

        ${PYTHON3_BIN} ${SCRIPTS_DIR}/generate_csc_mtx.py \
            ${PEAK_FILE} \
            ${SAMPLE_OUT}/outs/raw_mtx/barcodes.tsv \
            ${SAMPLE_OUT}/outs/peak_read_ov.tsv.gz \
            ${SAMPLE_OUT}/outs/raw_mtx/matrix.mtx

        awk 'BEGIN{OFS="\t"} {print $1"-"$2"-"$3,$1"-"$2"-"$3}' \
            ${SAMPLE_OUT}/outs/raw_mtx/peaks.bed \
            > ${SAMPLE_OUT}/outs/raw_mtx/features.tsv
    fi

    # ------------------------------------------------------------------
    # 步骤6: 5kb bin 矩阵
    # ------------------------------------------------------------------
    echo "--> 步骤6: 生成 5kb bin 矩阵..."
    local BIN_DIR="${SAMPLE_OUT}/outs/bin5k_mtx"
    mkdir -p "${BIN_DIR}"

    if [ ! -f "${BIN_DIR}/barcodes.tsv" ]; then
        if [ ! -f "${SAMPLE_OUT}/bins_5k.bed" ]; then
            ${BEDTOOLS_BIN} makewindows \
                -g ${SAMPLE_OUT}/genome.sizes \
                -w 5000 | \
                awk '$1 ~ /^chr([0-9]+|[XYM])$/' | \
                grep -v chrM | \
                awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$1"-"$2"-"$3}' \
                > ${SAMPLE_OUT}/bins_5k.bed
        fi

        zcat ${SAMPLE_OUT}/outs/fragments.tsv.gz | sed '/^#/d' | \
            awk 'BEGIN{OFS="\t"} $1 ~ /^chr([0-9]+|[XYM])$/ {print $1,$2,$3,$4}' | \
            sort -k1,1 -k2,2n | \
            ${BEDTOOLS_BIN} intersect \
                -a ${SAMPLE_OUT}/bins_5k.bed \
                -b - \
                -wo -sorted -g ${SAMPLE_OUT}/genome.sizes | \
            sort -k8,8 | \
            ${BEDTOOLS_BIN} groupby -g 8 -c 4 -o freqdesc | \
            gzip > ${SAMPLE_OUT}/outs/bin5k_read_ov.tsv.gz

        cut -f1-3 ${SAMPLE_OUT}/bins_5k.bed > ${BIN_DIR}/bins.bed
        zcat ${SAMPLE_OUT}/outs/bin5k_read_ov.tsv.gz | cut -f1 \
            > ${BIN_DIR}/barcodes.tsv

        ${PYTHON3_BIN} ${SCRIPTS_DIR}/generate_csc_mtx.py \
            ${SAMPLE_OUT}/bins_5k.bed \
            ${BIN_DIR}/barcodes.tsv \
            ${SAMPLE_OUT}/outs/bin5k_read_ov.tsv.gz \
            ${BIN_DIR}/matrix.mtx

        awk 'BEGIN{OFS="\t"} {print $1"-"$2"-"$3,$1"-"$2"-"$3}' \
            ${BIN_DIR}/bins.bed \
            > ${BIN_DIR}/features.tsv
    fi

    # ------------------------------------------------------------------
    # 步骤7: STAR EmptyDrops 过滤
    # ------------------------------------------------------------------
    echo "--> 步骤7: STAR EmptyDrops 过滤..."
    if [ ! -f "${SAMPLE_OUT}/outs/filtered_mtx/barcodes.tsv" ]; then
        ${STAR_BIN} --runMode soloCellFiltering \
            ${SAMPLE_OUT}/outs/raw_mtx \
            ${SAMPLE_OUT}/outs/filtered_mtx/ \
            --soloCellFilter EmptyDrops_CR 3000 0.99 10 45000 90000 500 0.01 20000 0.001 10000 1
    fi

    # ------------------------------------------------------------------
    # 步骤8: QC 指标
    # ------------------------------------------------------------------
    echo "--> 步骤8: QC 指标..."
    if [ ! -f "${SAMPLE_OUT}/outs/fragment_size_distribution.tsv" ]; then
        echo -e "size\tcount" > ${SAMPLE_OUT}/outs/fragment_size_distribution.tsv
        zcat ${SAMPLE_OUT}/outs/fragments.tsv.gz | sed '/^#/d' | \
            awk '{print $3-$2}' | sort -n | uniq -c | sort -b -k2,2n | \
            awk 'BEGIN{OFS="\t"} {print $2,$1}' \
            >> ${SAMPLE_OUT}/outs/fragment_size_distribution.tsv
    fi

    if [ ! -f "${SAMPLE_OUT}/outs/frip.txt" ]; then
        local TOTAL_FRAGS=$(zcat ${SAMPLE_OUT}/outs/fragments.tsv.gz | sed '/^#/d' | \
            awk '$1 ~ /^chr([0-9]+|[XYM])$/' | wc -l)
        local IN_PEAKS=$(zcat ${SAMPLE_OUT}/outs/fragments.tsv.gz | sed '/^#/d' | \
            awk 'BEGIN{OFS="\t"} $1 ~ /^chr([0-9]+|[XYM])$/ {print $1,$2,$3}' | \
            ${BEDTOOLS_BIN} intersect -a - -b ${PEAK_FILE} -u | wc -l)
        local FRIP=$(echo "scale=4; ${IN_PEAKS} / ${TOTAL_FRAGS}" | bc)
        echo -e "total_fragments\t${TOTAL_FRAGS}" > ${SAMPLE_OUT}/outs/frip.txt
        echo -e "fragments_in_peaks\t${IN_PEAKS}" >> ${SAMPLE_OUT}/outs/frip.txt
        echo -e "FRiP\t${FRIP}" >> ${SAMPLE_OUT}/outs/frip.txt
        echo "  FRiP: ${FRIP} (${IN_PEAKS}/${TOTAL_FRAGS})"
    fi

    if [ ! -f "${SAMPLE_OUT}/outs/filtered_mtx/metrics.csv" ]; then
        ${PYTHON3_BIN} ${SCRIPTS_DIR}/output_qc_mtx.py \
            ${SAMPLE_OUT}/outs/filtered_mtx/matrix.mtx \
            ${SAMPLE_OUT}/outs/filtered_mtx/barcodes.tsv \
            ${SAMPLE_OUT}/outs/filtered_mtx/metrics.csv
    fi

    echo ""
    echo "=== 完成: ${SAMPLE} - ${MARK} ==="
    echo "  输出目录:        ${SAMPLE_OUT}/"
    echo "  peaks:           ${PEAK_FILE}"
    echo "  fragments:       ${SAMPLE_OUT}/outs/fragments.tsv.gz"
    echo "  filtered matrix: ${SAMPLE_OUT}/outs/filtered_mtx/"
    echo "  BAM:             ${SAMPLE_OUT}/outs/dedup.bam"
    echo "  bigWig:          ${SAMPLE_OUT}/outs/dedup.CPM.bw"
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
