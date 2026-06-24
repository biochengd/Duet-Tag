#!/usr/bin/env python3
"""
CUT&TAG 多抗体样本拆分脚本
流程：
  1. 根据 R2 中的 index 序列将混合测序数据拆分为各抗体独立文件
  2. 对每个抗体文件执行 fastp 质控 + barcode 提取 + 合并
     输出已提取好 barcode 的 R1/R2/R3/I1 文件，供脚本2直接使用

内置 index → 抗体 映射表：
  TACAC → H3K27ac
  CTTGA → H3K27me3
  AGCTA → H3K4me3
  TGAGT → H3K4me1
  ACTGC → H3K36me3
  GATCA → H3K9me3
  AGATG → ATAC

输出目录结构：
  demultiplex_out/
    {prefix}_{mark}_R1.fq.gz          # demultiplex 原始输出
    {prefix}_{mark}_R2.fq.gz
    {prefix}_{mark}_extracted/         # barcode 提取后的文件（供脚本2使用）
      {prefix}_{mark}_S1_L001_R1_001.fastq.gz   # reads
      {prefix}_{mark}_S1_L001_R2_001.fastq.gz   # barcode
      {prefix}_{mark}_S1_L001_R3_001.fastq.gz   # reads (另一端)
      {prefix}_{mark}_S1_L001_I1_001.fastq.gz

使用方式：
  在原始数据目录下直接运行：
  python 1_demultiplex_and_barcodeextract.py
"""

import os
import time
import subprocess
import multiprocessing as mp
from collections import defaultdict

# ============================================================================
# 工具路径
# ============================================================================
FASTP_BIN      = "/s2/Liuyang/biosoft/fastp"
EXTRACT_BC_BIN = "/s2/SHARE/reference/scATAC/extract_bc3_right_ATAC"
RUSH_BIN       = "/s1/Liuyang/biosoft/rush"
PIGZ_BIN       = "/s1/Liuyang/biosoft/pigz-2.8/pigz"

# fastp 分片数（并行 barcode 提取）
SPLIT_NUM = 60

# ============================================================================
# 内置 index → 抗体 映射表
# ============================================================================
INDEX_TO_MARK = {
    "TACAC": "H3K27ac",
    "CTTGA": "H3K27me3",
    "AGCTA": "H3K4me3",
    "TGAGT": "H3K4me1",
    "ACTGC": "H3K36me3",
    "GATCA": "H3K9me3",
    "AGATG": "ATAC",
}

# ============================================================================
# 算法参数
# ============================================================================
FIXED_SEQUENCE       = "GTCTCGTGGGCTCGG"
FIXED_SEQ_MAX_ERR    = 3
INDEX_MAX_ERR        = 1
SEARCH_START         = 20
SEARCH_END           = 40
ADAPTER_REMOVE_START = 24
ADAPTER_REMOVE_END   = 44
BATCH_SIZE           = 200000   # 每批 reads 数（调大减少进程通信开销）


# ============================================================================
# pigz IO 封装
# ============================================================================
def open_fastq_read(path):
    """用 pigz 多线程解压，比 Python gzip 快 3-5x"""
    if path.endswith(".gz"):
        proc = subprocess.Popen(
            [PIGZ_BIN, "-dc", "-p", "4", path],
            stdout=subprocess.PIPE, text=True, bufsize=1 << 20
        )
        return proc.stdout, proc
    return open(path, "r"), None


def open_fastq_write(path):
    """用 pigz 多线程压缩写入"""
    proc = subprocess.Popen(
        [PIGZ_BIN, "-c", "-p", "2"],
        stdin=subprocess.PIPE,
        stdout=open(path, "wb"),
        bufsize=1 << 20
    )
    return proc.stdin, proc


# ============================================================================
# 工具函数
# ============================================================================
def hamming(a, b):
    return sum(x != y for x, y in zip(a, b))


def find_fixed_seq(sequence):
    # 先用 str.find 精确匹配（C 实现，比循环快很多）
    pos = sequence.find(FIXED_SEQUENCE, SEARCH_START, SEARCH_END + 15)
    if pos != -1:
        return pos
    # 精确匹配失败再做容错
    end = min(SEARCH_END, len(sequence) - 20)
    for i in range(SEARCH_START, end):
        if hamming(sequence[i:i+15], FIXED_SEQUENCE) <= FIXED_SEQ_MAX_ERR:
            return i
    return None


def extract_index(sequence):
    pos = find_fixed_seq(sequence)
    if pos is None or pos + 20 > len(sequence):
        return None
    return sequence[pos + 15: pos + 20]


def match_index(query, target_set):
    if query in target_set:
        return query
    if INDEX_MAX_ERR > 0:
        for t in target_set:
            if hamming(query, t) <= INDEX_MAX_ERR:
                return t
    return None


def trim_r2(seq, qual):
    if len(seq) < ADAPTER_REMOVE_END:
        return seq, qual
    return (seq[:ADAPTER_REMOVE_START] + seq[ADAPTER_REMOVE_END:],
            qual[:ADAPTER_REMOVE_START] + qual[ADAPTER_REMOVE_END:])


# ============================================================================
# chunk 处理（子进程）
# ============================================================================
def process_chunk(args):
    r1_lines, r2_lines, target_set = args
    buffers = defaultdict(lambda: {"r1": [], "r2": []})
    unmatched = {"r1": [], "r2": []}
    stats = {"total": 0, "adapter_found": 0, "matched": 0,
             "index_counts": defaultdict(int), "matched_counts": defaultdict(int)}

    for i in range(0, len(r1_lines), 4):
        if i + 3 >= len(r1_lines):
            break
        stats["total"] += 1
        r2_seq  = r2_lines[i + 1].rstrip("\n")
        r2_qual = r2_lines[i + 3].rstrip("\n")

        idx = extract_index(r2_seq)
        if idx is None:
            unmatched["r1"] += r1_lines[i:i+4]
            unmatched["r2"] += r2_lines[i:i+4]
            continue
        stats["adapter_found"] += 1
        stats["index_counts"][idx] += 1

        best = match_index(idx, target_set)
        if best is None:
            unmatched["r1"] += r1_lines[i:i+4]
            unmatched["r2"] += r2_lines[i:i+4]
            continue
        stats["matched"] += 1
        stats["matched_counts"][best] += 1

        r2_seq_t, r2_qual_t = trim_r2(r2_seq, r2_qual)
        buffers[best]["r1"] += r1_lines[i:i+4]
        buffers[best]["r2"] += [r2_lines[i], r2_seq_t + "\n", r2_lines[i+2], r2_qual_t + "\n"]

    return dict(buffers), unmatched, stats


# ============================================================================
# 主拆分函数
# ============================================================================
def demultiplex(r1_file, r2_file, output_dir, prefix, workers):
    target_set = set(INDEX_TO_MARK.keys())
    os.makedirs(output_dir, exist_ok=True)
    t0 = time.time()

    # 预开输出文件（pigz 压缩写入）
    out_files, out_procs = {}, {}
    for idx, mark in INDEX_TO_MARK.items():
        r1_fh, r1_proc = open_fastq_write(os.path.join(output_dir, f"{prefix}_{mark}_R1.fq.gz"))
        r2_fh, r2_proc = open_fastq_write(os.path.join(output_dir, f"{prefix}_{mark}_R2.fq.gz"))
        out_files[idx] = {"r1": r1_fh, "r2": r2_fh}
        out_procs[idx] = {"r1": r1_proc, "r2": r2_proc}
    # 未匹配 reads 输出
    unmatched_r1_fh, unmatched_r1_proc = open_fastq_write(os.path.join(output_dir, f"{prefix}_unmatched_R1.fq.gz"))
    unmatched_r2_fh, unmatched_r2_proc = open_fastq_write(os.path.join(output_dir, f"{prefix}_unmatched_R2.fq.gz"))

    total_stats = {"total": 0, "adapter_found": 0, "matched": 0,
                   "index_counts": defaultdict(int), "matched_counts": defaultdict(int)}

    def flush_results(results):
        for buffers, unmatched, stats in results:
            for idx, data in buffers.items():
                if data["r1"]:
                    out_files[idx]["r1"].write("".join(data["r1"]).encode())
                    out_files[idx]["r2"].write("".join(data["r2"]).encode())
            if unmatched["r1"]:
                unmatched_r1_fh.write("".join(unmatched["r1"]).encode())
                unmatched_r2_fh.write("".join(unmatched["r2"]).encode())
            total_stats["total"]         += stats["total"]
            total_stats["adapter_found"] += stats["adapter_found"]
            total_stats["matched"]       += stats["matched"]
            for k, v in stats["index_counts"].items():
                total_stats["index_counts"][k] += v
            for k, v in stats["matched_counts"].items():
                total_stats["matched_counts"][k] += v

    print(f"开始拆分，{workers} 个 CPU 进程 + pigz 多线程 IO")
    print(f"输出目录: {output_dir}  前缀: {prefix}")

    # pigz 多线程解压读取
    r1_in, r1_proc = open_fastq_read(r1_file)
    r2_in, r2_proc = open_fastq_read(r2_file)

    pool = mp.Pool(processes=workers)
    try:
        r1_buf, r2_buf, chunk_jobs = [], [], []
        for r1_line, r2_line in zip(r1_in, r2_in):
            r1_buf.append(r1_line)
            r2_buf.append(r2_line)
            if len(r1_buf) == BATCH_SIZE * 4:
                chunk_jobs.append(pool.apply_async(process_chunk, ((r1_buf, r2_buf, target_set),)))
                r1_buf, r2_buf = [], []
                if len(chunk_jobs) >= workers * 2:
                    flush_results([j.get() for j in chunk_jobs])
                    chunk_jobs = []
                    print(f"  已处理 {total_stats['total']:,} reads...")
        if r1_buf:
            chunk_jobs.append(pool.apply_async(process_chunk, ((r1_buf, r2_buf, target_set),)))
        if chunk_jobs:
            flush_results([j.get() for j in chunk_jobs])
    finally:
        pool.close()
        pool.join()
        r1_in.close(); r2_in.close()
        if r1_proc: r1_proc.wait()
        if r2_proc: r2_proc.wait()
        for idx in out_files:
            out_files[idx]["r1"].close()
            out_files[idx]["r2"].close()
        for idx in out_procs:
            out_procs[idx]["r1"].wait()
            out_procs[idx]["r2"].wait()
        unmatched_r1_fh.close(); unmatched_r2_fh.close()
        unmatched_r1_proc.wait(); unmatched_r2_proc.wait()

    # ============================================================================
    # 统计报告
    # ============================================================================
    print("\n" + "=" * 60)
    print("拆分完成！")
    print(f"总 reads:     {total_stats['total']:,}")
    print(f"固定序列检出: {total_stats['adapter_found']:,} ({total_stats['adapter_found']/max(total_stats['total'],1)*100:.1f}%)")
    print(f"成功匹配:     {total_stats['matched']:,} ({total_stats['matched']/max(total_stats['total'],1)*100:.1f}%)")
    print("\n各抗体 reads：")
    for idx, mark in INDEX_TO_MARK.items():
        count = total_stats["matched_counts"].get(idx, 0)
        print(f"  {idx} ({mark}): {count:,}  →  {prefix}_{mark}_R1.fq.gz")

    unmatched = total_stats["total"] - total_stats["matched"]
    report_path = os.path.join(output_dir, "demultiplex_report.txt")
    with open(report_path, "w") as f:
        f.write("=" * 60 + "\n")
        f.write("CUT&TAG 拆分统计报告\n")
        f.write("=" * 60 + "\n\n")
        f.write("[总体统计]\n")
        f.write(f"总 reads 数\t{total_stats['total']:,}\n")
        f.write(f"找到固定序列\t{total_stats['adapter_found']:,}\t({total_stats['adapter_found']/max(total_stats['total'],1)*100:.1f}%)\n")
        f.write(f"成功匹配 index\t{total_stats['matched']:,}\t({total_stats['matched']/max(total_stats['total'],1)*100:.1f}%)\n")
        f.write(f"未匹配 reads\t{unmatched:,}\t({unmatched/max(total_stats['total'],1)*100:.1f}%)\n\n")
        f.write("[各抗体/模态统计]\n")
        f.write(f"{'Index':<8}{'Mark':<12}{'Reads':>12}{'占总数%':>10}{'占匹配%':>10}\t输出文件\n")
        f.write("-" * 80 + "\n")
        for idx, mark in INDEX_TO_MARK.items():
            count = total_stats["matched_counts"].get(idx, 0)
            pct_t = count / max(total_stats["total"],   1) * 100
            pct_m = count / max(total_stats["matched"], 1) * 100
            f.write(f"{idx:<8}{mark:<12}{count:>12,}{pct_t:>9.1f}%{pct_m:>9.1f}%\t{prefix}_{mark}_R1.fq.gz\n")
        f.write("\n[发现的所有 index 分布（前30）]\n")
        f.write(f"{'Index':<8}{'Reads':>12}{'匹配到':>14}\n")
        f.write("-" * 40 + "\n")
        sorted_idx = sorted(total_stats["index_counts"].items(), key=lambda x: x[1], reverse=True)
        for raw_idx, count in sorted_idx[:30]:
            best = match_index(raw_idx, set(INDEX_TO_MARK.keys()))
            matched_mark = f"{raw_idx}→{INDEX_TO_MARK[best]}" if best else "未匹配"
            f.write(f"{raw_idx:<8}{count:>12,}{matched_mark:>14}\n")
        if len(sorted_idx) > 30:
            f.write(f"... 还有 {len(sorted_idx) - 30} 个其他 index\n")

        elapsed = time.time() - t0
        f.write(f"\n[运行时间]\n{elapsed:.1f} 秒 ({elapsed/60:.1f} 分钟)\n")

    print(f"\n统计报告: {report_path}")
    print("=" * 60)


# ============================================================================
# Barcode 提取（fastp + extract_bc + 合并）
# ============================================================================
def extract_barcodes(sample_name, r1_file, r2_file, output_dir):
    """
    对单个样本执行 fastp 质控 + barcode 提取 + 合并
    输出到 output_dir/{sample_name}_extracted/
    """
    extract_dir = os.path.join(output_dir, f"{sample_name}_extracted")
    clean_dir   = os.path.join(output_dir, f"{sample_name}_clean")
    os.makedirs(extract_dir, exist_ok=True)
    os.makedirs(clean_dir,   exist_ok=True)

    r1_out = os.path.join(extract_dir, f"{sample_name}_S1_L001_R1_001.fastq.gz")
    r2_out = os.path.join(extract_dir, f"{sample_name}_S1_L001_R2_001.fastq.gz")
    r3_out = os.path.join(extract_dir, f"{sample_name}_S1_L001_R3_001.fastq.gz")
    i1_out = os.path.join(extract_dir, f"{sample_name}_S1_L001_I1_001.fastq.gz")

    # 已全部完成则跳过
    if all(os.path.exists(f) for f in [r1_out, r2_out, r3_out, i1_out]):
        print(f"  [{sample_name}] barcode 已提取，跳过")
        return extract_dir

    print(f"  [{sample_name}] 步骤A: fastp 质控 + 分片...")
    fastp_r1 = os.path.join(clean_dir, f"{sample_name}_R1.clean.fq")
    fastp_r2 = os.path.join(clean_dir, f"{sample_name}_R2.clean.fq")
    if not os.path.exists(f"{clean_dir}/0001.{sample_name}_R1.clean.fq"):
        cmd = [
            FASTP_BIN,
            "-i", r1_file, "-I", r2_file,
            "-o", fastp_r1, "-O", fastp_r2,
            "--json", os.path.join(clean_dir, f"{sample_name}.json"),
            "--html", os.path.join(clean_dir, f"{sample_name}.html"),
            "-w", "12", "--length_required", "34",
            "-s", str(SPLIT_NUM),
        ]
        with open(os.path.join(output_dir, f"fastp_{sample_name}.log"), "w") as log:
            subprocess.run(cmd, stdout=log, stderr=log, check=True)

    print(f"  [{sample_name}] 步骤B: 提取 barcode（并行 {SPLIT_NUM} 片）...")
    first_chunk = os.path.join(clean_dir, f"0001.{sample_name}", "R1_filter.fq")
    if not os.path.exists(first_chunk):
        # 生成分片编号列表并并行处理
        chunks = [f"{i:04d}" for i in range(1, SPLIT_NUM + 1)]
        procs = []
        for chunk in chunks:
            r1_chunk = os.path.join(clean_dir, f"{chunk}.{sample_name}_R1.clean.fq")
            r2_chunk = os.path.join(clean_dir, f"{chunk}.{sample_name}_R2.clean.fq")
            out_prefix = os.path.join(clean_dir, f"{chunk}.{sample_name}")
            cmd = [EXTRACT_BC_BIN,
                   "--r1_fastq_path", r1_chunk,
                   "--r2_fastq_path", r2_chunk,
                   "-o", out_prefix]
            procs.append(subprocess.Popen(cmd))
        for p in procs:
            p.wait()

    print(f"  [{sample_name}] 步骤C: 合并文件...")
    import glob as _glob

    def merge_files(pattern, out_path):
        if not os.path.exists(out_path):
            files = sorted(_glob.glob(pattern))
            pigz = subprocess.Popen(
                [PIGZ_BIN, "-c", "-p", "2"],
                stdin=subprocess.PIPE,
                stdout=open(out_path, "wb"),
            )
            for f in files:
                with open(f, "rb") as fh:
                    pigz.stdin.write(fh.read())
            pigz.stdin.close()
            pigz.wait()

    merge_files(os.path.join(clean_dir, f"*{sample_name}", "R1_filter.fq"), r1_out)
    merge_files(os.path.join(clean_dir, f"*{sample_name}", "barcode.fq"),   r2_out)
    merge_files(os.path.join(clean_dir, f"*{sample_name}", "R2_filter.fq"), r3_out)
    merge_files(os.path.join(clean_dir, f"*{sample_name}", "I1.fq"),        i1_out)

    # 清理临时文件
    import shutil
    shutil.rmtree(clean_dir, ignore_errors=True)

    print(f"  [{sample_name}] barcode 提取完成 → {extract_dir}")
    return extract_dir


# ============================================================================
# 入口
# ============================================================================
def main():
    import glob

    r1_files = sorted(glob.glob("*_R1*.fq.gz") + glob.glob("*_R1*.fastq.gz"))
    r2_files = sorted(glob.glob("*_R2*.fq.gz") + glob.glob("*_R2*.fastq.gz"))

    if not r1_files or not r2_files:
        print("错误：当前目录下未找到 *_R1*.fq.gz 或 *_R2*.fq.gz 文件")
        return

    r1_file = r1_files[0]
    r2_file = r2_files[0]
    prefix     = os.path.basename(os.path.abspath("."))
    output_dir = os.path.join(".", "demultiplex_out")
    workers    = min(len(INDEX_TO_MARK), mp.cpu_count())

    print(f"R1: {r1_file}")
    print(f"R2: {r2_file}")
    print(f"样本前缀: {prefix}")

    # 第一步：demultiplex
    t0 = time.time()
    demultiplex(r1_file, r2_file, output_dir, prefix, workers)
    print(f"\n总运行时间（demultiplex）: {time.time() - t0:.1f} 秒")

    # 第二步：对每个 mark 提取 barcode
    print("\n" + "=" * 60)
    print("开始 barcode 提取（fastp + extract_bc + 合并）")
    print("=" * 60)

    t1 = time.time()
    for idx, mark in INDEX_TO_MARK.items():
        sample_name = f"{prefix}_{mark}"
        r1 = os.path.join(output_dir, f"{sample_name}_R1.fq.gz")
        r2 = os.path.join(output_dir, f"{sample_name}_R2.fq.gz")
        if not os.path.exists(r1):
            print(f"  [{sample_name}] 未找到 demultiplex 输出，跳过")
            continue
        extract_barcodes(sample_name, r1, r2, output_dir)

    print(f"\nBarcode 提取总耗时: {time.time() - t1:.1f} 秒")
    print("\n" + "=" * 60)
    print("全部完成！各 mark 的提取文件位于：")
    for idx, mark in INDEX_TO_MARK.items():
        extract_dir = os.path.join(output_dir, f"{prefix}_{mark}_extracted")
        if os.path.isdir(extract_dir):
            print(f"  {mark}: {extract_dir}")
    print("=" * 60)
    print("\n现在可以在 demultiplex_out/ 目录下运行脚本2进行比对。")


if __name__ == "__main__":
    main()
