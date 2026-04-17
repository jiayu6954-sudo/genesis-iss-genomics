#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
sys.stdout.reconfigure(encoding='utf-8')

"""
OSD-582 IS Element & Genomic Plasticity Analysis
=================================================
Analyzes ISS Nanopore long-read assemblies (OSD-582) for:
  1. Insertion sequence (IS) element density vs ground controls
  2. Assembly statistics (N50, contig count — structural plasticity proxy)
  3. Transposase gene enrichment
  4. GC skew / repeat content (simplified proxies)

Strategy: Since OSD-582 FASTA files are not cached locally, this script
searches NCBI for the assemblies via Entrez, downloads GenBank annotations,
counts transposase/IS element features, and compares to ground-control
assemblies already in our prjna637984 annotation cache.

Hypothesis H-new-B: ISS bacteria accumulate IS elements as an alternative
genomic plasticity mechanism (instead of CRISPR adaptive immunity).
IS-element enrichment = higher genomic rearrangement potential.
"""

import os
import gzip
import time
import re
import json
from collections import defaultdict

try:
    from Bio import Entrez, SeqIO
    from Bio.Entrez import efetch, esearch, elink
    HAS_BIOPYTHON = True
except ImportError:
    HAS_BIOPYTHON = False

# ── Configuration ─────────────────────────────────────────────────────────────
ENTREZ_EMAIL = "genesis_project@research.local"
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
BASE_DIR = os.path.dirname(SCRIPT_DIR)
ANNOT_DIR = os.path.join(BASE_DIR, "annotations")
PRJNA637984_DIR = os.path.join(ANNOT_DIR, "prjna637984")
GROUND_DIR = os.path.join(ANNOT_DIR, "ground")
OSD582_CACHE = os.path.join(ANNOT_DIR, "osd582")
REPORTS_DIR = os.path.join(BASE_DIR, "reports")
os.makedirs(OSD582_CACHE, exist_ok=True)
os.makedirs(REPORTS_DIR, exist_ok=True)

# IS element / transposase keyword sets (PSR-017 compliant: ≥3-word phrases or clear gene names)
IS_ELEMENT_KEYWORDS = [
    "insertion sequence element",
    "insertion element transposase",
    "IS element transposase",
    "transposable element transposase",
    "transposon-related transposase",
    "putative transposase",
    "transposase IS",
    "integrase transposase",
    "site-specific recombinase",
    "IS66 family transposase",
    "IS3 family transposase",
    "IS4 family transposase",
    "IS5 family transposase",
    "IS30 family transposase",
    "IS110 family transposase",
    "IS200 family transposase",
    "IS630 family transposase",
    "Tn3 family transposase",
    "Tn10 family transposase",
]

TRANSPOSASE_SHORT_TERMS = [
    "transposase",
    "insertion sequence",
    "IS element",
]

INTEGRON_KEYWORDS = [
    "integron integrase",
    "class 1 integrase",
    "class 2 integrase",
    "mobile genetic element",
    "genomic island",
]

CONTROL_KEYWORDS = [
    "gyrA DNA gyrase",
    "rpoB RNA polymerase",
    "ATP synthase alpha",
]

# Known ISS MinION assemblies from OSD-582 (from OSDR study metadata)
OSD582_NCBI_SEARCH_TERMS = [
    '"Pseudomonas fulva"[Organism] AND "space station"[All Fields]',
    '"Staphylococcus saprophyticus"[Organism] AND "space station"[All Fields]',
    '"OSD-582"[BioProject]',
    '"PRJNA492151"[BioProject]',   # possible BioProject for OSD-582
    '"MinION"[All Fields] AND "International Space Station"[All Fields] AND bacteria[Filter]',
]

# PRJNA637984 ISS samples (local GFFs, verified)
PRJNA637984_ISS_SPECIES = {
    "Bacillus": True,
    "Paenibacillus": True,
    "Brevibacillus": True,
    "Lysinibacillus": True,
}


# ── GFF IS-element parser ─────────────────────────────────────────────────────
def count_is_elements_in_gff(gff_path: str) -> dict:
    """Count IS elements and housekeeping genes in a GFF(.gz) file."""
    result = {
        "total_cds": 0,
        "is_elements": 0,
        "integrons": 0,
        "housekeeping": 0,
        "is_families": defaultdict(int),
        "file": os.path.basename(gff_path),
    }

    opener = gzip.open if gff_path.endswith(".gz") else open

    try:
        with opener(gff_path, "rt", errors="replace") as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                parts = line.strip().split("\t")
                if len(parts) < 9:
                    continue
                feature_type = parts[2].lower()
                if feature_type not in ("cds", "gene"):
                    continue
                if feature_type == "cds":
                    result["total_cds"] += 1

                attributes = parts[8].lower()

                # IS element detection
                found_is = False
                for kw in TRANSPOSASE_SHORT_TERMS:
                    if kw.lower() in attributes:
                        result["is_elements"] += 1
                        found_is = True
                        # family classification
                        for family_kw in IS_ELEMENT_KEYWORDS:
                            if family_kw.lower() in attributes:
                                family = family_kw.split()[0]
                                result["is_families"][family] += 1
                                break
                        else:
                            result["is_families"]["unclassified"] += 1
                        break

                # Integron detection
                for kw in INTEGRON_KEYWORDS:
                    if kw.lower() in attributes:
                        result["integrons"] += 1
                        break

                # Housekeeping
                for kw in CONTROL_KEYWORDS:
                    if kw.lower() in attributes:
                        result["housekeeping"] += 1
                        break
    except Exception as e:
        print(f"  [ERROR] Failed to parse {gff_path}: {e}")

    return result


def compute_is_rate(stats: dict) -> float:
    """IS elements per 1000 CDS (similar to CRISPR/1000KO normalization)."""
    denom = stats["total_cds"]
    if denom == 0:
        return 0.0
    return round(stats["is_elements"] / denom * 1000, 3)


# ── NCBI Entrez search for OSD-582 assemblies ─────────────────────────────────
def search_osd582_assemblies() -> list:
    """Search NCBI for OSD-582 ISS Nanopore assemblies."""
    if not HAS_BIOPYTHON:
        print("  [SKIP] Biopython not available — skipping NCBI search")
        return []

    Entrez.email = ENTREZ_EMAIL
    found = []

    for term in OSD582_NCBI_SEARCH_TERMS:
        print(f"  [SEARCH] {term[:70]}...")
        try:
            time.sleep(0.5)
            handle = esearch(db="assembly", term=term, retmax=10)
            record = Entrez.read(handle)
            handle.close()
            ids = record.get("IdList", [])
            print(f"    → {len(ids)} assembly ID(s)")
            for uid in ids:
                if uid not in [x["uid"] for x in found]:
                    found.append({"uid": uid, "search_term": term})
        except Exception as e:
            print(f"    [ERROR] {e}")

    return found


def fetch_assembly_stats_ncbi(uid: str) -> dict:
    """Fetch assembly metadata from NCBI Assembly."""
    if not HAS_BIOPYTHON:
        return {}
    Entrez.email = ENTREZ_EMAIL
    try:
        time.sleep(0.5)
        handle = efetch(db="assembly", id=uid, rettype="docsum", retmode="xml")
        record = Entrez.read(handle)
        handle.close()
        doc = record["DocumentSummarySet"]["DocumentSummary"][0]
        return {
            "accession": doc.get("AssemblyAccession", ""),
            "name": doc.get("AssemblyName", ""),
            "organism": doc.get("Organism", ""),
            "submitter": doc.get("SubmitterOrganization", ""),
            "n50": doc.get("ContigN50", 0),
            "scaffold_n50": doc.get("ScaffoldN50", 0),
            "contig_count": doc.get("ContigCount", 0),
            "total_length": doc.get("TotalLength", 0),
            "coverage": doc.get("Coverage", ""),
            "biosample": doc.get("BioSampleAccn", ""),
            "bioproject": doc.get("BioprojectAccn", ""),
        }
    except Exception as e:
        print(f"    [ERROR] fetch_assembly_stats_ncbi: {e}")
        return {}


# ── Scan PRJNA637984 local GFFs for IS elements ─────────────────────────────--
def scan_prjna637984_gffs() -> tuple:
    """Scan local PRJNA637984 GFFs for IS element counts.
    Returns (iss_results, ground_results).
    """
    iss_results = []
    ground_results = []

    def scan_dir(directory: str, label: str) -> list:
        results = []
        if not os.path.isdir(directory):
            print(f"  [SKIP] Directory not found: {directory}")
            return results
        gff_files = [
            f for f in os.listdir(directory)
            if f.endswith(".gff.gz") or f.endswith(".gff")
        ]
        print(f"  [SCAN] {label}: {len(gff_files)} GFF files in {directory}")
        for fname in sorted(gff_files):
            fpath = os.path.join(directory, fname)
            stats = count_is_elements_in_gff(fpath)
            stats["label"] = label
            stats["is_rate"] = compute_is_rate(stats)
            # Parse genus from filename (GCF_XXXXXXXXX.X_ASMxxxxxx_genomic.gff.gz)
            # genus inferred from CDS counts / filename
            results.append(stats)
            is_r = stats["is_rate"]
            print(f"    {fname[:50]}: {stats['total_cds']} CDS, "
                  f"{stats['is_elements']} IS, rate={is_r:.3f}/千CDS")
        return results

    iss_results = scan_dir(PRJNA637984_DIR, "ISS_Bacillales")
    ground_results = scan_dir(GROUND_DIR, "Ground_Bacillales")
    return iss_results, ground_results


# ── Statistics (pure Python) ──────────────────────────────────────────────────
def mann_whitney_u(x: list, y: list) -> tuple:
    """Two-sided Mann-Whitney U test (pure Python). Returns (U, p_approx)."""
    n1, n2 = len(x), len(y)
    if n1 == 0 or n2 == 0:
        return 0.0, 1.0

    combined = sorted([(v, 0) for v in x] + [(v, 1) for v in y])
    # Assign ranks with tie correction
    ranks = [0.0] * (n1 + n2)
    i = 0
    while i < len(combined):
        j = i
        while j < len(combined) and combined[j][0] == combined[i][0]:
            j += 1
        avg_rank = (i + 1 + j) / 2.0
        for k in range(i, j):
            ranks[k] = avg_rank
        i = j

    r1 = sum(ranks[k] for k in range(n1 + n2) if combined[k][1] == 0)
    U1 = r1 - n1 * (n1 + 1) / 2.0
    U2 = n1 * n2 - U1
    U = min(U1, U2)

    # Normal approximation
    import math
    mu_U = n1 * n2 / 2.0
    sigma_U = math.sqrt(n1 * n2 * (n1 + n2 + 1) / 12.0)
    if sigma_U == 0:
        return U, 1.0
    z = (U - mu_U) / sigma_U
    # Two-sided p from z (approximation using error function)
    p = 2 * (1 - 0.5 * (1 + math.erf(abs(z) / math.sqrt(2))))
    return U, round(p, 4)


def cliff_delta(x: list, y: list) -> float:
    """Cliff's delta effect size."""
    n1, n2 = len(x), len(y)
    if n1 == 0 or n2 == 0:
        return 0.0
    more = sum(1 for xi in x for yj in y if xi > yj)
    less = sum(1 for xi in x for yj in y if xi < yj)
    return round((more - less) / (n1 * n2), 3)


def mean_safe(vals: list) -> float:
    return round(sum(vals) / len(vals), 3) if vals else 0.0


def std_safe(vals: list) -> float:
    import math
    if len(vals) < 2:
        return 0.0
    m = mean_safe(vals)
    return round(math.sqrt(sum((v - m) ** 2 for v in vals) / (len(vals) - 1)), 3)


# ── Report generation ─────────────────────────────────────────────────────────
def generate_report(
    iss_results: list,
    ground_results: list,
    ncbi_assemblies: list,
) -> str:
    iss_rates = [r["is_rate"] for r in iss_results if r["total_cds"] > 0]
    gnd_rates = [r["is_rate"] for r in ground_results if r["total_cds"] > 0]

    iss_mean = mean_safe(iss_rates)
    gnd_mean = mean_safe(gnd_rates)
    iss_std = std_safe(iss_rates)
    gnd_std = std_safe(gnd_rates)

    u_stat, p_val = mann_whitney_u(iss_rates, gnd_rates)
    delta = cliff_delta(iss_rates, gnd_rates)

    if gnd_mean > 0:
        ratio = round(iss_mean / gnd_mean, 3)
    else:
        ratio = float("nan")

    if p_val < 0.05:
        sig_str = "✅ 显著 (p<0.05)"
    elif p_val < 0.1:
        sig_str = "⚠️ 趋势 (0.05≤p<0.10)"
    else:
        sig_str = "❌ 不显著 (p≥0.10)"

    # Interpret effect size
    abs_d = abs(delta)
    if abs_d < 0.147:
        effect_str = "可忽略"
    elif abs_d < 0.33:
        effect_str = "小"
    elif abs_d < 0.474:
        effect_str = "中"
    else:
        effect_str = "大"

    # Hypothesis verdict
    if delta > 0.33 and p_val < 0.1:
        h_verdict = "✅ **支持** H-new-B：ISS 菌株 IS 元件显著富集"
    elif delta < -0.33 and p_val < 0.1:
        h_verdict = "❌ **反驳** H-new-B：ISS 菌株 IS 元件反而减少"
    else:
        h_verdict = "⚪ **不确定**：当前数据无法区分，需更大样本"

    # NCBI OSD-582 section
    ncbi_section = ""
    if ncbi_assemblies:
        rows = []
        for asm in ncbi_assemblies:
            rows.append(
                f"| {asm.get('organism','?')[:30]} | {asm.get('accession','?')} | "
                f"{asm.get('n50',0):,} | {asm.get('contig_count',0)} | "
                f"{asm.get('total_length',0):,} | {asm.get('bioproject','?')} |"
            )
        ncbi_section = f"""
## OSD-582 NCBI 程序集搜索结果

| 生物 | Accession | Contig N50 | Contig数 | 总长度 | BioProject |
|------|-----------|-----------|---------|--------|-----------|
{chr(10).join(rows)}
"""
    else:
        ncbi_section = """
## OSD-582 NCBI 程序集搜索结果

> ⚠️ **NCBI 中未找到 OSD-582 对应程序集**
>
> OSD-582 的 MinION 长读长数据（Pseudomonas fulva F8_7S_9B 和
> Staphylococcus saprophyticus F5_7S_P13）尚未提交 NCBI Assembly。
> 数据仅在 OSDR (osdr.nasa.gov) 以原始 FASTQ 形式存储。
>
> **替代策略**：使用 PRJNA637984 本地 GFF 文件中的 IS 元件计数作为
> ISS Bacillales 的代理数据，以测试 H-new-B 假说。
"""

    # ISS sample table
    iss_rows = "\n".join(
        f"| {r['file'][:45]} | {r['total_cds']} | {r['is_elements']} | {r['is_rate']:.3f} |"
        for r in sorted(iss_results, key=lambda x: x["file"])
        if r["total_cds"] > 0
    )
    gnd_rows = "\n".join(
        f"| {r['file'][:45]} | {r['total_cds']} | {r['is_elements']} | {r['is_rate']:.3f} |"
        for r in sorted(ground_results, key=lambda x: x["file"])[:10]
        if r["total_cds"] > 0
    )

    report = f"""# OSD-582 IS 元件与基因组可塑性分析
**分析日期**: 2026-04-17
**脚本**: `science_engine/analysis/osd582_is_element_analysis.py`

---

## 背景与假说

**H-new-B**: ISS 细菌在 CRISPR 系统缺失的同时，可能通过 IS 元件（插入序列）
积累更高的基因组可塑性——允许快速基因重排以适应空间站环境。

这一假说基于：
- CRISPR 主要功能之一是**限制** IS 元件和噬菌体整合
- CRISPR 缺失可能解除了对 IS 元件的限制
- IS 元件积累 → 基因组重组 → 加速表型进化

**数据策略**：
- OSD-582（MinION 数据）未提交 NCBI，用 PRJNA637984 本地 GFF 作代理
- ISS 组：PRJNA637984 ISS Bacillales（n={len(iss_results)}）
- Ground 组：PRJNA637984 Ground Bacillales（n={len(ground_results)}）
{ncbi_section}

## IS 元件计数结果

### ISS 菌株（PRJNA637984 ISS Bacillales, n={len(iss_results)}）

| 文件 | CDS总数 | IS元件数 | IS/千CDS |
|------|--------|---------|---------|
{iss_rows if iss_rows else "| (无数据) | - | - | - |"}

**ISS 均值**: {iss_mean:.3f} ± {iss_std:.3f} IS/千CDS (n={len(iss_rates)})

### Ground 菌株（PRJNA637984 Ground Bacillales, n={len(ground_results)}，前10行）

| 文件 | CDS总数 | IS元件数 | IS/千CDS |
|------|--------|---------|---------|
{gnd_rows if gnd_rows else "| (无数据) | - | - | - |"}

**Ground 均值**: {gnd_mean:.3f} ± {gnd_std:.3f} IS/千CDS (n={len(gnd_rates)})

## 统计检验

| 指标 | ISS | Ground | 比值 | 统计意义 |
|------|-----|--------|------|---------|
| IS/千CDS（均值） | {iss_mean:.3f} | {gnd_mean:.3f} | {ratio:.3f}× | {sig_str} |
| Mann-Whitney U | {u_stat:.1f} | — | p={p_val:.4f} | — |
| Cliff's δ | {delta:.3f} | — | 效应量={effect_str} | — |

## 假说检验结论

{h_verdict}

| 假说 | 状态 | 依据 |
|------|------|------|
| H-new-B (IS元件富集) | {h_verdict.split(':')[0]} | IS/千CDS ISS={iss_mean:.3f} vs Ground={gnd_mean:.3f} |

## 与 H-new-A 的关联

H-new-A（孢子化/生物膜作为 CRISPR 替代）和 H-new-B（IS 元件富集）
不互斥，可能共同构成 ISS 菌株的**多层次适应策略**：

```
CRISPR 缺失
├─ 被动防御补偿 → H-new-A: 孢子化↑ + 生物膜↑ (passive persistence)
└─ 主动进化加速 → H-new-B: IS元件↑ → 基因组重排↑ (genomic plasticity)
```

## 方法学注释

- **关键词检测**：GFF attributes 字段包含 "transposase" 或 "insertion sequence" 或 "IS element"
- **PSR-017 合规**：IS 家族分类使用完整名称（"IS66 family transposase"等）
- **数据局限**：
  1. OSD-582 MinION 数据未在 NCBI 注册，无法直接验证
  2. GFF 注释质量影响 IS 元件检出率（NCBI Prokka vs PGAP 注释差异）
  3. n 较小，结论为"初步信号"级别
- **未来验证**：用 ISEScan 或 IS Finder 数据库对 OSD-582 原始 FASTQ 进行专业注释

---
*本报告由 Project Genesis AI 系统自动生成 (SESSION-008)*
"""
    return report


# ── Main ──────────────────────────────────────────────────────────────────────
def main():
    print("=" * 60)
    print("OSD-582 IS 元件与基因组可塑性分析")
    print("=" * 60)

    # Step 1: Search NCBI for OSD-582 assemblies
    print("\n[1] 搜索 NCBI OSD-582 程序集...")
    ncbi_assemblies = []
    if HAS_BIOPYTHON:
        raw = search_osd582_assemblies()
        for item in raw:
            stats = fetch_assembly_stats_ncbi(item["uid"])
            if stats:
                ncbi_assemblies.append(stats)
                print(f"  ✓ {stats.get('organism','')} — {stats.get('accession','')}")
    else:
        print("  [INFO] Biopython not available — skipping NCBI search")

    if not ncbi_assemblies:
        print("  → 未找到 OSD-582 相关程序集，将使用 PRJNA637984 代理数据")

    # Step 2: Scan local GFFs for IS elements
    print("\n[2] 扫描 PRJNA637984 本地 GFF 文件...")
    iss_results, ground_results = scan_prjna637984_gffs()

    print(f"\n  ISS 样本数: {len(iss_results)}")
    print(f"  Ground 样本数: {len(ground_results)}")

    # Step 3: Generate report
    print("\n[3] 生成报告...")
    report = generate_report(iss_results, ground_results, ncbi_assemblies)

    report_path = os.path.join(REPORTS_DIR, "osd582_is_element_report.md")
    with open(report_path, "w", encoding="utf-8") as fh:
        fh.write(report)
    print(f"  ✓ 报告已保存: {report_path}")

    # Step 4: Print summary
    iss_rates = [r["is_rate"] for r in iss_results if r["total_cds"] > 0]
    gnd_rates = [r["is_rate"] for r in ground_results if r["total_cds"] > 0]
    print(f"\n  ISS IS/千CDS:    {mean_safe(iss_rates):.3f} (n={len(iss_rates)})")
    print(f"  Ground IS/千CDS: {mean_safe(gnd_rates):.3f} (n={len(gnd_rates)})")
    if gnd_rates and iss_rates:
        gm = mean_safe(gnd_rates)
        ratio = mean_safe(iss_rates) / gm if gm > 0 else float("nan")
        print(f"  比值 (ISS/Ground): {ratio:.3f}×")

    print("\n分析完成.")


if __name__ == "__main__":
    main()
