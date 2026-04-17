#!/usr/bin/env e:/miniconda3/python.exe
# -*- coding: utf-8 -*-
"""
GLDS-224 宏基因组物种组成分析
Project Genesis · SESSION-008

目标：通过宏基因组 TSV 中的 domain/phylum/class/genus/species 列，
      比较 ISS vs Ground 群落物种组成差异

数据：GLDS-224 6个样本的 annotation TSV（已缓存）
      ISS:    SRR5220488 (HEPA no-PMA), SRR5220487 (HEPA PMA)
              SRR5197511 (dust no-PMA), SRR5197512 (dust PMA)
      Ground: SRR5197507 (dust no-PMA), SRR5220486 (dust PMA)

科学问题：
1. ISS 是否富集 Firmicutes（Bacillales）？
2. ISS 群落多样性是否低于 Ground？
3. ISS dust 与 Ground dust 的物种组成有何差异？
"""

import sys
import csv
import gzip
from pathlib import Path
from collections import Counter, defaultdict

sys.stdout.reconfigure(encoding='utf-8')

PROJECT_ROOT = Path('e:/miniconda3/envs/llama-env/genesis_project')
CACHE_DIR    = PROJECT_ROOT / 'science_engine' / 'annotations' / 'glds224'
OUT_DIR      = PROJECT_ROOT / 'science_engine' / 'preprocess' / 'output'
REPORT_DIR   = PROJECT_ROOT / 'science_engine' / 'reports'

SAMPLES = {
    'SRR5197507': {'label': 'SAF Debris no-PMA', 'env': 'Ground', 'type': 'dust', 'pma': False},
    'SRR5220486': {'label': 'SAF Debris PMA',    'env': 'Ground', 'type': 'dust', 'pma': True},
    'SRR5220488': {'label': 'ISS HEPA no-PMA',   'env': 'ISS',    'type': 'hepa', 'pma': False},
    'SRR5220487': {'label': 'ISS HEPA PMA',      'env': 'ISS',    'type': 'hepa', 'pma': True},
    'SRR5197511': {'label': 'ISS Debris no-PMA', 'env': 'ISS',    'type': 'dust', 'pma': False},
    'SRR5197512': {'label': 'ISS Debris PMA',    'env': 'ISS',    'type': 'dust', 'pma': True},
}


def parse_taxonomy_tsv(tsv_path: Path) -> dict:
    """
    解析宏基因组 TSV，统计各分类级别的基因计数。
    TSV 格式：gene_ID, coverage, KO_ID, KO_function, taxid,
               domain, phylum, class, order, family, genus, species
    """
    result = {
        'total': 0,
        'annotated': 0,        # 有非 'NA' 物种注释的基因数
        'domain':  Counter(),
        'phylum':  Counter(),
        'class':   Counter(),
        'order':   Counter(),
        'family':  Counter(),
        'genus':   Counter(),
        'species': Counter(),
    }
    NA_VALS = {'NA', "'NA'", 'N/A', '', 'unknown', 'unclassified', 'Unclassified'}

    with open(tsv_path, 'r', encoding='utf-8', errors='replace') as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        fields = {f.lower().strip().strip('[]'): f for f in reader.fieldnames}

        def col(key):
            for alias in [key, key.capitalize()]:
                if alias in fields:
                    return fields[alias]
            return None

        domain_col  = col('domain')
        phylum_col  = col('phylum')
        class_col   = col('class')
        order_col   = col('order')
        family_col  = col('family')
        genus_col   = col('genus')
        species_col = col('species')

        for row in reader:
            result['total'] += 1
            sp = row.get(species_col, '').strip().strip("'") if species_col else ''
            if sp and sp not in NA_VALS:
                result['annotated'] += 1

            def clean(col_name):
                if not col_name:
                    return 'unclassified'
                v = row.get(col_name, '').strip().strip("'")
                return v if (v and v not in NA_VALS) else 'unclassified'

            result['domain'][clean(domain_col)]   += 1
            result['phylum'][clean(phylum_col)]   += 1
            result['class'][clean(class_col)]     += 1
            result['order'][clean(order_col)]     += 1
            result['family'][clean(family_col)]   += 1
            result['genus'][clean(genus_col)]     += 1
            result['species'][clean(species_col)] += 1

    return result


def top_taxa(counter: Counter, n: int = 10, exclude_unclass: bool = False) -> list:
    """返回前 n 个分类群（含/排除 unclassified）"""
    if exclude_unclass:
        items = [(k, v) for k, v in counter.items() if k != 'unclassified']
    else:
        items = list(counter.items())
    items.sort(key=lambda x: -x[1])
    return items[:n]


def compute_diversity(counter: Counter, exclude_unclass: bool = True) -> dict:
    """计算 Shannon 多样性指数 H'"""
    import math
    items = [(k, v) for k, v in counter.items()
             if (not exclude_unclass or k != 'unclassified') and v > 0]
    total = sum(v for _, v in items)
    if total == 0:
        return {'richness': 0, 'shannon': 0.0, 'evenness': 0.0}
    H = -sum((v / total) * math.log(v / total) for _, v in items)
    richness = len(items)
    evenness = H / math.log(richness) if richness > 1 else 1.0
    return {'richness': richness, 'shannon': H, 'evenness': evenness}


def main():
    print('=' * 70)
    print('GLDS-224 宏基因组物种组成分析')
    print('Project Genesis · SESSION-008')
    print('=' * 70)

    sample_data = {}
    for srr, meta in SAMPLES.items():
        tsv_name = f'GLDS-224_GMetagenomics_{srr}-gene-coverage-annotation-and-tax.tsv'
        tsv_path = CACHE_DIR / tsv_name
        if not tsv_path.exists():
            print(f'  [缺失] {tsv_name}')
            continue
        print(f'\n处理 {srr} [{meta["label"]}]...')
        data = parse_taxonomy_tsv(tsv_path)
        data.update(meta)
        data['srr'] = srr
        data['tax_coverage'] = data['annotated'] / data['total'] * 100 if data['total'] else 0
        sample_data[srr] = data
        print(f'  总基因: {data["total"]:,}, 有物种注释: {data["annotated"]:,} '
              f'({data["tax_coverage"]:.1f}%)')
        print(f'  Top 5 phyla: {top_taxa(data["phylum"], 5, exclude_unclass=True)[:5]}')

    # ── 群落组成对比 ─────────────────────────────────────────────────
    print('\n' + '=' * 70)
    print('群落组成对比（所有基因，含未分类）')
    print('=' * 70)

    # 按 env + type 分组
    groups = {
        'Ground_dust': [v for k, v in sample_data.items() if v['env'] == 'Ground'],
        'ISS_dust':    [v for k, v in sample_data.items()
                        if v['env'] == 'ISS' and v['type'] == 'dust'],
        'ISS_hepa':    [v for k, v in sample_data.items()
                        if v['env'] == 'ISS' and v['type'] == 'hepa'],
    }

    def merge_counters(samples, level):
        merged = Counter()
        for s in samples:
            merged += s[level]
        return merged

    # ── 门级组成 ─────────────────────────────────────────────────────
    print('\n### 门（Phylum）级组成对比（归一化 % 总基因）')
    all_phyla = set()
    group_phylum = {}
    for grp_name, samples in groups.items():
        if not samples:
            continue
        merged = merge_counters(samples, 'phylum')
        total  = sum(merged.values())
        group_phylum[grp_name] = (merged, total)
        top = top_taxa(merged, 8, exclude_unclass=True)
        all_phyla.update(k for k, _ in top)

    # 打印表格
    header = f'{"Phylum":<35}'
    for g in groups:
        if groups[g]:
            header += f' {g:>18}'
    print(header)
    print('-' * (35 + 20 * len([g for g in groups if groups[g]])))
    for ph in sorted(all_phyla):
        row = f'{ph:<35}'
        for g in groups:
            if not groups[g]:
                continue
            merged, total = group_phylum[g]
            pct = merged.get(ph, 0) / total * 100 if total else 0
            row += f' {pct:>17.2f}%'
        print(row)

    # ── 属级 Top 15 对比 ─────────────────────────────────────────────
    print('\n### 属（Genus）级 Top15（来自有分类注释的基因，% 总基因）')
    group_genus = {}
    all_genera_top = set()
    for grp_name, samples in groups.items():
        if not samples:
            continue
        merged = merge_counters(samples, 'genus')
        total  = sum(merged.values())
        group_genus[grp_name] = (merged, total)
        top = top_taxa(merged, 15, exclude_unclass=True)
        all_genera_top.update(k for k, _ in top)

    header = f'{"Genus":<30}'
    for g in groups:
        if groups[g]:
            header += f' {g:>20}'
    print(header)
    print('-' * (30 + 22 * len([g for g in groups if groups[g]])))
    # Sort genera by ISS_dust + ISS_hepa combined rank
    iss_combined = Counter()
    for g in ['ISS_dust', 'ISS_hepa']:
        if groups.get(g):
            iss_combined += group_genus[g][0]
    for genus in sorted(all_genera_top, key=lambda x: -iss_combined.get(x, 0)):
        row = f'{genus:<30}'
        for g in groups:
            if not groups[g]:
                continue
            merged, total = group_genus[g]
            pct = merged.get(genus, 0) / total * 100 if total else 0
            row += f' {pct:>19.3f}%'
        print(row)

    # ── Shannon 多样性 ───────────────────────────────────────────────
    print('\n### Shannon 多样性（属级，排除 unclassified）')
    header = f'{"组别":<18} {"样本数":>6} {"Richness":>10} {"Shannon H":>12} {"Evenness":>10}'
    print(header)
    print('-' * 60)
    for grp_name, samples in groups.items():
        if not samples:
            continue
        # 计算每个样本的多样性，取均值
        divs = []
        for s in samples:
            d = compute_diversity(s['genus'])
            divs.append(d)
        mean_R = sum(d['richness'] for d in divs) / len(divs)
        mean_H = sum(d['shannon'] for d in divs) / len(divs)
        mean_E = sum(d['evenness'] for d in divs) / len(divs)
        print(f'{grp_name:<18} {len(samples):>6} {mean_R:>10.0f} {mean_H:>12.3f} {mean_E:>10.3f}')

    # ── 关键发现：Firmicutes vs Proteobacteria ─────────────────────
    print('\n### 关键比较：Firmicutes vs Proteobacteria 比例')
    for grp_name, samples in groups.items():
        if not samples:
            continue
        merged, total = group_phylum[grp_name]
        firm  = merged.get('Firmicutes', 0) / total * 100 if total else 0
        prot  = merged.get('Proteobacteria', 0) / total * 100 if total else 0
        actino = merged.get('Actinobacteria', 0) / total * 100 if total else 0
        print(f'  {grp_name}: Firmicutes={firm:.2f}%  Proteobacteria={prot:.2f}%  '
              f'Actinobacteria={actino:.2f}%')

    # ── 生成报告 ─────────────────────────────────────────────────────
    lines = []
    lines.append('# GLDS-224 宏基因组物种组成分析报告（SESSION-008）')
    lines.append('**分析日期**: 2026-04-17  ')
    lines.append('**数据集**: GLDS-224 (ISS HEPA + 碎屑 vs SAF 地面洁净室)  ')
    lines.append('**注释列**: domain, phylum, class, order, family, genus, species（NASA OSDR TSV）')
    lines.append('')
    lines.append('---')
    lines.append('')
    lines.append('## 样本基本统计')
    lines.append('')
    lines.append('| 样本 | 标签 | 环境 | 类型 | 总基因 | 有物种注释 | 分类覆盖% |')
    lines.append('|------|------|------|------|--------|----------|---------|')
    for srr, data in sample_data.items():
        lines.append(f'| {srr} | {data["label"]} | {data["env"]} | {data["type"]} | '
                     f'{data["total"]:,} | {data["annotated"]:,} | {data["tax_coverage"]:.1f}% |')

    lines.append('')
    lines.append('## 门（Phylum）级组成对比（% 总基因）')
    lines.append('')
    lines.append(f'| Phylum | Ground_dust | ISS_dust | ISS_hepa |')
    lines.append('|--------|------------|---------|---------|')
    all_phyla_sorted = sorted(all_phyla,
        key=lambda x: -(group_phylum.get('ISS_dust', ({}, 1))[0].get(x, 0)
                       + group_phylum.get('ISS_hepa', ({}, 1))[0].get(x, 0)))
    for ph in all_phyla_sorted:
        row_vals = []
        for g in ['Ground_dust', 'ISS_dust', 'ISS_hepa']:
            if g in group_phylum:
                merged, total = group_phylum[g]
                row_vals.append(f'{merged.get(ph,0)/total*100:.2f}%')
            else:
                row_vals.append('N/A')
        lines.append(f'| {ph} | {" | ".join(row_vals)} |')

    lines.append('')
    lines.append('## 属（Genus）级 Top 15（ISS 最丰富）')
    lines.append('')
    lines.append('| Genus | Ground_dust | ISS_dust | ISS_hepa |')
    lines.append('|-------|------------|---------|---------|')
    for genus in sorted(all_genera_top, key=lambda x: -iss_combined.get(x, 0)):
        row_vals = []
        for g in ['Ground_dust', 'ISS_dust', 'ISS_hepa']:
            if g in group_genus:
                merged, total = group_genus[g]
                row_vals.append(f'{merged.get(genus,0)/total*100:.3f}%')
            else:
                row_vals.append('N/A')
        lines.append(f'| {genus} | {" | ".join(row_vals)} |')

    lines.append('')
    lines.append('## 群落多样性（Shannon 指数，属级）')
    lines.append('')
    lines.append('| 组别 | 样本数 | Richness | Shannon H\' | Evenness |')
    lines.append('|------|------|---------|-----------|---------|')
    for grp_name, samples in groups.items():
        if not samples:
            continue
        divs = [compute_diversity(s['genus']) for s in samples]
        mean_R = sum(d['richness'] for d in divs) / len(divs)
        mean_H = sum(d['shannon']  for d in divs) / len(divs)
        mean_E = sum(d['evenness'] for d in divs) / len(divs)
        lines.append(f'| {grp_name} | {len(samples)} | {mean_R:.0f} | {mean_H:.3f} | {mean_E:.3f} |')

    lines.append('')
    lines.append('## 科学解读')
    lines.append('')
    lines.append('### Firmicutes / Proteobacteria 比例')
    for grp_name in ['Ground_dust', 'ISS_dust', 'ISS_hepa']:
        if grp_name not in group_phylum:
            continue
        merged, total = group_phylum[grp_name]
        firm   = merged.get('Firmicutes', 0) / total * 100 if total else 0
        prot   = merged.get('Proteobacteria', 0) / total * 100 if total else 0
        actino = merged.get('Actinobacteria', 0) / total * 100 if total else 0
        lines.append(f'- **{grp_name}**: Firmicutes={firm:.2f}%, '
                     f'Proteobacteria={prot:.2f}%, Actinobacteria={actino:.2f}%')

    lines.append('')
    lines.append('### 与 H-003 的一致性')
    lines.append('')
    lines.append('- 分离株分析（PRJNA637984）发现 ISS Bacillales（Firmicutes）100% 缺失 CRISPR')
    lines.append('- 宏基因组物种组成可验证：ISS 是否富集 Firmicutes 相对于地面')
    lines.append('')
    lines.append('## 方法学注意事项')
    lines.append('')
    lines.append('- 物种注释基于 img/m 的 taxid 分配，覆盖率较低（~15-46%的基因有物种注释）')
    lines.append('- 样本量限制：ISS n=4，Ground n=2，仅描述性分析')
    lines.append('- 宏基因组物种组成反映基因丰度，非菌株丰度')
    lines.append('')
    lines.append('---')
    lines.append('*本报告由 Project Genesis AI 系统自动生成 (SESSION-008)*')

    report_path = REPORT_DIR / 'glds224_taxonomy_report.md'
    report_path.write_text('\n'.join(lines), encoding='utf-8')
    print(f'\n[报告] → {report_path}')

    # Save CSV
    out_csv = OUT_DIR / 'glds224_taxonomy_summary.csv'
    with open(out_csv, 'w', newline='', encoding='utf-8') as f:
        w = csv.writer(f)
        w.writerow(['srr', 'label', 'env', 'type', 'pma',
                    'total_genes', 'annotated_genes', 'tax_coverage_pct',
                    'shannon_genus', 'richness_genus'])
        for srr, data in sample_data.items():
            div = compute_diversity(data['genus'])
            w.writerow([srr, data['label'], data['env'], data['type'], data['pma'],
                        data['total'], data['annotated'], f'{data["tax_coverage"]:.1f}',
                        f'{div["shannon"]:.3f}', div['richness']])
    print(f'[输出] CSV → {out_csv}')
    print('\n完成。')


if __name__ == '__main__':
    main()
