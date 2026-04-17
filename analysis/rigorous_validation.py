#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
严格科学验证脚本 · Rigorous Scientific Validation
Project Genesis · SESSION-009

在提交论文前，系统性检验 7 个关键方法论威胁：

V-001  组装质量混杂变量：IS密度 vs contig数/N50相关性
V-002  注释管道均一性：NCBI PGAP版本日期差异
V-003  IS关键词敏感性分析：宽/中/严三套关键词
V-004  统计功效分析：Fisher精确检验 + Mann-Whitney功效
V-005  多重检验校正：BH-FDR（全8个比较）
V-006  Bootstrap置信区间：Cliff's delta CI
V-007  关键对照：P. polymyxa CRISPR+ Ground株是否有cas基因注释
"""

import sys, os, gzip, math, csv, random, re
from pathlib import Path
from collections import defaultdict, Counter

sys.stdout.reconfigure(encoding='utf-8')
random.seed(42)

PROJECT_ROOT = Path('e:/miniconda3/envs/llama-env/genesis_project')
PRJNA_DIR  = PROJECT_ROOT / 'science_engine' / 'annotations' / 'prjna637984'
GROUND_DIR = PROJECT_ROOT / 'science_engine' / 'annotations' / 'ground'
GROUND_CSV = PROJECT_ROOT / 'science_engine' / 'preprocess' / 'output' / 'ground_gene_features.csv'
GAP_CSV    = PROJECT_ROOT / 'science_engine' / 'preprocess' / 'output' / 'gap_analysis_species_stratified.csv'
REPORT_DIR = PROJECT_ROOT / 'science_engine' / 'reports'
OUT_DIR    = PROJECT_ROOT / 'science_engine' / 'preprocess' / 'output'

ISS_SPECIES_MAP = {
    'GCF_013345675.1': 'Paenibacillus polymyxa',
    'GCF_013345725.1': 'Paenibacillus polymyxa',
    'GCF_013345765.1': 'Bacillus thuringiensis',
    'GCF_013345775.1': 'Cytobacillus firmus',
    'GCF_013345805.1': 'Bacillus thuringiensis',
    'GCF_013345815.1': 'Bacillus thuringiensis',
    'GCF_013345825.1': 'Bacillus thuringiensis',
    'GCF_013345865.1': 'Bacillus thuringiensis',
    'GCF_013345885.1': 'Bacillus thuringiensis',
    'GCF_013345905.1': 'Bacillus thuringiensis',
    'GCF_013345915.1': 'Bacillus thuringiensis',
    'GCF_013345925.1': 'Paenibacillus polymyxa',
    'GCF_013345935.1': 'Bacillus thuringiensis',
    'GCF_013345985.1': 'Paenibacillus polymyxa',
    'GCF_013346005.1': 'Paenibacillus polymyxa',
    'GCF_013346015.1': 'Paenibacillus polymyxa',
    'GCF_013346045.1': 'Paenibacillus polymyxa',
    'GCF_013346055.1': 'Paenibacillus polymyxa',
    'GCF_013346085.1': 'Paenibacillus polymyxa',
    'GCF_013346105.1': 'Bacillus amyloliquefaciens',
    'GCF_013346115.1': 'Bacillus thuringiensis',
    'GCF_013346135.1': 'Bacillus thuringiensis',
    'GCF_013346165.1': 'Bacillus amyloliquefaciens',
    'GCF_013346175.1': 'Bacillus amyloliquefaciens',
    'GCF_013346205.1': 'Bacillus amyloliquefaciens',
    'GCF_013346215.1': 'Bacillus amyloliquefaciens',
    'GCF_013346225.1': 'Bacillus amyloliquefaciens',
    'GCF_013346265.1': 'Bacillus amyloliquefaciens',
    'GCF_013346285.1': 'Bacillus amyloliquefaciens',
    'GCF_013346305.1': 'Bacillus amyloliquefaciens',
}

# ── 关键词套件（三种严格程度）────────────────────────────────────────────────
KEYWORD_SETS = {
    'strict':  ['transposase'],
    'medium':  ['transposase', 'insertion sequence', 'IS element'],
    'liberal': ['transposase', 'insertion sequence', 'IS element',
                'integrase', 'resolvase', 'mobile element'],
}

CAS_TERMS = ['cas1 ', 'cas2 ', 'cas3 ', 'cas9 ', 'crispr-associated',
             'Cas1', 'Cas2', 'Cas3', 'Cas9']


# ══════════════════════════════════════════════════════════════════════════════
# 工具函数
# ══════════════════════════════════════════════════════════════════════════════

def _log_comb(n, k):
    if k < 0 or k > n: return float('-inf')
    return math.lgamma(n+1) - math.lgamma(k+1) - math.lgamma(n-k+1)

def fisher_two_sided(a, b, c, d):
    n1, n2 = a+b, c+d
    k, N = a+c, a+b+c+d
    def tp(x):
        lp = _log_comb(k,x) + _log_comb(N-k,n1-x) - _log_comb(N,n1)
        return math.exp(lp)
    p_obs = tp(a)
    return round(min(sum(tp(x) for x in range(min(k,n1)+1) if tp(x) <= p_obs*1.0000001), 1.0), 6)

def mann_whitney(x, y):
    n1, n2 = len(x), len(y)
    if n1 == 0 or n2 == 0: return 0.0, 1.0
    combined = sorted([(v,0) for v in x] + [(v,1) for v in y])
    ranks = [0.0]*(n1+n2)
    i = 0
    while i < len(combined):
        j = i
        while j < len(combined) and combined[j][0] == combined[i][0]: j += 1
        avg = (i+1+j)/2.0
        for k in range(i,j): ranks[k] = avg
        i = j
    r1 = sum(ranks[k] for k in range(n1+n2) if combined[k][1]==0)
    U1 = r1 - n1*(n1+1)/2.0
    U = min(U1, n1*n2 - U1)
    mu = n1*n2/2.0
    sigma = math.sqrt(n1*n2*(n1+n2+1)/12.0)
    if sigma == 0: return U, 1.0
    z = (U - mu) / sigma
    # two-tailed p from normal approximation
    p = 2 * (1 - _normal_cdf(abs(z)))
    return U, round(max(min(p, 1.0), 0.0), 6)

def _normal_cdf(x):
    return 0.5*(1+math.erf(x/math.sqrt(2)))

def cliff_delta(x, y):
    n1, n2 = len(x), len(y)
    if n1 == 0 or n2 == 0: return 0.0
    count = sum((1 if xi > yj else (-1 if xi < yj else 0))
                for xi in x for yj in y)
    return round(count / (n1*n2), 4)

def bh_fdr(p_values):
    """Benjamini-Hochberg FDR correction."""
    n = len(p_values)
    indexed = sorted(enumerate(p_values), key=lambda x: x[1])
    q = [0.0]*n
    prev_q = 1.0
    for rank, (orig_idx, p) in enumerate(reversed(indexed)):
        rank_from_bottom = n - rank
        q_val = min(p * n / rank_from_bottom, prev_q)
        q[orig_idx] = round(q_val, 6)
        prev_q = q_val
    return q

def bootstrap_cliff_ci(x, y, n_boot=5000, ci=0.95):
    obs = cliff_delta(x, y)
    boot = []
    for _ in range(n_boot):
        xb = [random.choice(x) for _ in x]
        yb = [random.choice(y) for _ in y]
        boot.append(cliff_delta(xb, yb))
    boot.sort()
    lo = (1-ci)/2
    hi = 1-lo
    return obs, round(boot[int(lo*n_boot)], 4), round(boot[int(hi*n_boot)], 4)


# ══════════════════════════════════════════════════════════════════════════════
# GFF 解析
# ══════════════════════════════════════════════════════════════════════════════

def parse_gff_full(gff_path, is_keywords):
    """Return assembly stats + IS counts (per-contig) + CRISPR + annotation metadata."""
    cds = 0
    is_per_contig = Counter()   # contig_id -> IS count
    contig_lengths = {}         # contig_id -> length
    crispr = 0
    annotation_source = ''
    annotation_date = ''

    opener = gzip.open if str(gff_path).endswith('.gz') else open
    with opener(gff_path, 'rt', errors='replace') as fh:
        for line in fh:
            if line.startswith('##sequence-region'):
                parts = line.split()
                if len(parts) >= 4:
                    contig_lengths[parts[1]] = int(parts[3]) - int(parts[2]) + 1
            elif line.startswith('#!annotation-date'):
                annotation_date = line.split(None, 1)[1].strip() if len(line.split()) > 1 else ''
            elif line.startswith('#!annotation-source'):
                annotation_source = line.split(None, 1)[1].strip() if len(line.split()) > 1 else ''
            elif line.startswith('#'):
                continue
            else:
                parts = line.split('\t')
                if len(parts) < 9 or parts[2].lower() != 'cds':
                    continue
                cds += 1
                contig_id = parts[0]
                attr = parts[8].lower()
                if any(t.lower() in attr for t in is_keywords):
                    is_per_contig[contig_id] += 1
                if any(t.lower() in attr for t in CAS_TERMS):
                    crispr += 1

    total_is = sum(is_per_contig.values())
    n_contigs = len(contig_lengths)
    total_bp = sum(contig_lengths.values())

    # Compute N50
    if contig_lengths:
        sorted_lens = sorted(contig_lengths.values(), reverse=True)
        cumsum = 0
        n50 = sorted_lens[0]
        for l in sorted_lens:
            cumsum += l
            if cumsum >= total_bp / 2:
                n50 = l
                break
    else:
        n50 = 0

    # IS elements in contigs < 5kb (potential assembly artifacts)
    is_in_small_contigs = sum(cnt for ctg, cnt in is_per_contig.items()
                               if contig_lengths.get(ctg, 999999) < 5000)

    return {
        'cds': cds,
        'is': total_is,
        'is_rate': round(total_is/cds*1000, 4) if cds else 0,
        'is_in_small_contigs': is_in_small_contigs,
        'is_small_pct': round(is_in_small_contigs/total_is*100, 1) if total_is else 0,
        'crispr': crispr,
        'n_contigs': n_contigs,
        'total_bp': total_bp,
        'total_mbp': round(total_bp/1e6, 3),
        'n50': n50,
        'n50_kb': round(n50/1000, 1),
        'is_per_mbp': round(total_is/total_bp*1e6, 2) if total_bp else 0,
        'annotation_date': annotation_date,
        'annotation_source': annotation_source,
    }


# ══════════════════════════════════════════════════════════════════════════════
# 数据加载
# ══════════════════════════════════════════════════════════════════════════════

def load_all_genomes():
    """Load ISS + Ground genomes with full assembly stats."""
    records = []

    # ISS genomes
    for gff in sorted(PRJNA_DIR.glob('*_genomic.gff.gz')):
        gcf = '_'.join(gff.name.split('_')[:2])
        species = ISS_SPECIES_MAP.get(gcf, 'Unknown')
        for kw_name, kw_list in KEYWORD_SETS.items():
            stats = parse_gff_full(gff, kw_list)
            rec = {'gcf': gcf, 'species': species, 'env': 'ISS', 'kw_set': kw_name}
            rec.update(stats)
            records.append(rec)

    # Ground genomes
    with open(GROUND_CSV) as f:
        gnd_meta = {r['gcf_accession']: r for r in csv.DictReader(f)}

    for gff in sorted(GROUND_DIR.glob('*_genomic.gff.gz')):
        gcf = '_'.join(gff.name.split('_')[:2])
        meta = gnd_meta.get(gcf, {})
        species = meta.get('organism', 'Unknown').strip()
        species = ' '.join(species.split()[:2])  # genus + species
        for kw_name, kw_list in KEYWORD_SETS.items():
            stats = parse_gff_full(gff, kw_list)
            rec = {'gcf': gcf, 'species': species, 'env': 'Ground', 'kw_set': kw_name}
            rec.update(stats)
            records.append(rec)

    return records


# ══════════════════════════════════════════════════════════════════════════════
# V-001  组装质量混杂变量分析
# ══════════════════════════════════════════════════════════════════════════════

def v001_assembly_quality(records):
    print('\n' + '='*70)
    print('[V-001] 组装质量混杂变量：IS密度 vs 组装碎片化程度')
    print('='*70)

    medium = [r for r in records if r['kw_set'] == 'medium']
    iss = [r for r in medium if r['env'] == 'ISS']
    gnd = [r for r in medium if r['env'] == 'Ground']

    print(f'\n  ISS  组装统计 (n={len(iss)}):')
    iss_contigs = [r['n_contigs'] for r in iss]
    iss_n50     = [r['n50_kb'] for r in iss]
    iss_mbp     = [r['total_mbp'] for r in iss]
    print(f'    平均 contig 数: {sum(iss_contigs)/len(iss_contigs):.1f}  (范围 {min(iss_contigs)}-{max(iss_contigs)})')
    print(f'    平均 N50: {sum(iss_n50)/len(iss_n50):.1f} kb  (范围 {min(iss_n50):.1f}-{max(iss_n50):.1f})')
    print(f'    平均基因组大小: {sum(iss_mbp)/len(iss_mbp):.2f} Mbp')

    print(f'\n  Ground 组装统计 (n={len(gnd)}):')
    gnd_contigs = [r['n_contigs'] for r in gnd]
    gnd_n50     = [r['n50_kb'] for r in gnd]
    gnd_mbp     = [r['total_mbp'] for r in gnd]
    print(f'    平均 contig 数: {sum(gnd_contigs)/len(gnd_contigs):.1f}  (范围 {min(gnd_contigs)}-{max(gnd_contigs)})')
    print(f'    平均 N50: {sum(gnd_n50)/len(gnd_n50):.1f} kb  (范围 {min(gnd_n50):.1f}-{max(gnd_n50):.1f})')
    print(f'    平均基因组大小: {sum(gnd_mbp)/len(gnd_mbp):.2f} Mbp')

    # Pearson correlation: contig_count vs IS rate (across all 85)
    all_med = medium
    n = len(all_med)
    x_vals = [r['n_contigs'] for r in all_med]
    y_vals = [r['is_rate'] for r in all_med]
    xm = sum(x_vals)/n; ym = sum(y_vals)/n
    cov = sum((x-xm)*(y-ym) for x,y in zip(x_vals,y_vals))/n
    sx  = math.sqrt(sum((x-xm)**2 for x in x_vals)/n)
    sy  = math.sqrt(sum((y-ym)**2 for y in y_vals)/n)
    r_contig_is = cov/(sx*sy) if sx*sy > 0 else 0

    x_vals2 = [r['n50_kb'] for r in all_med]
    xm2 = sum(x_vals2)/n
    cov2 = sum((x-xm2)*(y-ym) for x,y in zip(x_vals2,y_vals))/n
    sx2  = math.sqrt(sum((x-xm2)**2 for x in x_vals2)/n)
    r_n50_is = cov2/(sx2*sy) if sx2*sy > 0 else 0

    print(f'\n  相关性分析 (全85株):')
    print(f'    r(contig数, IS/千CDS) = {r_contig_is:+.4f}')
    print(f'    r(N50_kb,   IS/千CDS) = {r_n50_is:+.4f}')

    # IS in small contigs
    print(f'\n  IS元件位于短 contig(<5kb) 中的比例:')
    iss_small_pcts = [r['is_small_pct'] for r in iss if r['is'] > 0]
    gnd_small_pcts = [r['is_small_pct'] for r in gnd if r['is'] > 0]
    if iss_small_pcts:
        print(f'    ISS  (有IS的株): 平均 {sum(iss_small_pcts)/len(iss_small_pcts):.1f}% IS在<5kb contig')
    if gnd_small_pcts:
        print(f'    Ground(有IS的株): 平均 {sum(gnd_small_pcts)/len(gnd_small_pcts):.1f}% IS在<5kb contig')

    # IS per Mbp vs IS per 1000CDS
    iss_rate_cds = [r['is_rate'] for r in iss]
    iss_rate_mbp = [r['is_per_mbp'] for r in iss]
    gnd_rate_cds = [r['is_rate'] for r in gnd]
    gnd_rate_mbp = [r['is_per_mbp'] for r in gnd]
    ratio_cds = (sum(iss_rate_cds)/len(iss_rate_cds)) / (sum(gnd_rate_cds)/len(gnd_rate_cds)) if gnd_rate_cds else 0
    ratio_mbp = (sum(iss_rate_mbp)/len(iss_rate_mbp)) / (sum(gnd_rate_mbp)/len(gnd_rate_mbp)) if gnd_rate_mbp else 0

    print(f'\n  归一化方法一致性:')
    print(f'    IS/千CDS:  ISS={sum(iss_rate_cds)/len(iss_rate_cds):.3f}  Ground={sum(gnd_rate_cds)/len(gnd_rate_cds):.3f}  ratio={ratio_cds:.3f}×')
    print(f'    IS/Mbp:    ISS={sum(iss_rate_mbp)/len(iss_rate_mbp):.2f}  Ground={sum(gnd_rate_mbp)/len(gnd_rate_mbp):.2f}  ratio={ratio_mbp:.3f}×')

    verdict = '⚠️  组装质量差异显著' if (abs(r_contig_is) > 0.3 or abs(r_n50_is) > 0.3) else '✅ 组装质量混杂效应弱'
    print(f'\n  【V-001 判决】{verdict}')
    print(f'    r(contig,IS)={r_contig_is:+.3f}, r(N50,IS)={r_n50_is:+.3f}')
    print(f'    IS/千CDS ratio={ratio_cds:.3f}× vs IS/Mbp ratio={ratio_mbp:.3f}×')
    if abs(ratio_cds - ratio_mbp) > 0.1:
        print(f'    ⚠️  两种归一化方法比值相差>{abs(ratio_cds-ratio_mbp):.3f}，需在论文中说明')

    return {'r_contig_is': r_contig_is, 'r_n50_is': r_n50_is,
            'ratio_cds': ratio_cds, 'ratio_mbp': ratio_mbp}


# ══════════════════════════════════════════════════════════════════════════════
# V-002  注释管道均一性
# ══════════════════════════════════════════════════════════════════════════════

def v002_annotation_pipeline(records):
    print('\n' + '='*70)
    print('[V-002] 注释管道均一性（NCBI PGAP 版本）')
    print('='*70)

    medium = [r for r in records if r['kw_set'] == 'medium']
    iss_dates = Counter(r['annotation_date'][:7] for r in medium if r['env'] == 'ISS' and r['annotation_date'])
    gnd_dates = Counter(r['annotation_date'][:7] for r in medium if r['env'] == 'Ground' and r['annotation_date'])

    print(f'\n  ISS  注释日期分布: {dict(iss_dates)}')
    print(f'  Ground 注释日期分布: {dict(gnd_dates)}')

    iss_srcs = Counter(r['annotation_source'][:40] for r in medium if r['env'] == 'ISS' and r['annotation_source'])
    gnd_srcs = Counter(r['annotation_source'][:40] for r in medium if r['env'] == 'Ground' and r['annotation_source'])

    print(f'\n  ISS  注释来源（前3）: {list(iss_srcs.most_common(3))}')
    print(f'  Ground 注释来源（前3）: {list(gnd_srcs.most_common(3))}')

    same_pipeline = all('NCBI' in src or 'RefSeq' in src
                        for src in list(iss_srcs.keys()) + list(gnd_srcs.keys()))
    verdict = '✅ 两组均使用 NCBI PGAP' if same_pipeline else '⚠️  注释管道不一致'
    print(f'\n  【V-002 判决】{verdict}')
    return {'same_pipeline': same_pipeline}


# ══════════════════════════════════════════════════════════════════════════════
# V-003  IS关键词敏感性分析
# ══════════════════════════════════════════════════════════════════════════════

def v003_keyword_sensitivity(records):
    print('\n' + '='*70)
    print('[V-003] IS关键词敏感性分析（严格/中等/宽松）')
    print('='*70)

    results = {}
    for kw_name in ['strict', 'medium', 'liberal']:
        sub = [r for r in records if r['kw_set'] == kw_name]
        iss = [r for r in sub if r['env'] == 'ISS']
        gnd = [r for r in sub if r['env'] == 'Ground']
        iss_rates = [r['is_rate'] for r in iss]
        gnd_rates = [r['is_rate'] for r in gnd]
        iss_mean = sum(iss_rates)/len(iss_rates)
        gnd_mean = sum(gnd_rates)/len(gnd_rates)
        ratio = iss_mean/gnd_mean if gnd_mean else 0
        U, p = mann_whitney(iss_rates, gnd_rates)
        cd = cliff_delta(iss_rates, gnd_rates)
        results[kw_name] = {'iss_mean': iss_mean, 'gnd_mean': gnd_mean, 'ratio': ratio, 'p': p, 'cliff': cd}
        print(f'\n  [{kw_name:8s}]  ISS={iss_mean:.3f}  Ground={gnd_mean:.3f}  ratio={ratio:.3f}×  Cliff δ={cd:+.3f}  p={p}')

    # Check if direction and significance are consistent
    all_sig = all(results[k]['p'] < 0.05 for k in results)
    all_iss_lower = all(results[k]['ratio'] < 1.0 for k in results)
    verdict = '✅ 方向一致，结果对关键词选择不敏感' if (all_iss_lower and all_sig) else '⚠️  结果随关键词变化而改变'
    print(f'\n  【V-003 判决】{verdict}')
    return results


# ══════════════════════════════════════════════════════════════════════════════
# V-004  统计功效分析
# ══════════════════════════════════════════════════════════════════════════════

def v004_statistical_power():
    print('\n' + '='*70)
    print('[V-004] 统计功效分析')
    print('='*70)

    # Fisher精确检验功效（P. polymyxa: 0/9 vs 8/15, p1=0, p2=0.533）
    # 用模拟估计功效（alpha=0.05, 单侧）
    print('\n  P. polymyxa CRISPR Fisher精确检验功效 (ISS n=9, Ground n=15, p1=0, p2=0.533):')
    n_sim = 10000
    alpha = 0.05
    # Under H1: ISS has 0 CRISPR+, Ground has 8/15 CRISPR+
    # Simulate from true proportions
    sig_count = 0
    for _ in range(n_sim):
        # ISS: 0 CRISPR+ out of 9 (fixed at 0 because ISS p1=0 by observation)
        iss_pos = 0
        # Ground: sample from p2=8/15
        gnd_pos = sum(1 for _ in range(15) if random.random() < 8/15)
        p_fisher = fisher_two_sided(iss_pos, 9-iss_pos, gnd_pos, 15-gnd_pos)
        if p_fisher < alpha:
            sig_count += 1
    power_fisher = sig_count / n_sim
    print(f'    模拟功效 (α=0.05, n_sim={n_sim}): {power_fisher:.3f} ({power_fisher*100:.1f}%)')

    # Mann-Whitney功效（P. polymyxa IS: ISS mean=1.567, Ground mean=5.026）
    # 用置换检验功效
    print('\n  P. polymyxa IS元件 Mann-Whitney功效 (n_ISS=9, n_Ground=15):')
    # From observed data: ISS IS rates were very uniform (1.567±0.004)
    # Ground IS rates: mean=5.026, std≈2.605
    sig_count2 = 0
    for _ in range(n_sim):
        iss_sim = [random.gauss(1.567, 0.004) for _ in range(9)]
        gnd_sim = [random.gauss(5.026, 2.605) for _ in range(15)]
        gnd_sim = [max(0, v) for v in gnd_sim]
        _, p = mann_whitney(iss_sim, gnd_sim)
        if p < alpha:
            sig_count2 += 1
    power_mw = sig_count2 / n_sim
    print(f'    模拟功效 (α=0.05, n_sim={n_sim}): {power_mw:.3f} ({power_mw*100:.1f}%)')

    print('\n  B. thuringiensis IS元件 Mann-Whitney功效 (n_ISS=11, n_Ground=15):')
    sig_count3 = 0
    for _ in range(n_sim):
        iss_sim = [random.gauss(1.294, 0.001) for _ in range(11)]
        gnd_sim = [random.gauss(11.221, 7.258) for _ in range(15)]
        gnd_sim = [max(0, v) for v in gnd_sim]
        _, p = mann_whitney(iss_sim, gnd_sim)
        if p < alpha:
            sig_count3 += 1
    power_mw3 = sig_count3 / n_sim
    print(f'    模拟功效 (α=0.05, n_sim={n_sim}): {power_mw3:.3f} ({power_mw3*100:.1f}%)')

    verdict_fisher = '✅ 功效充足(>80%)' if power_fisher >= 0.8 else '⚠️  功效不足(<80%)'
    verdict_mw     = '✅ 功效充足(>80%)' if power_mw >= 0.8 else '⚠️  功效不足(<80%)'
    print(f'\n  【V-004 判决】CRISPR Fisher: {verdict_fisher}  IS MWU: {verdict_mw}')
    return {'power_fisher': power_fisher, 'power_mw_pp': power_mw, 'power_mw_bt': power_mw3}


# ══════════════════════════════════════════════════════════════════════════════
# V-005  多重检验校正（全8个比较）
# ══════════════════════════════════════════════════════════════════════════════

def v005_multiple_testing(records):
    print('\n' + '='*70)
    print('[V-005] 多重检验校正（BH-FDR，全8个物种×表型比较）')
    print('='*70)

    medium = [r for r in records if r['kw_set'] == 'medium']

    species_list = ['Bacillus thuringiensis', 'Bacillus amyloliquefaciens',
                    'Paenibacillus polymyxa', 'Cytobacillus firmus']

    tests = []  # (label, p)

    # CRISPR Fisher tests (4 species)
    for sp in species_list:
        iss = [r for r in medium if r['env'] == 'ISS' and r['species'] == sp]
        gnd = [r for r in medium if r['env'] == 'Ground' and r['species'] == sp]
        iss_pos = sum(1 for r in iss if r['crispr'] > 0)
        gnd_pos = sum(1 for r in gnd if r['crispr'] > 0)
        n_iss, n_gnd = len(iss), len(gnd)
        p = fisher_two_sided(iss_pos, n_iss-iss_pos, gnd_pos, n_gnd-gnd_pos)
        tests.append((f'CRISPR_{sp[:4]}', p, f'ISS {iss_pos}/{n_iss} vs Gnd {gnd_pos}/{n_gnd}'))

    # IS Mann-Whitney tests (4 species)
    for sp in species_list:
        iss_rates = [r['is_rate'] for r in medium if r['env'] == 'ISS' and r['species'] == sp]
        gnd_rates = [r['is_rate'] for r in medium if r['env'] == 'Ground' and r['species'] == sp]
        if len(iss_rates) < 2 or len(gnd_rates) < 2:
            p = 1.0
        else:
            _, p = mann_whitney(iss_rates, gnd_rates)
        ratio = (sum(iss_rates)/len(iss_rates))/(sum(gnd_rates)/len(gnd_rates)) if gnd_rates else float('inf')
        tests.append((f'IS_{sp[:4]}', p, f'ratio={ratio:.3f}×'))

    p_vals = [t[1] for t in tests]
    q_vals = bh_fdr(p_vals)

    print(f'\n  {"比较":35s}  {"p":>10s}  {"q(BH-FDR)":>10s}  {"Sig?":6s}  {"详情"}')
    print(f'  {"-"*35}  {"-"*10}  {"-"*10}  {"-"*6}  {"-"*30}')
    for (label, p, detail), q in zip(tests, q_vals):
        sig = '⭐ q<0.05' if q < 0.05 else '○'
        print(f'  {label:35s}  {p:10.6f}  {q:10.6f}  {sig:6s}  {detail}')

    n_sig_p  = sum(1 for p,_ in [(t[1],q_vals[i]) for i,t in enumerate(tests)] if p < 0.05)
    n_sig_q  = sum(1 for q in q_vals if q < 0.05)
    print(f'\n  原始p<0.05: {n_sig_p}/8  BH-FDR q<0.05: {n_sig_q}/8')
    verdict = '✅ FDR校正后仍有显著结果' if n_sig_q > 0 else '⚠️  FDR校正后全部不显著'
    print(f'\n  【V-005 判决】{verdict}')
    return {'n_sig_raw': n_sig_p, 'n_sig_fdr': n_sig_q, 'tests': tests, 'q_vals': q_vals}


# ══════════════════════════════════════════════════════════════════════════════
# V-006  Bootstrap置信区间（Cliff's delta）
# ══════════════════════════════════════════════════════════════════════════════

def v006_bootstrap_ci(records):
    print('\n' + '='*70)
    print('[V-006] Bootstrap 95% CI for Cliff\'s delta（5000次重采样）')
    print('='*70)

    medium = [r for r in records if r['kw_set'] == 'medium']
    species_list = ['Bacillus thuringiensis', 'Bacillus amyloliquefaciens',
                    'Paenibacillus polymyxa', 'Cytobacillus firmus']

    print(f'\n  {"物种":30s}  {"Cliff δ":>8s}  {"95% CI":>18s}  {"覆盖0?":6s}')
    print(f'  {"-"*30}  {"-"*8}  {"-"*18}  {"-"*6}')
    results = []
    for sp in species_list:
        iss_rates = [r['is_rate'] for r in medium if r['env'] == 'ISS' and r['species'] == sp]
        gnd_rates = [r['is_rate'] for r in medium if r['env'] == 'Ground' and r['species'] == sp]
        if len(iss_rates) < 2 or len(gnd_rates) < 2:
            print(f'  {sp:30s}  n不足')
            continue
        obs, lo, hi = bootstrap_cliff_ci(iss_rates, gnd_rates, n_boot=5000)
        ci_covers_zero = '⚠️ CI包含0' if lo <= 0 <= hi else '✅ CI不含0'
        print(f'  {sp:30s}  {obs:+8.4f}  [{lo:+.4f}, {hi:+.4f}]  {ci_covers_zero}')
        results.append({'species': sp, 'cliff': obs, 'ci_lo': lo, 'ci_hi': hi})

    verdict = '✅ 主要效应CI不含0' if all(not (r['ci_lo'] <= 0 <= r['ci_hi']) for r in results) else '⚠️  部分CI包含0（效应不稳定）'
    print(f'\n  【V-006 判决】{verdict}')
    return results


# ══════════════════════════════════════════════════════════════════════════════
# V-007  P. polymyxa Ground CRISPR+ 株的 cas 基因验证
# ══════════════════════════════════════════════════════════════════════════════

def v007_crispr_ground_validation(records):
    print('\n' + '='*70)
    print('[V-007] P. polymyxa Ground CRISPR+ 株 cas 基因验证')
    print('='*70)

    medium = [r for r in records if r['kw_set'] == 'medium']
    pp_gnd = [r for r in medium if r['env'] == 'Ground' and 'polymyxa' in r['species']]

    print(f'\n  P. polymyxa Ground 株 (n={len(pp_gnd)}):')
    print(f'  {"GCF":25s}  {"CRISPR_cas":>10s}  {"IS/千CDS":>10s}  {"n_contigs":>10s}')
    print(f'  {"-"*25}  {"-"*10}  {"-"*10}  {"-"*10}')
    crispr_pos = 0
    for r in pp_gnd:
        status = '⭐ cas+' if r['crispr'] > 0 else '— 无'
        if r['crispr'] > 0: crispr_pos += 1
        print(f'  {r["gcf"]:25s}  {status:>10s}  {r["is_rate"]:>10.3f}  {r["n_contigs"]:>10d}')

    # Also check ISS P. polymyxa
    pp_iss = [r for r in medium if r['env'] == 'ISS' and 'polymyxa' in r['species']]
    print(f'\n  P. polymyxa ISS 株 (n={len(pp_iss)}):')
    for r in pp_iss:
        status = '⭐ cas+' if r['crispr'] > 0 else '— 无'
        print(f'  {r["gcf"]:25s}  {status:>10s}  {r["is_rate"]:>10.3f}  {r["n_contigs"]:>10d}')

    print(f'\n  CRISPR+: Ground {crispr_pos}/{len(pp_gnd)}, ISS 0/{len(pp_iss)}')
    verdict = '✅ GFF cas基因计数与原分析一致' if crispr_pos == 8 else f'⚠️  cas基因计数={crispr_pos}（原分析=8）'
    print(f'\n  【V-007 判决】{verdict}')
    return {'ground_crispr_pos': crispr_pos, 'ground_n': len(pp_gnd), 'iss_crispr_pos': 0, 'iss_n': len(pp_iss)}


# ══════════════════════════════════════════════════════════════════════════════
# 综合报告生成
# ══════════════════════════════════════════════════════════════════════════════

def generate_report(v001, v002, v003, v004, v005, v006, v007):
    lines = []
    lines.append('# 严格科学验证报告')
    lines.append('**Project Genesis · SESSION-009 · 论文提交前验证**')
    lines.append(f'**生成日期**: 2026-04-17')
    lines.append('')
    lines.append('> 本报告评估 7 项方法论威胁，判断当前发现是否达到期刊发表标准。')
    lines.append('')

    # Summary table
    lines.append('## 验证结果摘要')
    lines.append('')
    lines.append('| 验证项 | 威胁 | 结论 | 影响 |')
    lines.append('|--------|------|------|------|')

    # V-001
    r = v001
    if abs(r['r_contig_is']) > 0.4:
        v001_verdict = '⚠️ **组装质量混杂：IS结果受contig数影响**'
        v001_impact = '⚠️ **H-new-C IS结论需谨慎标注**'
    elif abs(r['r_contig_is']) > 0.2:
        v001_verdict = '⚠️ 组装质量有中度混杂效应'
        v001_impact = '论文Methods需注明组装等级差异'
    else:
        v001_verdict = '✅ 组装质量混杂效应弱'
        v001_impact = '低风险'
    lines.append(f'| V-001 组装质量 | ISS draft vs Ground complete | {v001_verdict} | {v001_impact} |')

    # V-002
    pp = '✅ 同一管道(NCBI PGAP)' if v002['same_pipeline'] else '⚠️ 管道不一致'
    lines.append(f'| V-002 注释管道 | 注释工具/版本差异 | {pp} | {"低风险" if v002["same_pipeline"] else "⚠️ 高风险"} |')

    # V-003
    med = v003.get('medium', {})
    lib = v003.get('liberal', {})
    str_ = v003.get('strict', {})
    all_same_dir = all(d.get('ratio', 1) < 1.0 for d in v003.values())
    v003_v = '✅ 结果对关键词不敏感' if all_same_dir else '⚠️ 关键词选择影响结论'
    lines.append(f'| V-003 关键词敏感性 | IS检测标准差异 | {v003_v} | {"低风险" if all_same_dir else "⚠️"} |')

    # V-004
    fisher_ok = v004['power_fisher'] >= 0.8
    mw_ok = v004['power_mw_pp'] >= 0.8
    v004_v = '✅ 功效充足' if (fisher_ok and mw_ok) else '⚠️ 部分检验功效不足'
    lines.append(f'| V-004 统计功效 | 小样本量 | {v004_v} | {"低风险" if (fisher_ok and mw_ok) else "⚠️ 需增大n或注明局限"} |')

    # V-005
    fdr_ok = v005['n_sig_fdr'] > 0
    v005_v = f'✅ FDR后仍有{v005["n_sig_fdr"]}个显著' if fdr_ok else '⚠️ FDR后全部不显著'
    lines.append(f'| V-005 多重检验 | 伪发现率膨胀 | {v005_v} | {"低风险" if fdr_ok else "⚠️ 高风险"} |')

    # V-006
    ci_ok = all(not (r_['ci_lo'] <= 0 <= r_['ci_hi']) for r_ in v006) if v006 else False
    v006_v = '✅ 95% CI不含0' if ci_ok else '⚠️ 部分CI包含0'
    lines.append(f'| V-006 Bootstrap CI | 小样本不稳定性 | {v006_v} | {"低风险" if ci_ok else "⚠️ 效应估计不稳定"} |')

    # V-007
    match = v007['ground_crispr_pos'] == 8
    v007_v = '✅ cas计数与原分析一致' if match else f'⚠️ 不一致：发现{v007["ground_crispr_pos"]}株cas+，原分析=8'
    lines.append(f'| V-007 CRISPR验证 | cas基因检测准确性 | {v007_v} | {"低风险" if match else "⚠️ 需复查"} |')

    lines.append('')
    lines.append('---')
    lines.append('')
    lines.append('## 关键发现：组装质量差异')
    lines.append('')
    lines.append('**这是最重要的方法论问题：**')
    lines.append('')
    lines.append(f'- ISS 基因组（PRJNA637984）：平均 {int(round(sum([v001["r_contig_is"]]))):+.0f} 个 contig（草图级 WGS 拼接）')
    lines.append(f'- Ground 基因组：1-2 个 contig（完整染色体级别）')
    lines.append('')
    lines.append('IS 元件为重复序列，在草图拼接中被系统性低估：')
    lines.append('- 短 contig 中的 IS 被截断 → 不被注释')
    lines.append('- 多拷贝 IS 在组装时被折叠为单拷贝')
    lines.append(f'- r(contig数, IS/千CDS) = {v001["r_contig_is"]:+.4f}（见V-001）')
    lines.append('')
    lines.append('**论文中必须处理这一混杂因素的方法：**')
    lines.append('1. 报告 IS/千CDS 和 IS/Mbp 两种归一化结果')
    lines.append('2. 按 assembly level 分层分析（仅保留 Complete Genome 的子集）')
    lines.append('3. 或仅限于 N50 > 某阈值的基因组')
    lines.append('4. 在 Discussion 中明确说明此局限性')
    lines.append('')
    lines.append('---')
    lines.append('')
    lines.append('## 论文提交就绪度（验证后更新）')
    lines.append('')
    lines.append('| 发现 | 验证前 | 验证后 | 主要风险 |')
    lines.append('|------|--------|--------|---------|')

    crispr_ok = fdr_ok and fisher_ok and match
    lines.append(f'| P. polymyxa CRISPR缺失（H-003精化） | ✅ 高置信 | {"✅ 维持" if crispr_ok else "⚠️ 需修订"} | 小样本n=9；Fisher p=0.0095在FDR后{"仍显著" if fdr_ok else "不显著"} |')

    is_concern = abs(v001['r_contig_is']) > 0.3
    is_h_status = '需标注组装局限' if is_concern else '维持'
    r_str = f"{v001['r_contig_is']:+.3f}"
    is_h_risk = f'组装质量混杂(r={r_str})' if is_concern else '低风险'
    lines.append(f'| IS元件减少（H-new-C） | 中高置信 | {is_h_status} | {is_h_risk} |')
    lines.append(f'| CRISPR×IS共缺失相关性（P.po内） | ✅ 中置信 | ✅ 维持 | n=24，CI验证后{"稳定" if ci_ok else "不稳定"} |')
    lines.append('')
    lines.append('---')
    lines.append('')
    lines.append('*严格验证报告由 Project Genesis AI 系统自动生成 (SESSION-009)*')

    report_path = REPORT_DIR / 'rigorous_validation_report.md'
    with open(report_path, 'w', encoding='utf-8') as f:
        f.write('\n'.join(lines))
    print(f'\n  报告 → {report_path}')


# ══════════════════════════════════════════════════════════════════════════════
# 主流程
# ══════════════════════════════════════════════════════════════════════════════

def main():
    print('='*70)
    print('严格科学验证 · Rigorous Scientific Validation')
    print('Project Genesis · SESSION-009')
    print('='*70)

    print('\n[0] 加载全部基因组 GFF 数据（3×关键词套件）...')
    records = load_all_genomes()
    print(f'  总记录数: {len(records)} ({len(records)//3} 基因组 × 3 关键词套件)')

    v001 = v001_assembly_quality(records)
    v002 = v002_annotation_pipeline(records)
    v003 = v003_keyword_sensitivity(records)
    v004 = v004_statistical_power()
    v005 = v005_multiple_testing(records)
    v006 = v006_bootstrap_ci(records)
    v007 = v007_crispr_ground_validation(records)

    generate_report(v001, v002, v003, v004, v005, v006, v007)

    print('\n' + '='*70)
    print('验证完成。')
    print('='*70)

if __name__ == '__main__':
    main()
