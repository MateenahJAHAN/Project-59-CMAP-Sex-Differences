"""Pilot, unfiltered screen, sampling audit and FAERS cohort comparison."""
import json
import hashlib
import logging
import os
import re
from pathlib import Path
import sys

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, spearmanr

REPO = Path(__file__).resolve().parents[1]
PROOF = REPO/'verification/source_replay/phases2_5'


def save_json(path, obj):
    path.write_text(json.dumps(obj, indent=2, allow_nan=False)+'\n')
    print(path.name, obj, flush=True)


def screen_inputs():
    sys.path.insert(0, str(PROOF))
    import replay
    replay.init()
    return replay


def pilot(out):
    r = screen_inputs()
    sub = r.known[r.known.cmap_name == 'vorinostat']
    matrix = r.expr(sub.sig_id)
    tests = r.mwu(matrix, sub)
    annot = sub.set_index('sig_id').donor_sex.loc[matrix.columns]
    for sex, name in [('M','male_mean'), ('F','female_mean')]:
        values = matrix.loc[:,annot == sex].to_numpy()
        tests[name] = [row.mean() for row in values]
    tests.to_csv(out/'gene_tests.csv.gz', index=False)
    archived = pd.read_csv(REPO/'Phase1_Vorinostat/vorinostat_all_genes_results.csv',dtype={'gene_id':str})
    both = tests.merge(archived, on='gene_id',suffixes=('_run','_archive'),validate='one_to_one')
    assert len(both) == len(tests) == len(archived) == 12328
    hits = tests[tests.p_adj < .05]
    strong = hits.loc[hits.difference.abs() > 1, 'gene_id']
    archived_strong = archived.loc[(archived.p_adj < .05)&(archived.difference.abs() > 1),'gene_id']
    save_json(out/'comparison.json', dict(signatures=len(sub), lines=sub.cell_iname.nunique(),
              gene_tests=len(tests), adjusted_hits=len(hits), above_effect_threshold=len(strong),
              changed_strong_gene_ids=sorted(set(strong)^set(archived_strong)),
              max_absolute_differences={col:float((both[col+'_run']-both[col+'_archive']).abs().max())
                   for col in ['p_value','p_adj','male_mean','female_mean','difference']}))


def unfiltered(out):
    r = screen_inputs()
    if r.gseapy.__version__ != '1.3.2':
        raise RuntimeError('Phase 6 requires GSEApy 1.3.2. Run fellowship/install.sh or use --phase6-python.')
    logging.disable(logging.WARNING)
    rows = []
    names = json.loads((REPO/'Phase5_PerSignature_Categories/reference/top50_names.json').read_text())
    for pos, drug in enumerate(names, 1):
        for lineage, sub in r.known[r.known.cmap_name == drug].groupby('cell_lineage'):
            if lineage == 'unknown':
                continue
            home = set(r.map5.get(lineage, []))
            reason = 'no_GTEx_home_tissue' if not home else 'too_few_signatures_or_DE_genes'
            values = {}
            ranked_count = None
            n_lines = sub.groupby('donor_sex').cell_iname.nunique().to_dict()
            if home and len(sub) >= 3:
                # The published Phase 6 replay presented columns in metadata order to pandas.
                idx = [r.positions[s] for s in sub.sig_id]
                arr = np.asarray(r.Z[sorted(idx)])[np.argsort(np.argsort(idx))]
                mean = pd.DataFrame(arr.T, index=r.gene_ids).mean(axis=1)
                strong = mean[mean.abs()>2]
                if len(strong):
                    de = pd.DataFrame({'ensembl':strong.index.map(r.gene_map),'score':strong.values})
                    de = de.sort_values('score',ascending=False).dropna(subset=['ensembl'])
                    rank = de.groupby('ensembl').score.mean().sort_values(ascending=False)
                    ranked_count = len(rank)
                    sets = r.minimum_sets(rank)
                    reason = 'gsea_empty'
                    if sets:
                        fit = r.gseapy.prerank(rnk=rank, gene_sets=sets, permutation_num=1000,
                            min_size=15, max_size=5000, seed=42, no_plot=True, outdir=None, threads=1)
                        values = fit.res2d.drop_duplicates('Term').set_index('Term').to_dict('index')
                        reason = 'gene_set_too_small'
            for tissue in sorted(r.union):
                item = values.get(tissue)
                is_home = tissue in home
                cat = 'D' if item is None else ('A' if is_home else 'B') if float(item['FDR q-val'])<.05 else 'C'
                rows.append(dict(drug=drug,lineage=lineage,tissue=tissue,is_home_tissue=is_home,
                    n_signatures=len(sub),n_M=n_lines.get('M',0),n_F=n_lines.get('F',0),n_ranked=ranked_count,
                    NES=item['NES'] if item else None,FDR=item['FDR q-val'] if item else None,
                    ES=item['ES'] if item else None,NOMp=item['NOM p-val'] if item else None,
                    tag=item['Tag %'] if item else None,category=cat,skip_reason='' if item else reason))
        print(f'Phase 6 {pos}/50 {drug}', flush=True)
    frame = pd.DataFrame(rows)
    target = out/'tissue_results.csv.gz'
    tmp = target.with_suffix('.tmp')
    frame.to_csv(tmp,index=False,compression='gzip')
    os.replace(tmp,target)
    ref = pd.read_csv(REPO/'verification/results/phase6_comparison_labels.tsv.gz', sep='\t',keep_default_na=False)
    keys = ['drug','lineage','tissue']
    both = frame.merge(ref,on=keys,suffixes=('_run','_release'),how='outer',indicator=True,validate='one_to_one')
    assert (both._merge == 'both').all() and len(both) == 39996
    delta = both[(both.category_run != both.category_release)|(both.skip_reason_run != both.skip_reason_release)]
    delta.to_csv(out/'differences_from_release.csv',index=False)
    scores = pd.read_csv(REPO/'Phase6_NoSampleFilter/results/tested_tissue_scores.tsv.gz',sep='\t')
    observed = frame.loc[frame.category!='D',keys+['NES','FDR']].merge(
        scores, on=keys, suffixes=('_run','_release'), how='outer', indicator=True, validate='one_to_one')
    matched = observed.loc[observed._merge=='both']
    max_difference = {col:float((matched[col+'_run']-matched[col+'_release']).abs().max()) for col in ['NES','FDR']}
    observed.loc[(observed._merge!='both')|
                 ~np.isclose(observed.NES_run, observed.NES_release, atol=1e-12,rtol=1e-12)|
                 ~np.isclose(observed.FDR_run, observed.FDR_release, atol=1e-12,rtol=1e-12)].to_csv(
                     out/'score_differences_from_release.csv',index=False)
    save_json(out/'comparison.json',dict(comparisons=len(both),categories=frame.category.value_counts().to_dict(),
              matched_scores=len(matched),unmatched_scores=len(observed)-len(matched),score_max_absolute_difference=max_difference,
              compared_with='released GSEApy 1.3.2 replay',changed_category_or_skip_reason=len(delta)))


def sampling(out):
    records = pd.read_parquet(Path(os.environ['P59_PREPARED'])/'signatures.parquet')
    keys = ['drug','lineage','tissue']
    current = pd.read_csv(out.parent/'phase6/tissue_results.csv.gz')[keys+['category']]
    ref = pd.read_csv(REPO/'verification/results/phase6_comparison_labels.tsv.gz', sep='\t')[keys+['category']]
    changes = pd.read_csv(REPO/'Phase6_NoSampleFilter/results/replay_category_changes.csv')
    earlier = ref.merge(changes, on=keys, how='left', validate='one_to_one')
    changed = earlier.archived_category.notna()
    assert (earlier.loc[changed,'category']==earlier.loc[changed,'replay_category']).all()
    earlier.loc[changed,'category'] = earlier.loc[changed,'archived_category']

    assignments = {}
    for name, rows in [('rerun',current),('historical',earlier)]:
        best = rows.groupby(['drug','lineage']).category.min()
        linked = records.join(best,on=['cmap_name','cell_lineage'])
        # The notebook assigns unmapped 'unknown' lineages to D.
        linked['category'] = linked.category.fillna('D')
        assignments[name] = linked.set_index('sig_id').category.sort_index()
        counts = linked.groupby(['category','donor_sex']).agg(signatures=('sig_id','size'),lines=('cell_iname','nunique'))
        counts.to_csv(out/(name+'_sampling.csv'))
        by_line = linked[linked.category=='A'].groupby(['cmap_name','cell_lineage','cell_iname','donor_sex']).size()
        stats = []
        for (drug,lin),g in by_line.groupby(level=[0,1]):
            n=int(g.sum()); f=g[g.index.get_level_values('donor_sex')=='F']; m=g[g.index.get_level_values('donor_sex')=='M']
            stats.append(dict(drug=drug,cell_lineage=lin,total_sigs=n,cell_lines=len(g),biggest_line=int(g.max()),
                         smallest_line=int(g.min()),effective_lines=n*n/float((g*g).sum()),pct_from_biggest=100*g.max()/n,
                         sigs_female=int(f.sum()), sigs_male=int(m.sum()), lines_female=len(f), lines_male=len(m)))
        stats=pd.DataFrame(stats)
        stats['inflation'] = (stats.cell_lines/stats.effective_lines.round(2)).round(2)
        stats['sigs_total'] = stats.total_sigs
        stats['both_sexes_present'] = stats[['lines_female','lines_male']].min(axis=1) >= 1
        stats['both_sexes_2plus'] = stats[['lines_female','lines_male']].min(axis=1) >= 2
        stats.to_csv(out/(name+'_category_A_pairs.csv'),index=False)
        if name=='historical':
            refdir = REPO/'Phase7_SexBalance_Audit/results'
            published_counts = pd.read_csv(refdir/'signatures_by_category_and_sex.csv').set_index('category')
            published_lines = pd.read_csv(refdir/'cell_lines_by_category_and_sex.csv').set_index('category')
            for (category,sex), row in counts.iterrows():
                suffix = 'female' if sex=='F' else 'male'
                assert row.signatures == published_counts.loc[category,'signatures_'+suffix]
                assert row.lines == published_lines.loc[category,'lines_'+suffix]
            for f in ['pseudoreplication_category_A.csv','category_A_by_drug_lineage_sex.csv']:
                orig=pd.read_csv(refdir/f)
                cmp=stats.merge(orig,on=['drug','cell_lineage'],validate='one_to_one',suffixes=('_run','_archive'))
                assert len(cmp)==len(orig)==len(stats)
                for col in set(orig)&set(stats)-{'drug','cell_lineage'}:
                    tol=.00501 if col in {'effective_lines','inflation'} else .0501 if col=='pct_from_biggest' else 0
                    assert np.allclose(cmp[col+'_run'],cmp[col+'_archive'],rtol=0,atol=tol),col
        save_json(out/(name+'_counts.json'),dict(signatures=len(linked),category_A_pairs=len(stats),
                  category_A_signatures=int(linked.category.eq('A').sum()),female_A=int(counts.loc[('A','F'),'signatures']),
                  male_A=int(counts.loc[('A','M'),'signatures'])))

    original = assignments['historical'].rename('sig_cat')
    checksum = hashlib.sha256(original.to_csv(sep='\t').encode()).hexdigest()
    expected = json.loads((REPO/'Phase7_SexBalance_Audit/results/archived_signature_labels.json').read_text())['sha256']
    if checksum != expected:
        raise ValueError('The reconstructed historical signature assignments differ from the October export')
    current = assignments['rerun']
    differences = original[current!=original].to_frame().join(current.rename('replay_category'))
    differences = differences.join(records.set_index('sig_id')[['cmap_name','cell_lineage','donor_sex']])
    changes = differences.groupby(['cmap_name','cell_lineage','donor_sex','sig_cat','replay_category']).size()
    changes.rename('signatures').to_csv(out/'changed_signature_counts.csv')
    save_json(out/'assignment_comparison.json', dict(historical_assignments=len(original),
        historical_sha256=checksum, changed_replay_assignments=len(differences)))


def faers(out):
    folder = REPO/'Phase8_SideEffects_FAERS'
    sheet = pd.read_csv(folder/'per_drug_events.csv').set_index('drug').sort_index()
    pattern = r'^\s*(\S+)\s+GROUP \d \([^)]+\)\s+([AB-])\s+(\d+)\s+(\d+)\s+(\d+)'
    matches = re.findall(pattern,(folder/'archived_cell_20.txt').read_text(),re.MULTILINE)
    rows = pd.DataFrame(matches, columns=['drug','phase6_best_enrichment','female_and_male_reports',
                                         'reaction_terms_tested','reaction_terms_flagged'])
    rows['phase6_best_enrichment'] = rows.phase6_best_enrichment.replace('-','none')
    for col in ['female_and_male_reports','reaction_terms_tested','reaction_terms_flagged']:
        rows[col] = pd.to_numeric(rows[col])
    pd.testing.assert_frame_equal(rows.set_index('drug').sort_index()[sheet.columns],sheet)

    ref=pd.read_csv(REPO/'Phase6_NoSampleFilter/results/signature_ABCD_per_drug.csv').set_index('drug').loc[rows.drug]
    groups=np.select([ref.A.to_numpy()>0,ref.B.to_numpy()>0],['A','B'],default='none')
    assert (groups == rows.phase6_best_enrichment).all()
    percent=(100*rows.reaction_terms_flagged/rows.reaction_terms_tested).round(1)
    hit=groups!='none'
    test=mannwhitneyu(percent[hit],percent[~hit],alternative='two-sided',method='exact')
    rho=spearmanr(np.log10(rows.female_and_male_reports),percent)
    save_json(out/'comparison.json',dict(source='archived October 2026 drug-level counts',drugs=len(rows),
         enriched_drugs=int(hit.sum()),comparison_drugs=int((~hit).sum()),
         total_tested_terms=int(rows.reaction_terms_tested.sum()),enriched_median_percent=float(percent[hit].median()),
         comparison_median_percent=float(percent[~hit].median()),U=float(test.statistic),p=float(test.pvalue),
         spearman_rho=float(rho.statistic)))


if __name__=='__main__':
    n=int(sys.argv[1])
    dest=Path(sys.argv[2])/f'phase{n}'
    dest.mkdir(parents=True, exist_ok=True)
    {1:pilot,6:unfiltered,7:sampling,8:faers}[n](dest)
