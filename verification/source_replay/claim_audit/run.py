"""Check published Part A claims from line-level measurements, source labels and DNA calls."""
import argparse
import csv
import hashlib
import json
import statistics
from collections import Counter, defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
RUN = None
PA = None
SUPP = None
EXPECTED_PATTERNS = ['Y only', 'XIST only', 'neither', 'both']
MAP = {(True, False): 'Y only', (False, True): 'XIST only', (False, False): 'neither', (True, True): 'both'}
SEX = {'Male': 'Y only', 'Female': 'XIST only'}



def load(path):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def unique(records, col):
    result = {r[col]: r for r in records if r[col]}
    assert len(result) == sum(bool(r[col]) for r in records), (col, len(result), len(records))
    return result


def save(path, records):
    with path.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, list(records[0]), delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(records)


def frac(matches, denominator):
    return f'{matches}/{denominator}'


class Audit:
    def __init__(self, out):
        self.out = out
        out.mkdir(parents=True, exist_ok=True)
        self.claims = []
        self.files = {}

    def source(self, path):
        self.files[path.relative_to(RUN).as_posix()] = hashlib.sha256(path.read_bytes()).hexdigest()
        return load(path)

    def check(self, name, published, actual, source):
        added = (name in {'LINCS donor-labelled lines with both consensus flags and atypical pattern',
                          'All resolved state names lacking one or both consensus flags'} or
                 (name.startswith('GDSC2/PRISM ') and 'single-marker low / measured' in name))
        basis = 'independent fresh cohort check' if added else 'original manuscript claim'
        self.claims.append(dict(claim=name, reference=str(published), calculated=str(actual), basis=basis,
                                verdict='matches' if str(published) == str(actual) else 'differs', source=source))

    def write(self):
        save(self.out/'claims.tsv', self.claims)
        (self.out/'source_sha256.json').write_text(json.dumps(self.files, indent=2, sort_keys=True)+'\n')


def enrich(lines, cuts):
    for r in lines:
        y, x = float(r['y_score']), float(r['XIST'])
        r['pattern'] = MAP[y >= cuts['y_score'], x >= cuts['XIST']]
    return lines


def bulk_refs(a):
    root = RUN/'bulk_RNA/results/C'
    bounds = a.source(root/'C_expression_cutoffs.tsv')
    cuts = {(r['expr_source'], r['variable']): float(r['cutoff_log2']) for r in bounds}
    records = a.source(root/'C_cellline_sex_labels.tsv')
    matches = unique(records, 'RRID')
    for r in matches.values():
        s = r['expressed_state_source']
        if s not in ('CCLE','CMP'):
            continue
        r['pattern'] = MAP[float(r[s.lower()+'_yscore']) >= cuts[s,'yscore'],
                           float(r[s.lower()+'_xist']) >= cuts[s,'xist']]
    return matches


def compare(a, which, lines, other, state_key):
    paired = []
    for r in lines:
        p = other.get(r['cellosaurus_id'])
        if not p or 'pattern' not in p:
            continue
        assay = p.get('expressed_state_source', 'EmeraldBay')
        quality = r['confidence'] in {'high','medium'} and p.get('confidence') in {'high','medium'}
        paired.append(dict(accession=r['cellosaurus_id'], name=r['cell_name'], sample=r.get('sample',''),
                           first_pattern=r['pattern'], second_pattern=p['pattern'],
                           pattern_agrees=r['pattern']==p['pattern'], first_state=r['state'], second_state=p[state_key],
                           state_agrees=r['state']==p[state_key], second_source=assay,
                           first_donor=r['annotation'], second_donor=p.get('annotated_sex') or p.get('annotation'),
                           first_confidence=r['confidence'], second_confidence=p.get('confidence',''),
                           both_at_least_medium=quality if which=='tahoe_emerald' else ''))
    if which=='tahoe_bulk':
        save(a.out/'tahoe_bulk_disagreements.tsv',[r for r in paired if not r['pattern_agrees'] or not r['state_agrees']])
    return paired


def check_comparison(a, group, expected, records, field):
    a.check(group, expected, frac(sum(r[field] for r in records),len(records)), field)


def dna(a, rows, refs, label, expected):
    pairs=[]
    for r in rows:
        assay = refs[r['cellosaurus_id']]['amelogenin_class']
        if assay not in {'X,Y', 'X only'}:
            continue
        expressed = r['pattern'] in {'Y only','both'}
        pairs.append(dict(accession=r['cellosaurus_id'],cell_name=r['cell_name'],source_sample=r.get('sample',''),
                          rna_pattern=r['pattern'],rna_Y_high=expressed,dna_amelogenin=assay,
                          match=expressed == (assay == 'X,Y')))
    grouped=defaultdict(set)
    for p in pairs:
        grouped[p['accession']].add(p['match'])
    assert all(len(v)==1 for v in grouped.values()), 'Within-line Y calls differ by experiment'
    a.check(label+' marker vs DNA (unique accessions)',expected,
            frac(sum(v=={True} for v in grouped.values()),len(grouped)), 'raw Cellosaurus Amelogenin class vs expression cutoff')
    return pairs


def do_external(a):
    src = PA/'vehicle_tables'
    f1cuts=a.source(src/'F1_cutoffs.tsv')
    cuts={k:{r['score']:float(r['cutoff']) for r in f1cuts if r['dataset']==k} for k in ('EmeraldBay','sci-Plex (all samples)')}
    e=enrich(a.source(src/'F1_emerald_states.tsv'),cuts['EmeraldBay'])
    s=enrich(a.source(src/'F1_sciplex_states.tsv'),cuts['sci-Plex (all samples)'])
    b=bulk_refs(a)
    strref=unique(a.source(PA/'str/amelogenin_by_record.tsv'),'accession')
    a.check('EmeraldBay endpoint vehicle cells',100253,sum(int(r['n_vehicle_cells']) for r in e),'F1')
    a.check('EmeraldBay line states XY, Y-lost, XIST-lost, XX-like, ambiguous','17,16,12,6,1',
            ','.join(str(sum(r['state']==k for r in e)) for k in ['XY','Y-lost','XIST-lost','XX-like','ambiguous']),'F1')
    a.check('EmeraldBay confidence high, medium, low','42,4,6',
            ','.join(str(sum(r['confidence']==k for r in e)) for k in ['high','medium','low']),'F1')
    for name,sub,expected in [('all',e,'28/50'),('medium and high',[r for r in e if r['confidence']!='low'],'24/44')]:
        lab=[r for r in sub if r['annotation'] in SEX]
        a.check('EmeraldBay atypical donor pattern '+name,expected,frac(sum(r['pattern']!=SEX[r['annotation']] for r in lab),len(lab)),'F1 scores and cutoffs')
    a.check('sci-Plex vehicle cells by experiment','2631,17578,3361',
            ','.join(str(sum(int(r['n_vehicle_cells']) for r in s if r['sample']==k)) for k in ('sciPlex2','sciPlex3','sciPlex4')),'F1')
    sp=compare(a,'sciplex_bulk',s,b,'expressed_state')
    check_comparison(a,'sci-Plex against bulk, patterns','5/6',sp,'pattern_agrees')
    check_comparison(a,'sci-Plex against bulk, state names','6/6',sp,'state_agrees')
    dna(a,e,strref,'emerald','47/52')
    dna(a,s,strref,'sciplex','3/3')

    bulk_str=a.source(PA/'str_fresh_bulk.tsv')
    for name,expected in [('CCLE','714/734'),('CMP','1007/1031')]:
        tested=[r for r in bulk_str if r['source']==name and r['comparable']=='True']
        a.check('Figure 3 '+name+' marker vs DNA',expected,
                frac(sum((r['rna_call']=='Y expressed')==(r['amelogenin']=='X,Y') for r in tested),len(tested)),
                'raw-derived bulk RNA and parsed DNA')
    both_bulk=[r for r in b.values() if r['cmp_state'] not in ('not_measured','') and r['ccle_state'] not in ('not_measured','')]
    a.check('S1 GDSC2/PRISM state-name agreement across the two bulk sources','729/755',
            frac(sum(r['cmp_state']==r['ccle_state'] for r in both_bulk),len(both_bulk)),'C source states')

    lin=a.source(RUN/'part_b/run/A/results/A/A_cellline_sex_states.tsv')
    both_lin=[r for r in lin if r['state_sanger']!='no_data' and r['state_ccle']!='no_data']
    a.check('S1 LINCS state-name agreement across the two sources','17/18',
            frac(sum(r['state_sanger']==r['state_ccle'] for r in both_lin),len(both_lin)),'A source states')
    a.check('S1 LINCS named lines with CCLE measurements',135,sum(r['state_ccle']!='no_data' for r in lin),'A source states')
    donors=[r for r in lin if r['cmap_donor_sex'] in ('M','F')]
    has_rna=[r for r in donors if r['state']!='no_data']
    called=[r for r in has_rna if r['state']!='ambiguous']
    a.check('LINCS named donor labels identical','141/141',frac(sum(r['cellosaurus_sex']=={'M':'Male','F':'Female'}[r['cmap_donor_sex']] for r in donors),len(donors)),'A')
    a.check('S10 donor-labelled names with external RNA regardless of call',128,len(has_rna),'A')
    a.check('LINCS donor-labelled resolved atypical states','61/128',frac(sum(r['state'] in {'Y-lost','XIST-lost','discordant'} for r in called),len(called)),'A')
    a.check('LINCS measured names, distinct accessions after alias merger','140/139',
            frac(sum(r['state']!='no_data' for r in lin),len({r['cellosaurus_ac'] for r in lin if r['state']!='no_data'})),'A')
    complete=[r for r in called if r['Yhi'] in ('True','False') and r['Xhi'] in ('True','False')]
    atypical=sum(MAP[r['Yhi']=='True',r['Xhi']=='True']!=SEX[{'M':'Male','F':'Female'}[r['cmap_donor_sex']]] for r in complete)
    a.check('LINCS donor-labelled lines with both consensus flags and atypical pattern','61/127',frac(atypical,len(complete)),'A')
    exceptions=[{k:r[k] for k in ['cell_iname','cellosaurus_ac','state','state_ccle','state_sanger','Yhi','Xhi',
                'Yhi_ccle','Xhi_ccle','Yhi_sanger','Xhi_sanger']} for r in lin if r['cell_iname'] in {'HCT116','U266B1'}]
    save(a.out/'lincs_exceptions.tsv',exceptions)
    a.check('All resolved state names lacking one or both consensus flags','U266B1',
            ','.join(sorted(r['cell_iname'] for r in called if r not in complete)),'A')
    return e,b,strref


def do_tahoe(a,root,e,b,d):
    def src(tag):
        return a.source(root/('partA_local_'+tag+'.tsv'))
    cuts={r['score']:float(r['cutoff']) for r in src('cutoffs') if r['level']=='lineplate'}
    t=enrich(src('states'),cuts)
    ti=unique(t,'cellosaurus_id')
    a.check('Figure 2 cutoffs Y, XIST, Yp, Yq','4.02,4.68,3.96,3.58',','.join(f'{cuts[k]:.2f}' for k in ['y_score','XIST','yp_score','yq_score']),'local cutoff table')
    a.check('Vehicle cells before filter',2205786,sum(int(r['n_vehicle_cells_before_filter']) for r in t),'local line counts')
    counts=[]
    for group,expected in [('Male',[18,1,14,0]),('Female',[0,4,12,0]),('Unlabelled',[0,1,0,0]),('All',[18,6,26,0])]:
        group_rows=[r for r in t if (group=='All' or (group=='Unlabelled' and r['annotation'] not in SEX) or r['annotation']==group)]
        values=[sum(r['pattern']==p for r in group_rows) for p in EXPECTED_PATTERNS]
        counts.append(dict(recorded_donor=group,total=len(group_rows),**dict(zip(EXPECTED_PATTERNS,values))))
        a.check('Table 1 counts '+group,','.join(map(str,expected)),','.join(map(str,values)),'line scores rethresholded')
    save(a.out/'table1.tsv',counts)
    lab=[r for r in t if r['annotation'] in SEX]
    a.check('Abstract donor atypical','27/49',frac(sum(r['pattern'] != SEX[r['annotation']] for r in lab),len(lab)),'line scores rethresholded')
    altcuts={r['score']:float(r['cutoff']) for r in src('cutoffs_linelevel') if r['level']=='line'}
    sensitivity=enrich(src('states_linelevel'),altcuts)
    sal=[r for r in sensitivity if r['annotation'] in SEX]
    a.check('Alternative line-level cutoffs, atypical labels','27/49',frac(sum(r['pattern']!=SEX[r['annotation']] for r in sal),len(sal)),'full line calls, line-level threshold refit')
    pie=unique(src('pie_mapping'),'cellosaurus_id')
    pl=[r for r in t if pie[r['cellosaurus_id']]['pie_sex'] in SEX]
    a.check('PIE recorded label, atypical expression','27/49',frac(sum(r['pattern']!=SEX[pie[r['cellosaurus_id']]['pie_sex']] for r in pl),len(pl)),'PIE label joined to score-derived pattern')
    sb=compare(a,'line_level_cutoffs_bulk',sensitivity,b,'expressed_state')
    check_comparison(a,'S9 external bulk state names under line-level cutoffs','46/49',sb,'state_agrees')
    male=[r for r in t if r['annotation']=='Male']
    arms=Counter((float(r['yp_score'])>=cuts['yp_score'])+(float(r['yq_score'])>=cuts['yq_score']) for r in male)
    a.check('Figure 2 males both, neither, one arm high','17,12,4',','.join(str(arms[k]) for k in (2,0,1)),'line arm scores')
    s6={'CFPAC-1':'5.37,2.12,2.93','NCI-H2030':'7.28,1.04,2.60','COLO 205':'1.46,4.40,3.67','HCT15':'3.57,5.64,5.12'}
    for name,reading in s6.items():
        row=next(r for r in t if r['cell_name']==name)
        a.check('Table S6 Yp, Yq, composite for '+name,reading,
                ','.join(f'{float(row[k]):.2f}' for k in ('yp_score','yq_score','y_score')),'per-line arm scores')
    line=next(r for r in t if r['cell_name']=='NCI-H2030')
    before,after=int(line['n_vehicle_cells_before_filter']),int(line['n_vehicle_cells'])
    a.check('NCI-H2030 marker-rich filter removed','22186/46472',frac(before-after,before),'line counts pre/post filter')
    p=src('lineplate_scores')
    a.check('Table 2 plate agrees with pooled call','685/700',frac(sum(r['plate_state']==ti[r['cell_line_id']]['state'] for r in p),len(p)),'line by plate state joined to pooled state')
    a.check('All plates agree in line','48/50',frac(sum(r['n_plates']==r['n_plates_agree'] for r in t),len(t)),'line table')
    tb=compare(a,'tahoe_bulk',t,b,'expressed_state')
    check_comparison(a,'Table 2 Tahoe vs bulk marker patterns','47/49',tb,'pattern_agrees')
    check_comparison(a,'Note S9 Tahoe vs bulk state names','46/49',tb,'state_agrees')
    te=compare(a,'tahoe_emerald',t,unique(e,'cellosaurus_id'),'state')
    check_comparison(a,'Figure 3 Tahoe vs EmeraldBay patterns, all','42/47',te,'pattern_agrees')
    check_comparison(a,'Figure 3 Tahoe vs EmeraldBay patterns, medium/high in both','40/41',[r for r in te if r['both_at_least_medium']],'pattern_agrees')
    ds=dna(a,t,d,'tahoe','48/50')
    male_neither=[r for r in ds if ti[r['accession']]['annotation']=='Male' and r['rna_pattern']=='neither']
    a.check('Male-labelled Tahoe lines with neither marker lacking AMELY','13/14',frac(sum(r['dna_amelogenin']=='X only' for r in male_neither),len(male_neither)),'DNA profiles vs rethresholded scores')
    amb=[r for r in src('ambient') if r['unit_type']=='plate' and r['cells']=='marker < 3 (foreign cells removed)']
    for marker,term,expected,places in [('Y','both arms low','0.012-0.033',3),('XIST','XIST-negative','0.0013-0.0040',4)]:
        aa=[r for r in amb if term in r['line_set']]
        a.check('Figure S1 '+marker+' plate bounds per 1000 UMI',expected,f'{min(float(r["per_1000"]) for r in aa):.{places}f}-{max(float(r["per_1000"]) for r in aa):.{places}f}','per-plate ambient table with foreign cells excluded')
        a.check('Figure S1 '+marker+' plate count',14,len(aa),'ambient')
    fo=[r for r in src('foreign_cells') if r['arm_state']=='Y lost']
    share=[float(r['share_foreign_Y'])*100 for r in fo]
    a.check('Figure S1 fully arm-low lines',28,len(fo),'foreign cells')
    a.check('Figure S1 range in percentage of marker-rich cells','0.05-0.33',f'{min(share):.2f}-{max(share):.2f}','foreign cells')
    a.check('Figure S1 median observed / expected, fully arm-low',56,round(statistics.median(float(r['obs_over_exp_Y']) for r in fo)),'foreign cells')
def do_mosaic(a,root):
    def src(tag):
        return a.source(root/('partA_local_'+tag+'.tsv'))
    depths=[r for r in src('depth_curve') if r['cell_name']=='A549' and r['marker']=='Y']
    a.check('Figure S2 A549 observed Y-zero rate above homogeneous Poisson in each depth bin',
            '10/10',frac(sum(float(r['zero_frac'])>float(r['homogeneous_poisson_zero_frac']) for r in depths),len(depths)),
            'depth-binned observed and modeled zero fractions')
    abc=src('mosaic_abcde')
    a.check('In-scope line-marker pairs',22,sum(r['scope']=='in scope' for r in abc),'mosaic fit rows')
    for name,expected in [('A549',60135),('NCI-H460',126343)]:
        r=next(r for r in abc if r['cell_name']==name and r['marker']=='Y')
        a.check('Figure S2 Y cells, '+name,expected,int(r['n_cells']),'local mosaic model counts')


def do_s1(a,path):
    tab=a.source(path)
    a.check('Table S1 number of distinct rows',1694,len({(r['cellosaurus_id'],r['source']) for r in tab}),'S1')
    a.check('Table S1 number of distinct accessions',1123,len({r['cellosaurus_id'] for r in tab}),'S1')
    a.check('Table S1 rows without RNA',21,sum(r['state']=='not_measured' for r in tab),'S1')
    for source,count in [('Tahoe-100M',50),('EmeraldBay',52),('sciPlex2',1),('sciPlex3',3),('sciPlex4',2),('LINCS2020',139),('GDSC2',969),('PRISM',478)]:
        a.check('Table S1 '+source,count,sum(r['source']==source for r in tab),'S1')


def do_a3(a,path):
    tested=[r for r in a.source(path) if r['used']=='True']
    misses=[f'{r["line"]}/{r["dataset"]}/{r["estimator"]}' for r in tested if r['line']=='A549' and float(r['abs_diff'])>0.15]
    a.check('Figure S2 A549 pairs exceeding tolerance','A549/sciPlex3/a3000',','.join(misses),'A3 pairs with usable estimates')

    tex=SUPP
    a.files['source_text/Table_S8.tex']=hashlib.sha256(tex.read_bytes()).hexdigest()
    block=tex.read_text().split('Table S8. Depth-matched shares, Tahoe against the other',1)[1]
    block=block.split(r'\endlastfoot',1)[1].split(r'\end{longtable}',1)[0]
    published={}
    for item in block.strip().split(r'\\'):
        if not item.strip(): continue
        parts=[' '.join(t.split()) for t in item.split('&')]
        line_name,marker,atlas=parts[0].replace(' (context)','').split(', ')
        estimator={'(a)':'a3000','(b)':'b2000','(c)':'c'}[parts[1][:3]]
        published[line_name,atlas,estimator]=(parts[2],*parts[3:7])
    assert len(published)==len(tested)
    for r in tested:
        k=(r['line'],r['dataset'],r['estimator'])
        rendered=(r['thinned_side'], *(f'{float(r[col]):.3f}' for col in ['share_tahoe','share_other','abs_diff']),
                  'yes' if r['agree']=='True' else 'no')
        a.check('Table S8 '+ '/'.join(k),'|'.join(published[k]),'|'.join(rendered),'A3 measured pair')


def main():
    global RUN,PA,SUPP
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run-root',required=True,type=Path,help='Raw reproduction directory with part_a, part_b, bulk_RNA and manuscript_review')
    parser.add_argument('--supplement',type=Path,help='Manuscript supplement TeX containing Table S8')
    parser.add_argument('--output',type=Path,default=HERE/'results')
    args=parser.parse_args()
    RUN=args.run_root.resolve()
    PA=RUN/'part_a'
    SUPP=args.supplement or RUN/'manuscript_review/baseline/manuscript/supplement.tex'
    audit=Audit(args.output)
    e,b,d=do_external(audit)
    dossier={}
    root=PA/'results'
    cuts={r['score']:float(r['cutoff']) for r in audit.source(root/'partA_local_cutoffs.tsv') if r['level']=='lineplate'}
    t=enrich(audit.source(root/'partA_local_states.tsv'),cuts)
    assert all((float(r['y_score'])>=cuts['y_score'])==(r['y_high']=='True')
           and (float(r['XIST'])>=cuts['XIST'])==(r['xist_high']=='True') for r in t)
    tb=compare(audit,'tahoe_bulk',t,b,'expressed_state')
    check_comparison(audit,'Table 2 Tahoe vs bulk marker patterns','47/49',tb,'pattern_agrees')
    check_comparison(audit,'Note S9 Tahoe vs bulk state names','46/49',tb,'state_agrees')
    dossier['Tahoe line states and four-way cutoffs']='complete'
    crows=audit.source(RUN/'bulk_RNA/results/C/C_cellline_sex_labels.tsv')
    ccuts={(r['expr_source'],r['variable']):float(r['cutoff_log2']) for r in audit.source(RUN/'bulk_RNA/results/C/C_expression_cutoffs.tsv')}
    for sex,neligible,nmeasured,n_neither,nmarker,marker in [('Male',596,579,256,272,'yscore'),('Female',302,300,122,125,'xist')]:
        rows=[r for r in crows if r['annotated_sex']==sex and r['lineage'] and r['sex_specific_lineage']=='False']
        measured=[r for r in rows if r['expressed_state_source'] in {'CCLE','CMP'}]
        def below(r,m):
            s=r['expressed_state_source']
            return float(r[s.lower()+'_'+m]) < ccuts[s,m]
        joint=sum(below(r,'yscore') and below(r,'xist') for r in measured)
        markers=sum(below(r,marker) for r in measured)
        audit.check('GDSC2/PRISM '+sex+' joint low / all eligible',frac(n_neither,neligible),frac(joint,len(rows)),'C selected source, raw cutoffs')
        audit.check('GDSC2/PRISM '+sex+' single-marker low / measured',frac(nmarker,nmeasured),frac(markers,len(measured)),'C selected source, raw cutoffs')
    if all((root/('partA_local_'+tag+'.tsv')).is_file() for tag in ['foreign_cells','pie_mapping','states_linelevel']):
        do_tahoe(audit,root,e,b,d)
        dossier['Tahoe plate and cutoff sensitivity, Figures 2, 3 and S1']='complete'
    else:
        dossier['Tahoe plate and cutoff sensitivity, Figures 2, 3 and S1']='awaiting downstream tables'
    if all((root/('partA_local_'+tag+'.tsv')).is_file() for tag in ['mosaic_abcde','depth_curve']):
        do_mosaic(audit,root)
        dossier['Figure S2A']='complete'
    else:
        dossier['Figure S2A']='awaiting fitted local model'
    s1=RUN/'manuscript_review/supplementary_table_S1.tsv'
    if s1.is_file() and (s1.parent/'table_s1_comparison.json').is_file():
        do_s1(audit,s1)
        dossier['Reassembled Table S1']='complete'
    else:
        dossier['Reassembled Table S1']='awaiting assembly'
    a3=PA/'external_checks/A3_pairs.tsv'
    if a3.is_file() and (a3.parent/'A3_summary.json').is_file():
        do_a3(audit,a3)
        dossier['Figure S2B and all Table S8 pairs']='complete'
    else:
        dossier['Figure S2B and all Table S8 pairs']='awaiting depth matching'
    audit.claims=list({r['claim']:r for r in audit.claims}.values())
    audit.write()
    (args.output/'coverage.json').write_text(json.dumps(dossier,indent=2)+'\n')
    print(dossier)


if __name__=='__main__':
    main()
