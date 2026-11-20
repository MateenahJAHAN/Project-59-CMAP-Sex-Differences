"""Build token -> chromosome table from the Ensembl 110 GTF (gene lines only).

Output: partA_local/gene_chrom.tsv with token_id, gene_symbol, ensembl_id, chrom, chrom_class
chrom_class: autosome, X, Y, MT, other (scaffold), unannotated (ID not in the GTF).
"""
import gzip, sys, time
import pandas as pd

OUT = __import__('os').environ.get('P59_PARTA_OUT', '$P59_DATA/tahoe/partA_local')
GTF = OUT + '/downloads/ensembl/Homo_sapiens.GRCh38.110.gtf.gz'
t0 = time.time()
rows = {}
with gzip.open(GTF, 'rt') as fh:
    for line in fh:
        if line.startswith('#'):
            continue
        f = line.split('\t', 9)
        if f[2] != 'gene':
            continue
        attr = f[8]
        i = attr.index('gene_id "') + 9
        gid = attr[i:attr.index('"', i)]
        assert gid not in rows, gid
        rows[gid] = f[0]
gtf = pd.Series(rows, name='chrom')
genes = pd.read_csv('$P59_DATA/tahoe/dmso/genes.tsv', sep='\t')
assert genes.token_id.is_unique and genes.ensembl_id.is_unique
genes['chrom'] = genes.ensembl_id.map(gtf)
auto = {str(i) for i in range(1, 23)}

def cls(c):
    if pd.isna(c):
        return 'unannotated'
    if c in auto:
        return 'autosome'
    if c in ('X', 'Y', 'MT'):
        return c
    return 'other'

genes['chrom_class'] = genes.chrom.map(cls)
genes.to_csv(OUT + '/gene_chrom.tsv', sep='\t', index=False)
print('GTF genes', len(gtf), 'Tahoe genes', len(genes))
print(genes.chrom_class.value_counts().to_string())
mk = ['XIST', 'TSIX', 'RPS4Y1', 'DDX3Y', 'KDM5D', 'UTY', 'EIF1AY', 'ZFY', 'USP9Y', 'NLGN4Y', 'TXLNGY', 'PRKY',
      'RPS4X', 'DDX3X', 'KDM5C', 'KDM6A', 'EIF1AX', 'ZFX', 'USP9X', 'NLGN4X', 'ACTB', 'GAPDH']
print(genes[genes.gene_symbol.isin(mk)][['gene_symbol', 'token_id', 'chrom']].to_string(index=False))
print('seconds', round(time.time() - t0, 1))
