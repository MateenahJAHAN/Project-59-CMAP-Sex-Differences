"""Canonical whole-cohort checks. Only complete-table digests leave the local analysis."""
import hashlib

Y_GENES=['RPS4Y1','DDX3Y','KDM5D','UTY','EIF1AY','ZFY','USP9Y','NLGN4Y']
MARKERS=['XIST',*Y_GENES]


def digest_rows(rows):
    digest=hashlib.sha256()
    for r in rows:
        digest.update(('\t'.join(str(v) for v in r)+'\n').encode('utf-8'))
    return digest.hexdigest()


def cohort_digest(frame):
    cols=['cell_iname','state','sample_id','row_index','n_available']
    return digest_rows(frame.sort_values('row_index')[cols].itertuples(index=False,name=None))


def median_digest(frame):
    return digest_rows((name,*(f'{float(r[k]):.4f}' for k in [*MARKERS,'L1000_Y_score']))
                       for name,r in frame.sort_index().iterrows())
