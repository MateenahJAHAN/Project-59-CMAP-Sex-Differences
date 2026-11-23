"""Encode the full Tahoe observation table without materializing its string columns in memory."""
import json
from pathlib import Path
import time

import numpy as np
import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.parquet as pq

import tp_common as spec


def main():
    start = time.time()
    index = Path(spec.DATA) / "index"
    source = pq.ParquetFile(spec.META["obs_metadata"])
    count = source.metadata.num_rows
    if count != 100_648_790:
        raise ValueError(f"Unexpected observation count {count}")
    columns = ["BARCODE_SUB_LIB_ID", "plate", "sample", "sublibrary", "cell_line", "drug", "pass_filter",
               "phase", "tscp_count", "gene_count", "pcnt_mito"]
    phase_counts = {phase: 0 for phase in spec.PHASES}
    filters = set()
    plates = {}
    keys = np.empty(count, dtype=np.int64)
    offset = full_count = 0
    temp = index / "obs_compact.parquet.part"
    writer = None
    try:
        for batch in source.iter_batches(batch_size=262_144, columns=columns):
            data = pa.Table.from_batches([batch])
            size = len(data)
            end = offset + size
            key = spec.encode_keys(data["BARCODE_SUB_LIB_ID"])
            keys[offset:end] = key
            plate = spec.plate_num(data["plate"]).astype(np.int8)
            for value in pc.unique(pa.array(plate)).to_pylist():
                plates[value] = None
            sub = data["sublibrary"].combine_chunks()
            assert pc.all(pc.starts_with(sub, "lib_")).as_py()
            sublib = pc.cast(pc.utf8_slice_codeunits(sub, 4), pa.int32()).to_numpy()
            assert np.array_equal(key // 1_000_000_000, sublib)
            pf = pc.cast(data["pass_filter"], pa.string()).combine_chunks()
            filters.update(pc.unique(pf).to_pylist())
            full = pc.equal(pf, "full")
            full_count += pc.sum(full).as_py()
            phase = pc.index_in(pc.cast(data["phase"], pa.string()), value_set=pa.array(spec.PHASES))
            if phase.null_count:
                raise ValueError("Unexpected cell-cycle phase")
            phase = pc.cast(phase, pa.int8())
            for idx, name in enumerate(spec.PHASES):
                phase_counts[name] += pc.sum(pc.equal(phase, idx)).as_py()
            out = pa.table({
                "obs_row": np.arange(offset, end, dtype=np.int64), "key": key, "plate": plate,
                "sample": spec.sample_num(data["sample"]).astype(np.int16), "sublibrary": sublib,
                "cell_line": pc.dictionary_encode(pc.cast(data["cell_line"], pa.string())),
                "drug": pc.dictionary_encode(pc.cast(data["drug"], pa.string())), "pass_full": full, "phase": phase,
                "tscp_count": pc.cast(data["tscp_count"], pa.int32()), "gene_count": pc.cast(data["gene_count"], pa.int32()),
                "pcnt_mito": data["pcnt_mito"],
            })
            if writer is None:
                writer = pq.ParquetWriter(temp, out.schema, compression="zstd")
            writer.write_table(out)
            offset = end
    finally:
        if writer is not None:
            writer.close()
    assert offset == count and filters == set(spec.PASS_FILTERS)
    keys.sort()
    assert np.all(keys[1:] != keys[:-1]), "duplicate barcode"
    temp.replace(index / "obs_compact.parquet")
    summary = dict(rows=count, pass_full=full_count, distinct_keys=count, phase_counts=phase_counts,
                   plates_in_obs_order=list(plates), seconds=round(time.time() - start, 1), utc=spec.utcnow())
    (index / "obs_compact_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(summary, flush=True)


if __name__ == "__main__":
    main()
