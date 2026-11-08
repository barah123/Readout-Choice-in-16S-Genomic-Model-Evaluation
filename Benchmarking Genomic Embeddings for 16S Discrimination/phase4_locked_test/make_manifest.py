"""Phase 5 unsealing checklist item 3: archive the accession list and split manifest hash.

Produces a deterministic record of exactly which sequences the locked test was built from,
so the result can be reproduced or audited later.
"""
import glob
import hashlib
import json
import os

rows = []
for f in sorted(glob.glob("raw_species_test/*.json")):
    for r in json.load(open(f)):
        rows.append({
            "block_id": r["block_id"], "species_query": r["species_query"],
            "is_target": r["is_target"], "accession": r["accession"],
            "isolate_group": r["isolate_group"], "slen": r["slen"],
            "seq_sha256": hashlib.sha256(r["sequence"].encode()).hexdigest(),
        })
rows.sort(key=lambda x: (x["block_id"], x["species_query"], x["accession"]))

with open("accession_manifest.tsv", "w") as fh:
    cols = ["block_id", "species_query", "is_target", "accession", "isolate_group",
            "slen", "seq_sha256"]
    fh.write("\t".join(cols) + "\n")
    for r in rows:
        fh.write("\t".join(str(r[c]) for c in cols) + "\n")

man_hash = hashlib.sha256(open("accession_manifest.tsv", "rb").read()).hexdigest()

# split manifest: the exact isolate->fold assignment is a deterministic function of
# (test_blocks.json, SEED=0, N_SPLITS=5), so hashing the block table pins the split.
blocks_hash = hashlib.sha256(open("test_blocks.json", "rb").read()).hexdigest()

meta = {
    "n_records": len(rows),
    "n_accessions": len({r["accession"] for r in rows}),
    "n_isolate_groups": len({r["isolate_group"] for r in rows}),
    "n_blocks": len({r["block_id"] for r in rows}),
    "accession_manifest_sha256": man_hash,
    "test_blocks_json_sha256": blocks_hash,
    "split_definition": "StratifiedKFold(n_splits=5, shuffle=True, random_state=0) over "
                        "isolate_group per block; deterministic given test_blocks.json",
}
json.dump(meta, open("manifest_hashes.json", "w"), indent=2)
print(json.dumps(meta, indent=2))
