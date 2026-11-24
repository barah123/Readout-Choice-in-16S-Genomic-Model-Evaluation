"""How much reference data actually exists per species, by source.

Evidence for the reference-depth limit. The subsampling study
(step2_reference_incompleteness.py) gives the error rate of a naive resolvability call as a
function of how many independent organisms per species the reference contains. This script
asks the complementary question: how many are actually available.

Counts, per species in the locked-test panels, from NCBI:
  - complete RefSeq genomes            (the only source that reliably reports a full 16S copy set)
  - complete + chromosome-level
  - assemblies at any level            (adds drafts, which lose rRNA copies to repeat collapse)
and compares against the 16S deposit count already in hand.

Rate-limited to stay inside NCBI's anonymous request ceiling.
Writes results/reference_depth_by_species.json.
"""
import json, os, subprocess, time, urllib.parse
from collections import defaultdict

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
OUT = os.path.abspath(os.path.join(HERE, "..", "results", "reference_depth_by_species.json"))
E = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"


def count(term):
    url = f"{E}?db=assembly&term={urllib.parse.quote(term)}&retmax=0&retmode=json"
    for _ in range(3):
        try:
            r = subprocess.run(["curl", "-sL", "--max-time", "25", url],
                               capture_output=True, text=True)
            return int(json.loads(r.stdout)["esearchresult"]["count"])
        except Exception:
            time.sleep(1.5)
    return None


def main():
    raw = os.path.join(ROOT, "phase4_locked_test", "raw_species_test")
    by_sp = defaultdict(set)
    for f in os.listdir(raw):
        if f.endswith(".json"):
            for r in json.load(open(os.path.join(raw, f))):
                by_sp[r["species_query"]].add(r["isolate_group"])
    species = sorted(by_sp)
    print(f"{len(species)} species", flush=True)

    out = {}
    for i, s in enumerate(species, 1):
        base = f'"{s}"[Organism] AND "latest refseq"[filter]'
        out[s] = {
            "complete": count(base + ' AND "Complete Genome"[Assembly Level]'),
            "complete_or_chromosome": count(
                base + ' AND ("Complete Genome"[Assembly Level] OR "Chromosome"[Assembly Level])'),
            "any_assembly_level": count(base),
            "deposits_16S": len(by_sp[s])}
        if i % 20 == 0:
            print(f"  {i}/{len(species)}", flush=True)
        time.sleep(0.5)

    ok = {k: v for k, v in out.items() if v["complete"] is not None}
    arr = lambda k: np.array([v[k] for v in ok.values()])
    summary = {"n_species": len(ok), "queried": "NCBI assembly, latest refseq",
               "caveat": "these species passed the gLM eligibility gate (>=10 isolates per "
                         "target, >=5 per neighbour), so every column is an upper bound on "
                         "what would be seen across bacteria generally"}
    for k in ("complete", "complete_or_chromosome", "any_assembly_level", "deposits_16S"):
        x = arr(k)
        summary[k] = {"median": float(np.median(x)), "q1": float(np.percentile(x, 25)),
                      "q3": float(np.percentile(x, 75)), "zero": int((x == 0).sum()),
                      "share_at_or_above": {str(t): float((x >= t).mean())
                                            for t in (2, 5, 12, 20, 63)}}
    json.dump({"summary": summary, "by_species": out}, open(OUT, "w"), indent=2)

    print(f"\n{'source':>26}{'median':>9}{'>=12 (10% err)':>17}")
    for k in ("complete", "complete_or_chromosome", "any_assembly_level", "deposits_16S"):
        print(f"{k:>26}{summary[k]['median']:>9.0f}"
              f"{summary[k]['share_at_or_above']['12']:>16.1%}")
    print(f"\nwrote {OUT}")


if __name__ == "__main__":
    main()
