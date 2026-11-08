"""Phase 2 step 1: fetch 16S records WITH isolate metadata attached.

Why this replaces fetch_seqs.py / fetch_staph.py
------------------------------------------------
Phase 1's audit (phase1_audit/audit_report.md, findings 1 and 4) established that the existing
fetch scripts keep only (species_group, accession, sequence) from a bare FASTA defline, so
genome/isolate identity is lost. That blocks protocol §12.1's requirement that copies/regions
from the same biological unit never cross a support/query split - the existing LOHO-CV holds
out identical *strings*, which is strictly weaker: two records from the same isolate differing
by one sequencing-error base land in different haplotype groups and can straddle the split.

NCBI esummary exposes a structured `strain` field (and a `subname`/`subtype` pair carrying
host/country/isolation_source/collection_date), verified available for these records. This
script pairs efetch (sequence) with esummary (metadata) so every row carries a real isolate
key rather than a regex guess off the title string.

Isolate grouping policy (declared here, before any model is run, per protocol §12.1)
------------------------------------------------------------------------------------
isolate_group = (species, strain) when a non-empty, non-placeholder strain is present.
Otherwise the record is its own singleton group keyed by accession, and flagged
`isolate_inferred=False` so a sensitivity analysis can exclude them. Placeholder strain values
(e.g. bare culture-collection prefixes like "ATCC" with no number) are treated as MISSING, not
as a shared group - lumping every "ATCC" record into one isolate would be a worse error than
splitting them, since they are genuinely different deposits.

Output: one JSON per panel with full per-record metadata, written to this directory.
Does not overwrite the Phase 0-3a caches in glm_prototype_bundle (1)/.
"""
import json
import re
import sys
import time
import urllib.parse
import urllib.request

BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
RETMAX = 300          # see CLAUDE.md: retmax=40 under-samples full-length genes badly
MIN_LEN, MAX_LEN = 1350, 1600


def _get(url, tries=6):
    for i in range(tries):
        try:
            with urllib.request.urlopen(url) as r:
                return r.read().decode()
        except urllib.error.HTTPError as e:
            if e.code in (429, 500, 502, 503):
                time.sleep(3 * (i + 1))
                continue
            raise
        except urllib.error.URLError:
            time.sleep(3 * (i + 1))
    raise RuntimeError("giving up after retries: " + url)


def esearch(term, retmax=RETMAX):
    url = (f"{BASE}/esearch.fcgi?db=nuccore&retmax={retmax}&retmode=json"
           f"&term={urllib.parse.quote(term)}")
    return json.loads(_get(url))["esearchresult"]["idlist"]


def esearch_count(term):
    url = f"{BASE}/esearch.fcgi?db=nuccore&retmax=0&retmode=json&term={urllib.parse.quote(term)}"
    return int(json.loads(_get(url))["esearchresult"]["count"])


def efetch_fasta(ids):
    url = (f"{BASE}/efetch.fcgi?db=nuccore&rettype=fasta&retmode=text"
           f"&id={','.join(ids)}")
    return _get(url)


def esummary(ids):
    url = f"{BASE}/esummary.fcgi?db=nuccore&retmode=json&id={','.join(ids)}"
    return json.loads(_get(url))["result"]


# A strain token must contain at least one digit OR be >=4 chars and not a bare
# culture-collection prefix; otherwise we treat it as missing rather than merging
# unrelated deposits under a shared placeholder.
_COLLECTION_PREFIXES = {"ATCC", "DSM", "NCTC", "CECT", "JCM", "KCTC", "CIP",
                        "LMG", "NBRC", "CGMCC", "BCRC", "NRRL", "CCUG", "VTT"}


def clean_strain(raw):
    if not raw:
        return None
    s = raw.strip()
    if not s or s.lower() in {"-", "na", "n/a", "none", "unknown", "missing"}:
        return None
    if s.upper() in _COLLECTION_PREFIXES:   # bare prefix, no catalogue number
        return None
    return s


def parse_subname(rec):
    """esummary gives parallel 'subtype' (pipe-separated keys) and 'subname' (values)."""
    keys = (rec.get("subtype") or "").split("|")
    vals = (rec.get("subname") or "").split("|")
    return {k: v for k, v in zip(keys, vals) if k}


def fetch_panel(groups, out_path, panel_name):
    records = []
    for grp, term in groups.items():
        total = esearch_count(term)
        time.sleep(0.4)
        ids = esearch(term)
        time.sleep(0.4)
        print(f"{grp}: NCBI total hits={total}, retrieved ids={len(ids)}", file=sys.stderr)
        if not ids:
            continue

        # efetch + esummary in matched sub-batches of 100
        for i in range(0, len(ids), 100):
            chunk = ids[i:i + 100]
            fasta = efetch_fasta(chunk)
            time.sleep(0.4)
            summ = esummary(chunk)
            time.sleep(0.4)

            # index esummary by accession.version for joining to FASTA deflines
            by_acc = {}
            for uid, rec in summ.items():
                if uid == "uids" or not isinstance(rec, dict):
                    continue
                acc = rec.get("accessionversion")
                if acc:
                    by_acc[acc] = rec

            for block in fasta.split(">"):
                if not block.strip():
                    continue
                header, _, seq = block.partition("\n")
                seq = "".join(seq.split()).upper().replace("U", "T")
                if not (MIN_LEN <= len(seq) <= MAX_LEN):
                    continue
                if not set(seq) <= set("ACGTN"):
                    continue
                acc = header.split()[0]
                rec = by_acc.get(acc, {})
                sub = parse_subname(rec)
                strain = clean_strain(rec.get("strain") or sub.get("strain"))
                records.append({
                    "panel": panel_name,
                    "species_group": grp,
                    "accession": acc,
                    "organism": rec.get("organism"),
                    "taxid": rec.get("taxid"),
                    "strain": strain,
                    "isolate_group": f"{grp}|{strain}" if strain else f"{grp}|ACC:{acc}",
                    "isolate_inferred": bool(strain),
                    "host": sub.get("host"),
                    "country": sub.get("country"),
                    "isolation_source": sub.get("isolation_source"),
                    "collection_date": sub.get("collection_date"),
                    "createdate": rec.get("createdate"),
                    "slen": len(seq),
                    "title": header.strip(),
                    "sequence": seq,
                })
        print(f"{grp}: kept {sum(1 for r in records if r['species_group'] == grp)} "
              f"full-length records", file=sys.stderr)

    with open(out_path, "w") as f:
        json.dump(records, f)
    print(f"wrote {out_path}: {len(records)} records", file=sys.stderr)
    return records


def q(species):
    return (f"{species}[organism] AND 16S ribosomal RNA[title] "
            f"NOT genome[title] NOT contig[title] NOT whole genome[title]")


if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "all"

    if which in ("all", "staph"):
        fetch_panel({
            "S_aureus": q("Staphylococcus aureus"),
            "S_epidermidis": q("Staphylococcus epidermidis"),
            "S_haemolyticus": q("Staphylococcus haemolyticus"),
            "S_saprophyticus": q("Staphylococcus saprophyticus"),
            "S_hominis": q("Staphylococcus hominis"),
            "S_lugdunensis": q("Staphylococcus lugdunensis"),
            "S_capitis": q("Staphylococcus capitis"),
        }, "staph_panel_meta.json", "staph_CoNS")

    if which in ("all", "cereus"):
        fetch_panel({
            "B_anthracis": q("Bacillus anthracis"),
            "B_cereus": q("Bacillus cereus"),
            "B_thuringiensis": q("Bacillus thuringiensis"),
            "B_weihenstephanensis": q("Bacillus weihenstephanensis"),
        }, "cereus_panel_meta.json", "cereus_group")
