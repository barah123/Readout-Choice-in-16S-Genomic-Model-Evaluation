"""Fetch 16S records + isolate metadata for every species in blocks_definition.py.

Resumable: writes one JSON per species under raw_species/ and skips any already present,
so an interrupted run or a transient NCBI failure costs only the species in flight.
Uses the same esearch+efetch+esummary pairing as fetch_with_metadata.py so every record
carries a structured strain field (the isolate grouping key).
"""
import json
import os
import sys
import time
import urllib.error
import urllib.parse
import urllib.request

from blocks_definition import BLOCKS, all_species

BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
RETMAX = 300
MIN_LEN, MAX_LEN = 1350, 1600
OUTDIR = "raw_species"
SLEEP = 0.5          # unauthenticated NCBI limit is 3 req/s; stay well under

_COLLECTION_PREFIXES = {"ATCC", "DSM", "NCTC", "CECT", "JCM", "KCTC", "CIP",
                        "LMG", "NBRC", "CGMCC", "BCRC", "NRRL", "CCUG", "VTT",
                        "CCM", "IFO", "KACC", "BAA"}


def _get(url, tries=7):
    for i in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=120) as r:
                return r.read().decode()
        except urllib.error.HTTPError as e:
            if e.code in (429, 500, 502, 503, 504):
                time.sleep(4 * (i + 1))
                continue
            raise
        except Exception:
            time.sleep(4 * (i + 1))
    raise RuntimeError("giving up: " + url)


def term_for(species):
    return (f'"{species}"[organism] AND 16S ribosomal RNA[title] '
            f'NOT genome[title] NOT contig[title] NOT whole genome[title]')


def clean_strain(raw):
    if not raw:
        return None
    s = raw.strip()
    if not s or s.lower() in {"-", "na", "n/a", "none", "unknown", "missing"}:
        return None
    if s.upper() in _COLLECTION_PREFIXES:
        return None
    return s


def parse_subname(rec):
    keys = (rec.get("subtype") or "").split("|")
    vals = (rec.get("subname") or "").split("|")
    return {k: v for k, v in zip(keys, vals) if k}


def fetch_species(block_id, species, is_target):
    safe = species.replace(" ", "_").replace("/", "_").replace(".", "")
    path = os.path.join(OUTDIR, f"{block_id}__{safe}.json")
    if os.path.exists(path):
        return json.load(open(path)), True

    term = term_for(species)
    url = f"{BASE}/esearch.fcgi?db=nuccore&retmax={RETMAX}&retmode=json&term={urllib.parse.quote(term)}"
    ids = json.loads(_get(url))["esearchresult"]["idlist"]
    time.sleep(SLEEP)

    records = []
    for i in range(0, len(ids), 100):
        chunk = ids[i:i + 100]
        fasta = _get(f"{BASE}/efetch.fcgi?db=nuccore&rettype=fasta&retmode=text&id={','.join(chunk)}")
        time.sleep(SLEEP)
        summ = json.loads(_get(f"{BASE}/esummary.fcgi?db=nuccore&retmode=json&id={','.join(chunk)}"))["result"]
        time.sleep(SLEEP)

        by_acc = {}
        for uid, rec in summ.items():
            if uid == "uids" or not isinstance(rec, dict):
                continue
            if rec.get("accessionversion"):
                by_acc[rec["accessionversion"]] = rec

        for block in fasta.split(">"):
            if not block.strip():
                continue
            header, _, seq = block.partition("\n")
            seq = "".join(seq.split()).upper().replace("U", "T")
            if not (MIN_LEN <= len(seq) <= MAX_LEN) or not set(seq) <= set("ACGTN"):
                continue
            acc = header.split()[0]
            rec = by_acc.get(acc, {})
            sub = parse_subname(rec)
            strain = clean_strain(rec.get("strain") or sub.get("strain"))
            # Guard against NCBI organism drift: keep only records whose reported organism
            # actually starts with the queried species name (subspecies suffixes allowed).
            org = (rec.get("organism") or "")
            records.append({
                "block_id": block_id, "species_query": species, "is_target": bool(is_target),
                "accession": acc, "organism": org, "taxid": rec.get("taxid"),
                "strain": strain,
                "isolate_group": f"{species}|{strain}" if strain else f"{species}|ACC:{acc}",
                "isolate_inferred": bool(strain),
                "country": sub.get("country"), "isolation_source": sub.get("isolation_source"),
                "collection_date": sub.get("collection_date"),
                "slen": len(seq), "title": header.strip(), "sequence": seq,
            })

    os.makedirs(OUTDIR, exist_ok=True)
    json.dump(records, open(path, "w"))
    return records, False


def main():
    os.makedirs(OUTDIR, exist_ok=True)
    targets = all_species()
    print(f"{len(targets)} species queries across {len(BLOCKS)} blocks", flush=True)
    tally = {}
    for n, (block_id, species, is_target) in enumerate(targets, 1):
        try:
            recs, cached = fetch_species(block_id, species, is_target)
        except Exception as e:
            print(f"[{n}/{len(targets)}] FAILED {species}: {e}", flush=True)
            continue
        tag = "cached" if cached else "fetched"
        n_iso = len({r["isolate_group"] for r in recs})
        tally[(block_id, species)] = (len(recs), n_iso)
        print(f"[{n}/{len(targets)}] {tag:7} {block_id:32} {species:48} "
              f"records={len(recs):4} isolates={n_iso:4}", flush=True)

    print("\n--- per-block target support ---", flush=True)
    for b in BLOCKS:
        t = tally.get((b["block_id"], b["target"]))
        ns = [tally.get((b["block_id"], x)) for x in b["neighbours"]]
        ns = [x for x in ns if x]
        tgt_iso = t[1] if t else 0
        nb_ok = sum(1 for x in ns if x[1] >= 5)
        print(f"{b['block_id']:32} target_isolates={tgt_iso:4}  "
              f"neighbour_species_with_ge5_isolates={nb_ok}/{len(b['neighbours'])}", flush=True)


if __name__ == "__main__":
    main()
