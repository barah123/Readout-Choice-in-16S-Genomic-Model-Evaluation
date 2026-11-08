"""Fetch full-length 16S records for the S. aureus / CoNS panel from NCBI.

Target: Staphylococcus aureus (Panel A, plan line 42: "V1 best" per Chakravorty et al. 2007).
Near neighbors: the six coagulase-negative staphylococci most commonly used as the comparison
set in 16S-based Staphylococcus differentiation studies (S. epidermidis, S. haemolyticus,
S. saprophyticus, S. hominis, S. lugdunensis, S. capitis). No outgroup fetched yet - this run
is only for the haplotype pre-check (MASTER_PLAN.md §4 item 2); add an outgroup when/if this
target clears the gate and moves to the full B1-B4 + gLM comparison.
"""
import json, time, sys, urllib.parse, urllib.request

BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

def _get(url, tries=6):
    for i in range(tries):
        try:
            with urllib.request.urlopen(url) as r:
                return r.read().decode()
        except urllib.error.HTTPError as e:
            if e.code == 429:
                time.sleep(3 * (i + 1))
                continue
            raise
    raise RuntimeError("giving up after retries: " + url)

def esearch(term, retmax=300):
    url = f"{BASE}/esearch.fcgi?db=nuccore&retmax={retmax}&retmode=json&term={urllib.parse.quote(term)}"
    return json.loads(_get(url))["esearchresult"]["idlist"]

def efetch(ids):
    url = f"{BASE}/efetch.fcgi?db=nuccore&rettype=fasta&retmode=text&id={','.join(ids)}"
    return _get(url)

groups = {
    "S_aureus": "Staphylococcus aureus[organism] AND 16S ribosomal RNA[title] NOT genome[title] NOT contig[title] NOT whole genome[title]",
    "S_epidermidis": "Staphylococcus epidermidis[organism] AND 16S ribosomal RNA[title] NOT genome[title] NOT contig[title] NOT whole genome[title]",
    "S_haemolyticus": "Staphylococcus haemolyticus[organism] AND 16S ribosomal RNA[title] NOT genome[title] NOT contig[title] NOT whole genome[title]",
    "S_saprophyticus": "Staphylococcus saprophyticus[organism] AND 16S ribosomal RNA[title] NOT genome[title] NOT contig[title] NOT whole genome[title]",
    "S_hominis": "Staphylococcus hominis[organism] AND 16S ribosomal RNA[title] NOT genome[title] NOT contig[title] NOT whole genome[title]",
    "S_lugdunensis": "Staphylococcus lugdunensis[organism] AND 16S ribosomal RNA[title] NOT genome[title] NOT contig[title] NOT whole genome[title]",
    "S_capitis": "Staphylococcus capitis[organism] AND 16S ribosomal RNA[title] NOT genome[title] NOT contig[title] NOT whole genome[title]",
}

records = []  # (group, header, seq)
for grp, term in groups.items():
    ids = esearch(term, retmax=300)
    time.sleep(0.4)
    print(grp, "n_ids:", len(ids), file=sys.stderr)
    if not ids:
        continue
    # efetch in sub-batches of 100 ids - a single request for 300 ids risks timing out / 429
    fasta_chunks = []
    for i in range(0, len(ids), 100):
        fasta_chunks.append(efetch(ids[i:i + 100]))
        time.sleep(0.4)
    fasta = "".join(fasta_chunks)
    n_kept = 0
    for chunk in fasta.split(">"):
        if not chunk.strip():
            continue
        header, _, seq = chunk.partition("\n")
        seq = "".join(seq.split())
        if 1350 <= len(seq) <= 1600 and set(seq.upper()) <= set("ACGTUN"):
            records.append((grp, header.strip(), seq.upper().replace("U", "T")))
            n_kept += 1
    print(grp, "n_full_length_kept:", n_kept, file=sys.stderr)

print(json.dumps(records), file=open("staph_group_16s.json", "w"))
from collections import Counter
print(Counter(r[0] for r in records), file=sys.stderr)
