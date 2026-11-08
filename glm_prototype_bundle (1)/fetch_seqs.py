"""Fetch full-length 16S records for the B. cereus group prototype from NCBI."""
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

def esearch(term, retmax=25):
    url = f"{BASE}/esearch.fcgi?db=nuccore&retmax={retmax}&retmode=json&term={urllib.parse.quote(term)}"
    return json.loads(_get(url))["esearchresult"]["idlist"]

def efetch(ids):
    url = f"{BASE}/efetch.fcgi?db=nuccore&rettype=fasta&retmode=text&id={','.join(ids)}"
    return _get(url)

groups = {
    "B_anthracis": "Bacillus anthracis[organism] AND 16S ribosomal RNA[title] NOT genome[title] NOT contig[title] NOT whole genome[title]",
    "B_cereus": "Bacillus cereus[organism] AND 16S ribosomal RNA[title] NOT genome[title] NOT contig[title] NOT whole genome[title]",
    "B_thuringiensis": "Bacillus thuringiensis[organism] AND 16S ribosomal RNA[title] NOT genome[title] NOT contig[title] NOT whole genome[title]",
    "B_weihenstephanensis": "Bacillus weihenstephanensis[organism] AND 16S ribosomal RNA[title] NOT genome[title] NOT contig[title]",
    "B_subtilis": "Bacillus subtilis[organism] AND 16S ribosomal RNA[title] NOT genome[title] NOT contig[title] NOT whole genome[title]",
}

records = []  # (group, header, seq)
for grp, term in groups.items():
    ids = esearch(term, retmax=40)
    print(grp, "n_ids:", len(ids), file=sys.stderr)
    fasta = efetch(ids)
    for chunk in fasta.split(">"):
        if not chunk.strip():
            continue
        header, _, seq = chunk.partition("\n")
        seq = "".join(seq.split())
        if 1350 <= len(seq) <= 1600 and set(seq.upper()) <= set("ACGTUN"):
            records.append((grp, header.strip(), seq.upper().replace("U", "T")))
    time.sleep(0.5)

# E. coli coordinate reference: J01695 (rrsA of E. coli K-12, the classical numbering source)
ec = efetch(["J01695"])
for chunk in ec.split(">"):
    if chunk.strip():
        header, _, seq = chunk.partition("\n")
        records.append(("E_coli_ref", header.strip(), "".join(seq.split()).upper().replace("U", "T")))

print(json.dumps(records), file=open("cereus_group_16s.json", "w"))
from collections import Counter
print(Counter(r[0] for r in records), file=sys.stderr)
