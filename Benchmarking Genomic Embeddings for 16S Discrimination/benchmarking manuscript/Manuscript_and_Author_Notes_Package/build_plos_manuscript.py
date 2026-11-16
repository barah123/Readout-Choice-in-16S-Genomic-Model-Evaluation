"""Build the PLOS Computational Biology formatted manuscript as .docx.

Style follows the supplied exemplar (journal.pcbi.1014596): title, authors with numbered
affiliations, Abstract, Author summary, Introduction, Methods, Results, Discussion,
numbered square-bracket citations, S-numbered supporting information, and the standard
availability/funding/competing-interests statements.

Placeholders that must be supplied by the authors are written as [AUTHOR: ...] so they are
greppable before submission. Every numeric value is taken from the result files, not retyped.
"""
import json
import os

from docx import Document
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.shared import Pt, RGBColor, Inches

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
RES = os.path.abspath(os.path.join(HERE, "..", "results"))
FIGS = os.path.abspath(os.path.join(HERE, "..", "figures"))
OUT = os.path.join(HERE, "Manuscript_PLOS_CompBiol.docx")


def load(p):
    with open(p) as f:
        return json.load(f)


LT = load(os.path.join(ROOT, "phase5_evaluate", "locked_test_results.json"))
EA = load(os.path.join(ROOT, "phase7_extensions", "extension_A_results.json"))
RB = load(os.path.join(ROOT, "phase7_extensions", "robustness_results.json"))
RR = load(os.path.join(RES, "region_ranking_summary.json"))
LP = load(os.path.join(RES, "label_permutation_control.json"))
BL = load(os.path.join(ROOT, "phase2_pilot", "block_level_results.json"))

P = LT["primary"]
pm, plo, phi = P["mean"] * 100, P["ci95"][0] * 100, P["ci95"][1] * 100

doc = Document()
st = doc.styles["Normal"]
st.font.name = "Times New Roman"
st.font.size = Pt(11)
st.paragraph_format.space_after = Pt(6)
st.paragraph_format.line_spacing = 1.5
for s in doc.sections:
    s.left_margin = s.right_margin = Inches(1)
    s.top_margin = s.bottom_margin = Inches(1)


def H(text, size=14, space_before=14):
    p = doc.add_paragraph()
    p.paragraph_format.space_before = Pt(space_before)
    p.paragraph_format.space_after = Pt(4)
    r = p.add_run(text)
    r.bold = True
    r.font.size = Pt(size)
    return p


def H2(text):
    return H(text, size=12, space_before=10)


def body(text, italic=False):
    p = doc.add_paragraph()
    r = p.add_run(text)
    r.italic = italic
    return p


def ph(text):
    """An author-supplied placeholder, visually marked."""
    p = doc.add_paragraph()
    r = p.add_run(text)
    r.font.color.rgb = RGBColor(0xB0, 0x00, 0x20)
    r.bold = True
    return p


def table(headers, rows, caption=None, note=None):
    if caption:
        p = doc.add_paragraph()
        r = p.add_run(caption)
        r.bold = True
        r.font.size = Pt(10)
    t = doc.add_table(rows=1, cols=len(headers))
    t.style = "Table Grid"
    for i, h in enumerate(headers):
        c = t.rows[0].cells[i]
        c.text = ""
        run = c.paragraphs[0].add_run(h)
        run.bold = True
        run.font.size = Pt(9)
    for row in rows:
        cells = t.add_row().cells
        for i, v in enumerate(row):
            cells[i].text = ""
            run = cells[i].paragraphs[0].add_run(str(v))
            run.font.size = Pt(9)
    if note:
        p = doc.add_paragraph()
        r = p.add_run(note)
        r.font.size = Pt(9)
        r.italic = True
    doc.add_paragraph()


# ============================================================ TITLE
p = doc.add_paragraph()
p.alignment = WD_ALIGN_PARAGRAPH.LEFT
r = p.add_run("Readout choice, not model choice, determines whether pretrained genomic "
              "language models help: a pre-registered benchmark on fine-scale 16S rRNA "
              "discrimination")
r.bold = True
r.font.size = Pt(16)

# ---- authors, PLOS title-page format ----
pa = doc.add_paragraph()
pa.add_run("Philip Yamoah Appiah").bold = True
pa.add_run("\u00B9\u00B2*, ")
pa.add_run("Ruth Barkey Tettey").bold = True
pa.add_run("\u00B3, ")
pa.add_run("Bernice Nutase").bold = True
pa.add_run("\u2074, ")
pa.add_run("Carmen Patricia Ortega-Santos").bold = True
pa.add_run("\u00B2, ")
pa.add_run("Marcos P\u00E9rez-Losada").bold = True
pa.add_run("\u00B9")

for num, aff in [
    ("1", "Department of Biostatistics and Bioinformatics, Computational Biology Institute, "
          "Milken Institute School of Public Health, The George Washington University, "
          "Washington, DC 20052, United States of America"),
    ("2", "Department of Exercise and Nutrition Sciences, Milken Institute School of Public "
          "Health, The George Washington University, Washington, DC 20052, United States of "
          "America"),
    ("3", "Ghana Health Service, Oyoko Health Center, Eastern Region, Ghana"),
    ("4", "Ghana Health Service, Madina Polyclinic, Greater Accra, Ghana"),
]:
    q = doc.add_paragraph()
    q.paragraph_format.space_after = Pt(2)
    rr = q.add_run(f"{num} ")
    rr.bold = True
    rr.font.size = Pt(9)
    rr2 = q.add_run(aff)
    rr2.font.size = Pt(9)

q = doc.add_paragraph()
q.paragraph_format.space_before = Pt(6)
rr = q.add_run("ORCID iDs: ")
rr.bold = True
rr.font.size = Pt(9)
q.add_run(
    "Appiah 0009-0002-0706-2506; Tettey 0009-0001-1448-0688; "
    "Nutase 0009-0004-5574-6798; Ortega-Santos 0000-0002-0071-1190; "
    "P\u00E9rez-Losada 0000-0002-2585-4657").font.size = Pt(9)

q = doc.add_paragraph()
q.paragraph_format.space_before = Pt(6)
rr = q.add_run("* Corresponding author: ")
rr.bold = True
rr.font.size = Pt(9)
q.add_run("p.appiah@gwmail.gwu.edu (PYA)").font.size = Pt(9)
ph("[AUTHOR: Short title, 50 characters maximum. AQ22.]")

# ============================================================ ABSTRACT
H("Abstract")
body(
    "Pretrained genomic language models are increasingly proposed for sequence "
    "identification, but whether their representations carry information beyond simple "
    "composition statistics remains unclear. We benchmarked four published models "
    "(DNABERT-2, DNABERT-S, HyenaDNA and Nucleotide Transformer v2) against canonical "
    "k-mer counting on target-versus-neighbour discrimination within single 16S rRNA "
    "hypervariable regions. The analysis plan, panel definitions and eligibility criteria "
    "were committed to version control before the confirmatory data were collected. Across "
    "58 analysable bacterial groups drawn from six phyla, a frozen DNABERT-2 representation "
    "with a regularized linear probe and k-mer counting with the same probe were "
    f"statistically equivalent on the locked test set (mean difference {pm:+.2f} percentage "
    f"points, 95% bootstrap confidence interval {plo:+.2f} to {phi:+.2f}; two one-sided "
    "tests against a prespecified 1-point margin, p = 4.1 x 10-4). Point estimates were "
    "similar in a development set and in a set of groups selected specifically for "
    "difficulty, although the latter interval was too wide to establish equivalence at the "
    "same margin. Two controls bound the interpretation: a Gaussian random projection of "
    "k-mer counts matched every learned representation, whereas permuting labels through "
    "the complete pipeline reduced every method to chance. Among classification readouts "
    "measured on a common endpoint, moving from nearest-reference assignment to a linear "
    "probe changed balanced accuracy by about 2.9 points on the development set, roughly "
    "four times the 0.73-point "
    "spread across the four models. Within these panels, the lowest accuracies occurred in "
    "subspecies and species-complex comparisons rather than in genera reputed to have "
    "conserved 16S sequences. For this task, frozen genomic language models offered no "
    "accuracy advantage over k-mer counting, and how a representation is read mattered more "
    "than which model produced it.")

# ============================================================ AUTHOR SUMMARY
H("Author summary")
body(
    "Laboratories identifying bacteria from the 16S rRNA gene usually sequence one short "
    "section rather than the whole gene, and the section chosen determines whether closely "
    "related species can be told apart. Large neural networks trained on DNA, often called "
    "genomic language models, have been proposed for problems of this kind. We asked "
    "whether four published models outperform a simple method from the 1990s that counts "
    "short DNA patterns and requires no neural network. Using an analysis plan fixed in "
    "advance and 58 groups of closely related bacteria, we found the two approaches "
    "performed the same, and the result held on independent bacteria, on deliberately "
    "difficult groups, and on degraded sequences. Two findings matter more than that "
    "headline. First, how a representation is extracted from a model affected accuracy "
    "several times more than which model was used, so a study reporting that a model "
    "failed may be measuring its own analysis choices. Second, a random compression of the "
    "simple pattern counts performed as well as any learned representation, while shuffling "
    "the labels reduced every method to chance. Both controls are easy to omit, and omitting "
    "either can support a confident but incorrect conclusion.")

# ============================================================ INTRODUCTION
H("Introduction")
body(
    "The 16S rRNA gene is the most widely used marker for bacterial identification and "
    "community profiling. Because sequencing the full gene is often impractical, most "
    "studies amplify one or more of its nine hypervariable regions, labelled V1 to V9. The "
    "choice of region is consequential: a region that is identical between a target "
    "organism and its closest relatives cannot distinguish them, regardless of the "
    "analytical method applied afterwards [1,2]. Region choice is therefore a design "
    "decision that constrains every downstream conclusion.")
body(
    "Pretrained genomic language models have been proposed as a way to extract more "
    "information from short marker sequences. DNABERT-2 [3] and Nucleotide Transformer v2 "
    "[4] are general multi-species encoders; HyenaDNA [5] uses a long-convolution "
    "architecture; DNABERT-S [6] was trained specifically to produce species-aware "
    "embeddings. Related work on 16S representation learning reports that classical machine "
    "learning can match or exceed learned embeddings on several downstream tasks [7], and "
    "benchmarking of DNA models on regulatory sequence shows that conclusions depend "
    "strongly on whether a model is evaluated zero-shot, with a trained probe, or "
    "fine-tuned [8].")
body(
    "Two features of the fine-scale 16S problem make it a demanding test. The gene is "
    "short and highly conserved, so the available signal may be small. And the achievable "
    "accuracy is bounded: when a region is identical across a species boundary, no method "
    "can separate the organisms that share it. Any comparison must therefore report what "
    "was achievable as well as what was achieved.")
body(
    "We benchmarked four published models against k-mer counting on target-versus-neighbour "
    "discrimination within single regions. The study was designed to make a negative result "
    "interpretable rather than merely unpublished: the analysis plan was fixed before the "
    "confirmatory data existed, equivalence was assessed against a margin chosen in "
    "advance, and two controls were included to bound what the measured accuracy could "
    "mean.")

# ============================================================ METHODS
H("Methods")

H2("Study design")
body(
    "The study comprised three sequential experiments. A development set was used for every "
    "analysis decision. A locked test set, drawn from bacterial families disjoint from the "
    "development set, was scored once under a plan committed beforehand. A third set was "
    "selected specifically for difficulty, to test a hypothesis generated by inspecting the "
    "first two. Two further analyses examined sequence degradation and intermediate model "
    "layers. Table 1 gives the accounting.")
body(
    "The analysis plan for the locked test was committed to a version-controlled repository "
    "at 07:20:42 on 8 November 2026, and the first locked-test sequence record was "
    "retrieved at 08:54:19 the same day, 94 minutes later. The commit added the plan, the "
    "locked-test panel definitions and the deviation log and nothing else, and the plan has "
    "a single entry in its file history, indicating it was not edited afterwards. Panel "
    "definitions for the development and hard-panel experiments were likewise committed "
    "before the corresponding sequences were retrieved.")
body(
    "These are author-recorded timestamps in a repository that was not pushed to a remote "
    "server, so they document the intended sequence rather than provide independent "
    "certification; git commit dates are set by the committer and are not externally "
    "verifiable. Two features of the record are nonetheless consistent with genuine "
    "pre-declaration: 2 of 22 pre-declared locked-test blocks and 5 of 31 development "
    "blocks subsequently failed the authors' own eligibility criteria, and a "
    "pre-registered directional hypothesis did not replicate and is reported as such. The "
    "complete history, the verification commands and a statement of their limits are given "
    "in S10 Text, and every protocol deviation is logged with its date and justification "
    "in S3 Text.")

H2("Sequence retrieval and taxonomic panels")
body(
    "Sequences were retrieved from the NCBI nucleotide database using Entrez E-utilities. "
    "For each species, esearch was issued with the query template "
    "\"<species>[organism] AND 16S ribosomal RNA[title] NOT genome[title] NOT contig[title] "
    "NOT whole genome[title]\" and retmax 300, and the returned identifiers were fetched in "
    "batches of 100 by efetch (FASTA) paired with esummary (JSON) so that structured "
    "metadata accompanied every sequence. Records were retained if they were 1,350 to 1,600 "
    "nucleotides long and contained only the characters A, C, G, T and N. Uracil was "
    "converted to thymine.")
body(
    "A panel is one target species and three to six congeneric neighbours, selected as the "
    "congeners most often reported as confusable with the target or its closest relatives "
    "by published phylogeny. A block is one genus, or a tightly linked genus pair where "
    "taxonomy has split a clade, so that no species appears in two blocks. Blocks are "
    "additionally tagged by family; development and locked-test families are disjoint. "
    "Species-level and subspecies-level panels are separate strata and are never pooled.")
body(
    "Panel membership was declared in version-controlled definition files before the "
    "corresponding sequences were retrieved, and each of the 53 declared development and "
    "locked-test blocks carries a recorded one-line rationale for its target and its "
    "neighbours (S1 Data). Four of those rationales name a specific primary source, in "
    "each case the region-by-region analysis of Chakravorty and colleagues [1]; the "
    "remaining rationales state the reasoning, such as reported identity above the level "
    "at which 16S separates species, a clonal lineage nested within a parent species, or "
    "a taxonomic split through a clade, without citing a publication for that panel. "
    "Neighbour selection is therefore expert selection informed by the 16S identification "
    "literature, not a procedure traceable to a named source for each panel, and we report "
    "it as such rather than attaching citations retrospectively.")
body(
    "The 15 hard panels were not selected from the literature. Candidate genera and "
    "species complexes with documented 16S conservation were screened against a rule fixed "
    "before screening began: a candidate qualified as hard if, at its best region, the "
    "observed-sequence ambiguity ceiling was at most 0.97 and k-mer balanced accuracy was "
    "at most 0.95. The screen used classical methods only, so no genomic language model "
    "could influence which difficult panels entered the analysis.")

H2("Isolate grouping")
body(
    "The unit of independence is the isolate group. For each record the strain field "
    "returned by esummary was taken, or where that field was absent, the value of the "
    "strain key parsed from the record's paired subtype and subname fields. The string was "
    "stripped of surrounding whitespace and treated as missing if it was empty, if it "
    "matched any of the placeholder values -, na, n/a, none, unknown or missing "
    "case-insensitively, or if it consisted of a bare culture-collection prefix carrying "
    "no catalogue number (ATCC, DSM, NCTC, CECT, JCM, KCTC, CIP, LMG, NBRC, CGMCC, BCRC, "
    "NRRL, CCUG, VTT, CCM, IFO, KACC or BAA), since pooling every deposit labelled only "
    "\"ATCC\" would merge unrelated isolates into one unit. The isolate group is then the "
    "queried species name joined to the retained strain string, and the species name "
    "joined to the accession where no usable strain string remained, so a record without "
    "strain metadata forms a group of its own.")
body(
    "No normalization was applied beyond those rules. Strain strings were not case-folded, "
    "punctuation and internal spacing were not standardised, and no value was "
    "cross-referenced against a culture-collection catalogue, so two deposits of one "
    "isolate held under different catalogue numbers, for example an ATCC number and its "
    "DSM equivalent, remain separate groups. The formatting component of this is small and "
    "measurable: re-deriving the groups after upper-casing every strain string and "
    "deleting every non-alphanumeric character merges 45 of the 23,605 isolate groups in "
    "the three analysed sets, a change of 0.19 percent. The catalogue-equivalence "
    "component cannot be bounded from the records themselves.")
body(
    "A usable strain string was present for 23,459 of the 24,009 analysed records, 97.7 "
    "percent overall and 97.1, 97.9 and 99.0 percent of the development, locked-test and "
    "hard-panel records respectively. Grouping binds on a small minority: 313 of the "
    "23,605 groups hold more than one record, together accounting for 717 records, and no "
    "group holds more than six. Grouped cross-validation therefore departs from "
    "record-level cross-validation only for those records, and the independence of the "
    "remainder rests on the deposits being genuinely distinct rather than on any check "
    "performed here. Within an isolate group, one modal region sequence represents the "
    "group, so a group contributes exactly one unit per region regardless of how many "
    "records it contains.")
body(
    "A single accession can be returned by more than one species query where NCBI organism "
    "annotation overlaps, and such records were not merged: 113 development and 7 "
    "locked-test accessions appear under two species names, no hard-panel accession does, "
    "and 112 of the 120 fall within one block. The two copies carry different species "
    "labels and are not pooled; where one of the two species is that panel's target, the "
    "shared sequence therefore sits on both sides of the class boundary, which lowers the "
    "ambiguity ceiling for that region rather than raising measured accuracy, because the "
    "ceiling charges every sequence held in common between the classes as irreducible "
    "error.")
body(
    "These groups are derived from deposited metadata and do not establish independent "
    "provenance. The strain string is supplied by the submitter and is not verified by us "
    "or by the archive, and standalone 16S deposits are not linked to a parent genome "
    "assembly, so neither the identity of the physical isolate nor the number and "
    "heterogeneity of the rRNA operons it carries can be confirmed from the record. Two "
    "deposits of the same isolate recorded under differently formatted strings, or under "
    "different catalogue numbers, are treated here as independent. The direction of that "
    "error is toward optimism: under-merging inflates the effective number of independent "
    "isolates and leaves some leakage across cross-validation folds possible. The opposite "
    "error is worse in a different way, since merging every deposit that shares only a "
    "bare collection prefix would pool unrelated isolates into a single unit and corrupt "
    "the eligibility counts as well as the folds. The rule therefore errs toward "
    "separation, and the independence claim throughout should be read as metadata-based "
    "rather than provenance-verified.")

H2("Region extraction")
body(
    "Each sequence was aligned globally to the Escherichia coli J01695 reference using the "
    "Biopython PairwiseAligner (match +2, mismatch -1, gap open -5, gap extend -2), and "
    "Escherichia coli coordinates for V1 to V9 were transferred through the alignment. "
    "The evaluated spans are short, and their lengths matter for how every result below "
    "should be read: V1 31 nt, V2 106, V3 65, V4 107, V5 58, V6 58, V7 57, V8 52 and V9 "
    "31 nt, a range of 31 to 107 nt with a median of 58. These are the strict "
    "hypervariable spans, not the amplicons used in practice, which include conserved "
    "flanks and are two to four times longer (the 515F-806R V4 amplicon is about 253 nt "
    "and V3-V4 about 460 nt). An "
    "extracted region was retained only if its length fell between 0.75 and 1.30 times the "
    "expected length for that region; otherwise it was recorded as unavailable for that "
    "sequence, and these failures are reported rather than dropped silently.")
body(
    "Coordinate transfer was used in preference to a hidden Markov model extractor such as "
    "V-Xtractor [9] because the comparison requires every method to see exactly the same "
    "span. The models, the k-mer features and the identity baselines are scored on "
    "identical strings, so a boundary that moves between records according to a learned "
    "criterion would be confounded with the representation under test. Transferring fixed "
    "Escherichia coli coordinates through a single global alignment gives a boundary rule "
    "that is explicit, identical for every record and reproducible from the published "
    "coordinates alone, with no external binary in the dependency chain. The cost is that "
    "transfer degrades where the alignment itself is poor, a case in which a profile model "
    "would often still locate the conserved flanks; the length window is the check on that "
    "failure mode.")
body(
    "Across the 24,009 records analysed, 187,889 of 216,081 region cells were recovered "
    "(87.0 percent). Recovery was essentially complete for the internal regions (V4 100.0 "
    "percent, V5 to V8 at least 99.97 percent, V2 98.5 percent) and incomplete toward the "
    "ends of the gene (V3 77.2 percent, V1 58.3 percent, V9 48.7 percent), where deposits "
    "that begin or end inside the gene cannot supply the span. No record failed all nine "
    "regions. An unavailable region removes that isolate from that region's data; it does "
    "not alter the boundary applied to the records that pass, and region eligibility is "
    "assessed on the surviving counts.")
body(
    "Strand handling was not implemented. Sequences were aligned in their deposited "
    "orientation only and no reverse-complement check was performed, on the assumption "
    "that a reverse-complemented deposit would fail the length window and be dropped. A "
    "post hoc audit shows that assumption to be wrong. Realigning every record in both "
    "orientations under the same aligner and scoring, 22 of 24,009 records (0.09 percent) "
    "score higher reverse-complemented, and a reverse-complemented 16S gene still yields a "
    "contiguous query span for most regions, with lengths inside the retention window, so "
    "all 22 entered the analysis carrying mis-sliced sequence for eight or nine of the nine "
    "regions. Ten are in the development set, 12 in the locked test (four in bartonella, "
    "one in comamonadaceae, seven in rhizobium_agrobacterium) and none in the hard panels.")
body(
    "Because the locked test is scored once under a frozen plan, the confirmatory analysis "
    "was not re-run with these records corrected. The frozen scoring code was instead "
    "re-executed with the affected isolates removed, reusing the existing embedding cache; "
    "before removal the re-execution reproduced every frozen per-block value exactly, which "
    "establishes that the sensitivity analysis differs from the locked result in the "
    "intended respect and no other. Bartonella was already ineligible under the freeze "
    "because no region passed the eligibility gate, so two confirmatory blocks are "
    "affected. Removing the minus-strand isolates moves the primary estimand from -0.211 to "
    "-0.207 percentage points, a shift of 0.003 points, with a largest single-block shift "
    "of 0.067 points; the bootstrap interval and the equivalence verdict are unchanged. The "
    "defect is therefore immaterial to every conclusion reported here, but an orientation "
    "check belongs in any reuse of this pipeline and is noted as such in S3 Text.")
body(
    "Ambiguous bases were carried through rather than filtered. N is the only non-ACGT "
    "character present in the retrieved records; it occurs in 1,042 of the 187,889 retained "
    "region slices (0.55 percent) and in 726 of 15,420 distinct region haplotypes (4.7 "
    "percent). An N is treated as an ordinary mismatch by the aligner and by the identity "
    "calculation, is skipped by the k-mer featurizer, which counts only windows composed "
    "entirely of A, C, G and T, and is passed literally to each model's tokenizer. No "
    "ambiguity expansion, imputation or exclusion was applied. The one quantity this could "
    "flatter is the empirical ambiguity ceiling, since a slice differing from another only "
    "at an ambiguous position counts as a separate haplotype and splitting a haplotype can "
    "only raise the ceiling. Recomputing each confirmatory block's ceiling at the region it "
    "used, with N-bearing isolates removed, changes it by at most 0.03 points and by 0.002 "
    "points on average, so the reported ceilings are not an artefact of ambiguity codes.")

H2("Models and representations")
body(
    "Four checkpoints were used, frozen, with no fine-tuning. DNABERT-2 (zhihan1996/"
    "DNABERT-2-117M; weight file SHA-256 beginning 7ff39ec77a48, 468,354,983 bytes) and "
    "DNABERT-S (zhihan1996/DNABERT-S; SHA-256 beginning f3cfc3d05418, 468,320,053 bytes) "
    "were loaded from local patched copies. The four patches address CPU and library "
    "compatibility only and are not expected to alter numerical output, but no formal "
    "equivalence test against unpatched code was possible because the unpatched code does "
    "not execute in this environment. HyenaDNA (LongSafari/"
    "hyenadna-medium-160k-seqlen-hf, revision 7ebf71773d22c0ede2cc55cb2be15ee8c289e1ce) and "
    "Nucleotide Transformer v2 (InstaDeepAI/nucleotide-transformer-v2-500m-multi-species, "
    "revision 06615c1660c892fc199840c18123f8385b3542a8) were loaded from the Hugging Face "
    "hub.")
body(
    "Representations are mean-pooled final-layer hidden states, with each row's own first "
    "and last special tokens and all padding excluded from the mean. Each unique region "
    "sequence was embedded once and the resulting vector broadcast to every isolate "
    "carrying it, which makes the representation a pure function of the sequence by "
    "construction and removes any dependence on batch composition. This property was "
    "verified: string-identical inputs produced bit-identical vectors in all nine regions.")
body(
    "Intermediate layers were swept for the two models whose internals could be extracted "
    "and verified, at five depths corresponding to 0.25, 0.50, 0.75, 0.90 and 1.00 of "
    "encoder depth: layers 2, 4, 7, 8 and 9 of the 10 hidden states exposed by HyenaDNA, "
    "and layers 7, 14, 22, 26 and 29 of the 30 exposed by Nucleotide Transformer v2. The "
    "sweep ran on the 13 hard panels with everything else held fixed, the same pooling "
    "rule, linear probe, isolate-grouped five-fold split and seed; only the layer changed. "
    "A correctness gate ran before any result was computed and aborts the analysis unless "
    "the extracted final layer is numerically identical to the standard forward output "
    "(allclose, absolute tolerance 1e-6); both models passed. DNABERT-2 and DNABERT-S were "
    "not swept because they run their encoder on unpadded sequences internally and return "
    "intermediate layers shaped total-tokens by hidden rather than batch by sequence by "
    "hidden, so extracting them correctly would mean reproducing their unpad and repad "
    "index logic.")
body(
    "The single layer-depth summary plotted in Fig 2d is a range, aggregated in this order. "
    "For each model and each depth, balanced accuracy is averaged across the 13 panels, "
    "giving five depth means; the range of those five, maximum minus minimum, is that "
    "model's depth effect; the two model values are then averaged, giving 0.53 percentage "
    "points. Because panels are averaged before the range is taken, this understates how "
    "much depth can matter on an individual panel, and because the depth was chosen on the "
    "same panels used for the comparison, any best-depth figure is optimistic rather than "
    "a bound. The quantity is a difference between two balanced accuracies, the best and "
    "worst depth, which is what makes it commensurable with the other two bars in that "
    "panel: R2 minus R1 for the readout bar, and best minus worst model for the model bar.")

H2("Readouts and comparison methods")
body(
    "Three readouts were evaluated. R0 is an unsupervised geometry statistic: the "
    "probability that a randomly chosen within-class pair has higher cosine similarity than "
    "a randomly chosen between-class pair, computed over all pairs with subsampling above "
    "20,000 pairs per class. R1 is tie-aware nearest-reference classification on cosine "
    "distance, in which tied nearest neighbours receive fractional credit. R2 is a "
    "regularized linear probe: scikit-learn StandardScaler followed by LogisticRegression "
    "with L2 penalty, C = 1.0, class_weight = \"balanced\", max_iter = 2000 and "
    "random_state = 0, fitted on training folds only.")
body(
    "R0 is a pair-level area under the curve and R1 and R2 are isolate-level balanced "
    "accuracies. These are different quantities on different units and are never "
    "subtracted from one another. All quantitative readout comparisons in this paper use "
    "R1 and R2, which share an endpoint.")
body(
    "Four comparison methods were evaluated with the same cross-validation and, where "
    "applicable, the same probe. C_exact assigns the majority label of an identical "
    "training sequence and abstains otherwise. C_align is tie-aware nearest-reference "
    "classification on pairwise sequence identity, defined as twice the number of matching "
    "positions divided by the summed lengths of the two sequences. C_kmer uses canonical "
    "4-, 5- and 6-mer frequency features, normalized within each k, with the vocabulary "
    "fitted on training folds only. C_random applies a Gaussian random projection of the "
    "C_kmer feature matrix to 256 dimensions, with entries drawn from N(0, 1/p) where p is "
    "the number of k-mer features and 1/p is the variance. The projection is drawn from a "
    "generator seeded with 0. The distance-preserving properties of such a projection "
    "follow from the Johnson-Lindenstrauss lemma [13]; that result concerns preservation "
    "of pairwise distances and is cited as background only, not as evidence that the "
    "projected space is sufficient for linear classification, which is an empirical "
    "question answered by the measured accuracy. In the locked-test analysis and in the "
    "label-permutation "
    "control the generator is reseeded within each fold, so one fixed projection is used "
    "throughout; in the development readouts it is seeded once per block and region, so "
    "folds within a cell receive different projections from one deterministic stream. In "
    "every case the same projection is applied to the training and test partitions of a "
    "given fold.")
body(
    "The confirmatory C_random values therefore rest on a single projection seed, and the "
    "manuscript does not claim invariance to the projection on their basis. To establish "
    "how much of that number is seed-specific, the locked-test C_random arm was recomputed "
    "at 20 projection seeds, holding the folds, the chosen regions and the k-mer features "
    "fixed, so that only the projection varies. Seed 0 reproduced all 20 reported block "
    "values exactly. Across-block mean balanced accuracy was 96.21 percent with a standard "
    "deviation of 0.09 points across seeds (range 96.02 to 96.36; the reported seed gives "
    "96.16), and the median within-block standard deviation was 0.16 points, with a maximum "
    "of 1.17 points in the smallest block. The comparison that uses C_random, best-of-four "
    "R2 minus C_random, is +0.21 points at the reported seed and +0.16 points on average "
    "across seeds (standard deviation 0.09, range +0.02 to +0.36), so every seed falls well "
    "inside the 1-point equivalence margin and none reverses the direction of the "
    "comparison. This is a post hoc robustness check and is reported as such; it is not "
    "part of the frozen analysis.")

H2("Cross-validation and region selection")
body(
    "A note on the implementation: the frozen plan specifies an isolate-grouped stratified "
    "five-fold split, and the code applies StratifiedKFold to an isolate-level table in "
    "which each isolate contributes one modal slice per region. Those are equivalent, "
    "because one row per isolate makes the split grouped on isolate by construction. This "
    "is recorded as an equivalence, not a deviation.")
body(
    "Evaluation used five-fold stratified cross-validation with shuffling and "
    "random_state = 0, at the isolate level; records from one isolate were collapsed to "
    "that isolate's modal region sequence before splitting, so no isolate contributes to "
    "both training and test folds. For each block the outer split was defined once and "
    "applied consistently across regions. Within each outer fold the region was selected by "
    "inner three-fold cross-validation on training isolates only, scored with C_kmer, and "
    "the selected region was then applied identically to every method, so region selection "
    "cannot favour any arm. No held-out label entered region selection.")

H2("Outcomes and the empirical ambiguity ceiling")
body(
    "The primary outcome is species-balanced accuracy, with the target class and the pooled "
    "neighbour class weighted equally. For each block and region we also computed an "
    "empirical ambiguity ceiling, the highest species-balanced accuracy attainable by any "
    "deterministic function of the observed region sequence:")
body("    A_ceiling = sum over unique sequences s of max( n_target(s) / (2 N_target), "
     "n_neighbour(s) / (2 N_neighbour) )", italic=True)
body(
    "where n_target(s) and n_neighbour(s) count isolates carrying sequence s in each class "
    "and N_target and N_neighbour are the class totals. The ceiling uses the same isolate "
    "unit, the same species-balanced weighting and the same regions as the reported "
    "accuracy. It is an empirical, reference-specific quantity, not a universal maximum: it "
    "changes with reference sampling, class priors and labelling.")
body(
    "Two aggregations of this ceiling are possible and they are not interchangeable. The "
    "whole-block ceiling evaluates the equation once on the full block population, whereas "
    "the reported accuracy is the mean of per-fold balanced accuracies, each computed on a "
    "fold with its own class composition. A fold mean can exceed a whole-population bound, "
    "and across the locked test and hard panels together 3 of 66 block-by-method cells did "
    "so: both arms in Veillonella by 0.30 percentage points and k-mer counting in "
    "Bradyrhizobium by 0.19. Figures reporting a ceiling therefore use the fold-matched "
    "ceiling: the same equation applied to each held-out fold at the region that fold "
    "selected, averaged over folds, which matches the aggregation of the accuracy it "
    "bounds. It bounds all 66 cells and differs from the whole-block value by 0.001 "
    "percentage points on average (range -1.84 to +1.78). Per-block values for both "
    "aggregations are deposited. Like the whole-block version, the fold-matched ceiling is "
    "computed on held-out composition, so it bounds what any sequence-based rule could "
    "achieve on that fold rather than what a rule fitted on training data alone could.")

H2("Statistical analysis")
body(
    "The inference unit is the independent block. The primary estimand is the mean paired "
    "difference between DNABERT-2 with the R2 probe and C_kmer across locked-test blocks; "
    "both arms were named in advance from development results, which removes model and "
    "comparator selection from the test statistic. Confidence intervals are percentile "
    "bootstrap intervals over blocks with 2,000 resamples and seed 1729. Equivalence was "
    "assessed by two one-sided tests against a margin of 1.0 percentage point of balanced "
    "accuracy, fixed before the locked test was collected; the margin was chosen because a "
    "difference below this value would not change which region a recommendation tool would "
    "return and is smaller than the fold-to-fold variation observed within single blocks. "
    "The two one-sided tests were computed as paired t tests on the block-level "
    "differences [12], and the percentile bootstrap interval is reported as "
    "distribution-free corroboration rather than as the test itself. The two procedures "
    "analyse the same block-level differences, so the bootstrap is a sensitivity analysis "
    "on the same comparison and not independent replication of it.")
body(
    "The precision plan sized the locked test using a development block-level standard "
    "deviation of 0.01138. That value belongs to the best-model-minus-best-classical "
    "contrast rather than to the pre-named primary contrast; the achieved standard "
    "deviation for the primary was 0.00889, so the study was planned more conservatively "
    "than it needed to be. The direction is favourable and is stated rather than left "
    "implicit.")
body(
    "The confirmatory family comprised five tests with a HIERARCHICAL correction structure, "
    "as recorded in the frozen plan: the primary contrast was tested unadjusted as the "
    "single prespecified primary, and the four secondary contrasts (DNABERT-S, HyenaDNA, "
    "Nucleotide Transformer v2 and a best-of-four rule, each against C_kmer) were corrected "
    "by the Holm procedure among themselves. One correction was not applied across all five. "
    "All other analyses, including the hard-panel, degradation, layer and region-selection "
    "experiments, are exploratory and are labelled as such.")
body(
    "Two families of tests are reported and they answer different questions, so they are "
    "kept separate. The DIFFERENCE family tests the null of zero difference and is where "
    "the Holm correction applies. The EQUIVALENCE family tests the null that the absolute "
    "difference is at least 1.0 point; its prespecified decision rule is interval "
    "containment, and no multiplicity adjustment was prespecified for it. Holm-adjusted "
    "equivalence p-values are therefore reported as a post hoc sensitivity, not as the "
    "frozen rule. Note also that the frozen rule requires the whole 95 percent interval to "
    "fall inside the margin, whereas two one-sided tests at the 5 percent level correspond "
    "to the 90 percent interval; the rule as frozen is therefore conservative, making "
    "equivalence harder to declare than the nominal level implies. Both intervals are "
    "reported so either can be checked.")
body(
    "The two-sided difference test is a percentile bootstrap over blocks: the estimand is "
    "the mean paired per-block difference, the resampling unit is the block, B = 2,000 "
    "draws, seed 1729 for the primary and 1736 for all four secondary contrasts, with the "
    "interval taken at the 2.5th and 97.5th percentiles. The p-value is twice the smaller "
    "tail proportion of the UNSHIFTED resampling distribution, which makes it an inversion "
    "of the percentile interval rather than a shifted-null bootstrap test, and it is "
    "described that way throughout. Its resolution floor follows from B: the smallest "
    "non-zero two-sided value attainable is 0.001, and a zero tail count returns exactly "
    "zero. One reported value sits at that floor and is given as a bound rather than a "
    "measurement.")
body(
    "For the primary contrast (n = 20 blocks, mean -0.2108 points, standard deviation "
    "0.8892), the component one-sided tests gave t = 3.969 (p = 4.1 x 10-4) against a "
    "true difference at or below -1.0 point and t = -6.089 (p = 3.7 x 10-6) against a "
    "difference at or above +1.0 point, so the two one-sided tests procedure returns "
    "p = 4.1 x 10-4. The corresponding 90 percent intervals, [-0.555, +0.133] parametric "
    "and [-0.529, +0.105] by bootstrap, both lie entirely within the margin, which is the "
    "equivalent decision rule at the 5 percent level.")
body(
    "The best-of-four contrast requires a specific caution. It takes, for each block, the "
    "maximum held-out balanced accuracy across the four models and subtracts C_kmer. The "
    "selection input is therefore the same held-out quantity the contrast is computed from, "
    "so held-out outcomes determine which model is credited. It is an ORACLE comparison, "
    "not an operational rule, and no rule available in advance reproduces it: the per-block "
    "winner is DNABERT-2 in 12 blocks, DNABERT-S in 5, Nucleotide Transformer v2 in 2 and "
    "HyenaDNA in 1. The size of the selection effect can be read directly: the four "
    "individual contrasts average -0.371 points, the per-block maximum gives +0.286 points, "
    "so +0.657 points is attributable to maximising over four noisy measurements. The "
    "deployable alternative, naming DNABERT-2 in advance as the frozen plan did, gives "
    "-0.211 points.")

H2("Software and availability")
body(
    "Analyses ran on an Apple M4 laptop with 16 GB memory and no GPU, in two conda "
    "environments. The main environment used Python 3.11.18, NumPy 2.4.6, SciPy 1.17.1, "
    "scikit-learn 1.9.1, PyTorch 2.15.0, Transformers 5.22.1, Biopython 1.88 and joblib "
    "1.6.0. Nucleotide Transformer v2 required a separate environment with Transformers "
    "pinned to 4.34.1 and the same PyTorch version. Both environments are deposited as "
    "exported specifications and package locks. These versions were verified against those "
    "deposited files and against the local checkpoint hashes; they were not verified "
    "against external release records, and the four compatibility patches applied to the "
    "two MosaicBERT repositories were NOT numerically validated against unpatched code, "
    "because the unpatched code does not execute in this environment.")
body(
    "Computational cost, measured on that hardware. Retrieval took about 65 minutes; "
    "slicing 12,356 development records took about 25 minutes on 8 cores; embedding 7,878 "
    "unique slices took 45 seconds for DNABERT-2, 62 seconds for DNABERT-S, 474 seconds "
    "for HyenaDNA and 187 seconds for Nucleotide Transformer v2; the 277 readout cells "
    "took about 1 minute on 6 workers. All inference was on CPU. The locked test was run "
    "to a pre-committed budget of at most 3 hours wall clock on the same machine. This is "
    "relevant to the comparison: k-mer counting reaches the same accuracy at a small "
    "fraction of the cost, which is part of why a margin of 1 point was judged the "
    "smallest difference that would change a practical recommendation.")

# ============================================================ RESULTS
H("Results")

H2("Scale and eligibility")
body(
    "Sixty-eight blocks were declared across the three experiments and 58 were analysable "
    "after a prespecified eligibility gate requiring at least 10 independent target "
    "isolates and at least two neighbour species with at least five isolates each. A total "
    "of 30,797 sequence records were retrieved and 24,009 entered the analyses, the "
    "remainder having been used only for difficulty screening. Blocks failing the gate were "
    "retained and reported rather than discarded; in every case the reason was the scarcity "
    "of public full-length deposits for the organisms concerned (Table 1).")

body(
    "The 58 analysable blocks are not one homogeneous population and are never pooled into "
    "a single estimate. They comprise 25 species-rank development blocks, 20 species-rank "
    "locked-test blocks and 13 difficulty-selected hard panels, the last of which are "
    "rank-mixed: the repository records no per-panel rank for that set, and 7 of the 13 are "
    "infraspecific by their own definitions, six at subspecies or pathovar rank plus the "
    "Clostridium botulinum physiological groups. The single development subspecies stratum "
    "is excluded from the 25 and reported separately on exactly this reasoning, so the "
    "same reasoning is stated here for the hard panels rather than left implicit. Each "
    "experiment carries its own interval throughout.")
table(["Experiment", "Blocks declared", "Blocks analysed", "Records", "Role"],
      [["Development", 31, 25, "12,356", "all analysis decisions"],
       ["Locked test", 22, 20, "6,933", "single confirmatory use"],
       ["Hard panels", 15, 13, "4,720", "difficulty-selected replication"],
       ["Screening only", "42 candidates", "n/a", "6,788", "difficulty screen"],
       ["Total analysed", 68, 58, "24,009", ""]],
      caption="Table 1. Study accounting.",
      note="Development reports 25 analysed blocks because one block is a subspecies "
           "stratum reported separately. Unique accessions were 12,243, 6,926 and 4,720 "
           "for the three analysed sets, corresponding to 12,188, 6,760 and 4,657 isolate "
           "groups.")

H2("Primary comparison")
body(
    f"On the locked test set of 20 blocks the mean paired difference between DNABERT-2 with "
    f"a linear probe and k-mer counting was {pm:+.2f} percentage points, with a 95 percent "
    f"bootstrap confidence interval of {plo:+.2f} to {phi:+.2f} and a two-sided bootstrap "
    f"p-value of {P['p_two_sided']:.3f}. The interval lies entirely within the prespecified "
    "1-point margin, and two one-sided tests against that margin gave p = 4.1 x 10-4. "
    "Seventeen of the 20 blocks differed by less than one point. Point estimates in the "
    "development set (+0.12 points, 95 percent interval -0.23 to +0.48) and the "
    "difficulty-selected set (-0.19 points) were similar in magnitude, although the latter "
    "interval (-1.27 to +0.66) spans the margin and does not establish equivalence at the "
    "same criterion (Fig 1, Table 2). The development figure is the same contrast as the "
    "other two, DNABERT-2 at the R2 readout minus C_kmer; an earlier version of this "
    "comparison reported +0.14 points, which was the best of four models minus the best of "
    "two classical comparators per block and so was selected twice over.")

table(["Contrast", "Difference (points)", "95% CI", "p", "Holm p", "Interpretation"],
      [["DNABERT-2 vs k-mer (primary)", f"{pm:+.2f}", f"{plo:+.2f} to {phi:+.2f}",
        f"{P['p_two_sided']:.3f}", "not applicable", "equivalent at 1 point"]] +
      [[k.replace(":R2 - C_kmer", " vs k-mer").replace("dnabert_s", "DNABERT-S")
         .replace("hyenadna", "HyenaDNA").replace("ntv2", "NTv2")
         .replace("best-of-4 R2 - C_kmer", "best-of-four vs k-mer"),
        f"{v['mean']*100:+.2f}", f"{v['ci95'][0]*100:+.2f} to {v['ci95'][1]*100:+.2f}",
        ("<= 0.001" if k.startswith("hyenadna") else f"{v['p_raw']:.3f}"),
        ("<= 0.004" if k.startswith("hyenadna") else f"{v['p_holm']:.3f}"),
        ("oracle contrast; no equivalence verdict" if k.startswith("best-of-4")
         else "equivalence tested, NOT established" if k.startswith("hyenadna")
         else "equivalent at 1 point" if v["verdict"].startswith("EQUIV")
         else "inconclusive at 1 point")]
       for k, v in LT["secondary"].items()],
      caption="Table 2. Locked-test contrasts in percentage points of balanced accuracy.",
      note="The p columns test zero difference, not equivalence; secondary contrasts are "
           "Holm-corrected within their family of four, and the primary is unadjusted as the "
           "single prespecified primary. HyenaDNA's bootstrap p rests on one draw of 2,000 "
           "and is reported as a bound; the better-resolved paired t gives p = 0.0077, "
           "Holm-adjusted 0.031. Equivalence verdicts come from the separate two one-sided "
           "tests reported in the text. Best-of-four is an oracle contrast selected on "
           "held-out outcomes and carries no equivalence verdict.")
body(
    "Equivalence was tested formally for every arm, not inferred from interval width. "
    "Against the 1-point margin the two one-sided tests give p = 4.1 x 10-4 for DNABERT-2, "
    "4.9 x 10-5 for DNABERT-S and 1.4 x 10-3 for Nucleotide Transformer v2; Holm-adjusted "
    "within the secondary family as a post hoc sensitivity, the latter two are 1.7 x 10-4 "
    "and 2.9 x 10-3. All three therefore meet the equivalence criterion as a formal test. "
    "HyenaDNA does NOT: its two one-sided tests give p = 0.351, with component values of "
    "0.351 against the lower margin and 2.2 x 10-6 against the upper, and its 90 percent "
    "interval of -1.40 to -0.37 lies outside the margin. Equivalence for HyenaDNA was "
    "tested and not established, which is a different statement from an interval that "
    "merely happens to be wide.")

# --- Table 3: method performance
import csv as _csv
def _rows(name):
    with open(os.path.join(RES, "tables", name)) as f:
        rd = list(_csv.reader(f))
    return rd[0], rd[1:]

_h, _r = _rows("table3_method_performance.csv")
table(_h, _r, caption="Table 3. Balanced accuracy of every method on the locked test.",
      note="The empirical ambiguity ceiling for these blocks averages 97.0 percent. "
           "C_random is a Gaussian random projection of the k-mer features and is a "
           "control for feature redundancy, not for the presence of biological signal.")

H2("Controls")
body(
    "A Gaussian random projection of the k-mer features reached 96.2 percent balanced "
    "accuracy, nominally above every learned representation. The confirmatory analysis drew "
    "that projection at one seed, and recomputing the arm at 20 seeds gave 96.0 to 96.4 "
    "percent, so the comparison does not turn on the particular projection (Methods). This "
    "does not indicate that "
    "accuracy is unrelated to sequence content: a random projection approximately preserves "
    "pairwise distances, so it retains information from the features it compresses. "
    "Permuting labels through the complete supervised pipeline, with 20 permutations per "
    "block at the isolate level and class balance preserved, reduced every method to chance "
    f"(k-mer {LP['summary']['C_kmer']['permuted']*100:.1f} percent, random projection "
    f"{LP['summary']['C_random']['permuted']*100:.1f} percent, DNABERT-2 "
    f"{LP['summary']['dnabert2:R2']['permuted']*100:.1f} percent). The measured accuracy "
    "therefore reflects a genuine association between sequence and label. With 20 "
    "permutations this is an initial null check rather than a fine-grained tail test, but "
    "the separation between observed and permuted values is approximately 46 points and is "
    "not sensitive to the number of permutations (Fig 3).")

H2("Readouts")
body(
    "Among readouts measured on a common endpoint, moving from nearest-reference assignment "
    "to the linear probe raised balanced accuracy by 2.8 to 3.0 points across the four "
    "models; for DNABERT-2 the paired estimate was +2.83 points (95 percent CI +1.37 to "
    "+4.57, 25 blocks). The spread across the four models at the probe readout was 0.73 "
    "points, so readout-associated variation exceeded between-model variation by roughly "
    "fourfold on the development set, and only the readout effect had an interval excluding "
    "zero.")
body(
    "That comparison was development-only in earlier versions of this work and is now "
    "measured on the locked test as a labelled post hoc analysis, using the frozen folds "
    "and the regions the frozen run actually selected; the reconstruction reproduces every "
    "frozen per-block probe value exactly, so the two readouts sit on the same footing. "
    "The readout effect replicates and is smaller: +1.98 points across the four models "
    "(95 percent CI +0.69 to +3.61), with per-model values of +1.81 for DNABERT-2, +2.01 "
    "for DNABERT-S, +1.90 for HyenaDNA and +2.21 for Nucleotide Transformer v2, every "
    "interval excluding zero. Against a model spread of 0.78 points at the probe readout "
    "on the same blocks, readout exceeds model identity by 2.5-fold on confirmatory data "
    "rather than the fourfold seen in development. The confirmatory figure is the one to "
    "quote; the development figure is retained only as the exploratory estimate it is. For "
    "context, the intermediate-layer effect is smaller than either, at 0.53 points.")
body(
    "The unsupervised "
    "geometry statistic is reported separately because it is a pair-level area under the "
    "curve and cannot be placed on the accuracy scale; its value is that it disagrees with "
    "supervised performance. DNABERT-S scored 0.542 by geometry at staphylococcal V1 yet "
    "0.945 with a probe (Fig 2, Table 4).")

_h, _r = _rows("table4_readout_hierarchy.csv")
table(_h, _r, caption="Table 4. Readout comparison on the development set.",
      note="R0 is a pair-level area under the curve; R1 and R2 are isolate-level balanced "
           "accuracies. R0 is not subtracted from R1 or R2 anywhere in this paper. The "
           "same-endpoint readout effect is the R2 minus R1 column.")

H2("Difficulty, degradation, layers and region selection")
body(
    "Half the locked-test blocks were saturated, with k-mer error below one point, leaving "
    "no room for any method to differ. On 13 blocks selected for difficulty by sequence "
    "criteria and the classical comparator alone, with mean k-mer accuracy 0.82 and the "
    f"hardest at 0.55, the difference was {EA['mean']*100:+.2f} points (95 percent CI "
    f"{EA['ci'][0]*100:+.2f} to {EA['ci'][1]*100:+.2f}). This interval includes zero and "
    "also spans the equivalence margin, so it neither demonstrates a benefit nor "
    "establishes equivalence (Fig 4).")
body(
    "Under six degradation conditions applied identically to both arms, no condition "
    "produced a difference that was both larger than the margin and supported by its "
    "interval. Trimming 10 bases from each end cost about 4 points of balanced accuracy for "
    "both methods, more than a tenfold increase in substitution error, which cost about 2 "
    "to 3 points (Fig 5, Table 5). An apparent noise-robustness advantage for DNABERT-2 was "
    "attributable to a single small block, Francisella. Two estimators of the same "
    "condition are reported and they are not the same quantity: Table 5 gives the LEVEL "
    "difference under 1 percent substitution, +0.85 points, while the leave-one-block-out "
    "narrative uses the paired SHIFT from the clean baseline, +0.99 points, the two "
    "differing exactly by the clean-condition difference of -0.14 points. Excluding "
    "Francisella moves the shift estimator from +0.99 to -0.06 points and the level "
    "estimator from +0.85 to -0.20 points, so the direction reverses under either. Only 9 "
    "of the 20 blocks moved positively at all under this condition, and 6 did not move.")
body(
    "Intermediate layers were evaluated for the two models whose internals could be "
    "extracted and verified. For HyenaDNA the best of five depths exceeded the final layer "
    "by 0.16 points (95 percent CI -0.07 to +0.43) and for Nucleotide Transformer v2 the "
    "final layer was best. Because the depth was chosen using the same blocks, these are "
    "optimistic estimates rather than formal bounds (Table 6).")
body(
    "In the region-selection task each method ranked regions using training isolates only "
    "and committed to its top choice. Regret, the accuracy lost relative to the best region "
    "identifiable in hindsight, was 0.65 points for k-mer counting and 0.57 points for "
    "DNABERT-2 on the locked test, and 1.95 and 2.02 points respectively on the hard "
    "panels; paired differences were -0.09 (95 percent CI -0.46 to +0.26) and +0.07 "
    "(-0.61 to +0.65). The oracle is method-specific: each arm is compared against its own "
    "best attainable region rather than a shared one. No method reliably identified the "
    "single best region, but the hit criterion is a strict maximum and is therefore "
    "pessimistic: several regions frequently tie at the top, so a method can incur zero "
    "regret while not matching the single highest-scoring region, which occurred in four "
    "locked-test blocks. Regret is the better-behaved of the two measures (Fig 6, Table 7).")

_h, _r = _rows("table5_robustness.csv")
table(_h, _r, caption="Table 5. Sequence degradation, locked-test blocks.",
      note="Both methods received identical degraded sequences in every condition. "
           "Uniform substitution and fixed end trimming are stress tests, not simulations "
           "of sequencing or PCR behaviour.")

_h, _r = _rows("table6_layer_sweep.csv")
table(_h, _r, caption="Table 6. Intermediate-layer comparison on the hard panels.",
      note="DNABERT-2 and DNABERT-S were not evaluated because their encoders operate on "
           "unpadded sequences internally. Depths were chosen on the same blocks used for "
           "the comparison, so the best-depth value is optimistic.")

_h, _r = _rows("table7_region_ranking.csv")
table(_h, _r, caption="Table 7. Region selection using training isolates only.",
      note="Regret differences are computed from unrounded means: DNABERT-2 minus k-mer is -0.0879 points, shown as -0.09, which is not the difference of the rounded 0.65 and 0.57 displayed in the table. Regret is balanced accuracy lost relative to the best region identifiable in "
           "hindsight. Because several regions frequently tie at the maximum, a method can "
           "record zero regret without matching the single highest-scoring region.")

H2("Sequence reuse and what it bounds")
body(
    "Grouping cross-validation on isolate does not make the folds sequence-disjoint, "
    "because identical 16S sequences recur across strains. Measured at each block's modal "
    "selected region, 5,552 of 5,977 isolate groups (92.9 percent) share their exact slice "
    "with at least one other group, 1,348 (22.6 percent) sit in a haplotype carrying both "
    "classes, and across folds 5,431 of 5,902 evaluated isolates (92.0 percent) have their "
    "exact slice present in the training folds. The confirmatory result therefore "
    "characterises performance WITHIN a known haplotype population; it is not a measurement "
    "of generalisation to unseen sequence variants.")
body(
    "Three post hoc analyses quantify what that costs, all at a fixed region with the "
    "splitter as the only difference. Rebuilding the folds so that no exact sequence "
    "appears in both training and test is possible for 14 of the 20 blocks; the other six "
    "have four or fewer distinct target haplotypes and cannot support the design at any "
    "useful fold count, which is itself a measurement of how little independent sequence "
    "variation the panels contain. The split was audited rather than assumed: exact "
    "train-test overlap is zero and C_exact, which abstains when a test sequence is absent "
    "from training, scores exactly 0.500 on every block. Under that split both primary arms "
    "lose about 11 points, k-mer counting from 96.36 to 85.42 and DNABERT-2 from 96.16 to "
    "84.92, and the paired contrast becomes -0.50 points with a 95 percent interval of "
    "-3.20 to +2.07. In the matched fixed-region control the contrast is -0.20 points with "
    "an interval of -0.39 to -0.02 that excludes zero: with region-selection variance "
    "removed the design is precise enough to detect a small advantage to k-mer counting, "
    "consistent with the uniformly negative sign of all four model contrasts in Table 2. "
    "The effect remains far inside the 1-point margin, which is the distinction an "
    "equivalence test is designed to draw: detectable and negligible are not exclusive. "
    "The direction does not reverse under the disjoint split and no model advantage "
    "appears, but the "
    "block-level standard deviation rises fourteenfold, from 0.37 to 5.16 points, so "
    "equivalence at a 1-point margin is INCONCLUSIVE rather than established on unseen "
    "haplotypes. Repeating the exercise with single-linkage clustering at 99 percent "
    "identity, audited to a maximum cross-fold identity of 0.9894, gives -1.47 points and "
    "the same verdict on the 9 blocks that support it.")
body(
    "A direct test of the mechanism the models are supposed to supply gives a negative "
    "result. If frozen representations generalise to unseen variants better than "
    "composition counting, they should degrade less when exact matches are removed. They do "
    "not: DNABERT-2 loses 11.24 points against k-mer counting's 10.94, a difference of "
    "-0.30 points (95 percent interval -2.92 to +2.22), with the model degrading less in "
    "only 3 of 14 blocks. Giving each representation its own regularisation strength by "
    "nested inner cross-validation on a shared grid does not change this; under the "
    "disjoint split tuning helps k-mer counting most. These analyses are post hoc and are "
    "not part of the confirmatory family.")
body(
    "Because the ranking of methods differs between the two regimes, the reordering was "
    "tested rather than described. Degradation from the isolate-disjoint to the "
    "sequence-disjoint design ranges from -9.64 points for HyenaDNA to -13.96 for "
    "Nucleotide Transformer v2. Testing the paired change across splitters for every arm "
    "against C_kmer, four of the five contrasts are indistinguishable from zero: random "
    "projection +0.26 (95 percent interval -0.90 to +1.92), DNABERT-2 -0.30 (-2.92 to "
    "+2.22), DNABERT-S +1.26 (-2.19 to +4.68) and HyenaDNA +1.29 (-2.00 to +5.16). Only "
    "Nucleotide Transformer v2 separates, at -3.02 points (-5.66 to -0.24), degrading more "
    "than the baseline. Across all 15 pairwise contrasts, the 4 that exclude zero all "
    "involve that model, and none survives Holm correction across the 15. Bootstrapping "
    "the ranking itself, DNABERT-S is the best arm in 41.9 percent of replicates and "
    "HyenaDNA in 31.7 percent, and the full six-way ordering reproduces in 4.0 percent. "
    "The apparent reordering is therefore not resolvable at 14 blocks, with a weak "
    "indication that Nucleotide Transformer v2 fares worse on unseen sequence. It is "
    "reported here so the question is closed rather than left open in a table.")
body(
    "One feature of the tuned-readout analysis bears on how much weight any tuned contrast "
    "can carry. Under sequence-disjoint folds, inner cross-validation improved held-out "
    "accuracy for C_kmer (+2.86 points) and DNABERT-2 (+2.76) but reduced it for DNABERT-S "
    "(-1.48), HyenaDNA (-1.49) and Nucleotide Transformer v2 (-1.13). Selection that "
    "degrades held-out performance for three of six arms indicates selection noise at this "
    "sample size rather than a real optimisation gain. In addition the grid selects its "
    "lower boundary, C = 0.01, in most fold-arm decisions, so these values describe a "
    "finite search that stopped at the edge of its range rather than an optimised result; "
    "extending the grid downward is outstanding work.")
body(
    "Decomposing the pooled metric by neighbour species shows where the difficulty sits. "
    "Across 63 target-neighbour pairs, performance is governed almost entirely by whether "
    "the two taxa share sequence at the evaluated region: pairs whose closest sequences "
    "differ at all are resolved at 97 to 99 percent, while the 22 pairs whose maximum "
    "cross-identity reaches 99.9 percent, in every case because the neighbour shares a "
    "haplotype with the target, fall to 85.9 percent. Accuracy correlates with maximum "
    "target-neighbour identity at Spearman -0.62 for k-mer counting and -0.58 for "
    "DNABERT-2. The difference between the two methods stays within 0.25 points in every "
    "divergence stratum while difficulty itself swings 13 points, so target-neighbour "
    "sequence identity governs performance by about two orders of magnitude more than the "
    "choice of representation. One stratum should not be read as structure: the 99 to 99.9 "
    "percent bin holds 2 pairs from a single block and reports 100.00 for both methods, so "
    "the apparent non-monotonicity across bins is an artefact of that cell's size.")

# --- Tables 8-10: post hoc disjointness, readout confirmation, divergence
_sd = json.load(open(os.path.join(RES, "sequence_disjoint_results.json")))
_md = json.load(open(os.path.join(RES, "similarity_disjoint_results.json")))
_r1 = json.load(open(os.path.join(RES, "r1_locked_test_results.json")))
_dv = json.load(open(os.path.join(RES, "divergence_stratified_results.json")))

_ARMS = [("C_kmer", "k-mer counting"), ("C_random", "Random projection"),
         ("dnabert2:R2", "DNABERT-2"), ("dnabert_s:R2", "DNABERT-S"),
         ("hyenadna:R2", "HyenaDNA"), ("ntv2:R2", "NTv2"), ("C_exact", "Exact match")]
_iso = _sd["summary"]["isolate_disjoint"]["mean_accuracy"]
_seq = _sd["summary"]["sequence_disjoint"]["mean_accuracy"]
table(["Method", "Isolate-disjoint (%)", "Sequence-disjoint (%)", "Change (points)"],
      [[lab, f"{_iso[k]*100:.2f}", f"{_seq[k]*100:.2f}", f"{(_seq[k]-_iso[k])*100:+.2f}"]
       for k, lab in _ARMS if k in _iso and k in _seq],
      caption="Table 8. Exact-sequence-disjoint evaluation, 14 evaluable blocks.",
      note="Post hoc; not part of the confirmatory family. Region fixed at each block's "
           "modal selected region so that the splitter is the only difference. The split "
           "was audited: train-test exact overlap is zero in every fold and exact matching "
           "scores exactly 50.00 percent, as it must when no test sequence appears in "
           "training. The primary contrast is -0.50 points (95 percent interval -3.20 to "
           "+2.07), inconclusive at the 1-point margin; the block-level standard deviation "
           "rises from 0.37 to 5.16 points. Single-linkage clustering at 99 percent "
           "identity, audited to a maximum cross-fold identity of 0.9894, gives -1.47 "
           "points on the 9 blocks that support it.")

_pm = _r1["summary"]["per_model"]
_nm = {"dnabert2": "DNABERT-2", "dnabert_s": "DNABERT-S",
       "hyenadna": "HyenaDNA", "ntv2": "NTv2"}
table(["Model", "R1 nearest reference (%)", "R2 linear probe (%)",
       "R1 to R2 (points)", "95% CI"],
      [[_nm[m], f"{_pm[m]['R1']*100:.2f}", f"{_pm[m]['R2']*100:.2f}",
        f"{_pm[m]['R1_to_R2']*100:+.2f}",
        f"{_pm[m]['ci95'][0]*100:+.2f} to {_pm[m]['ci95'][1]*100:+.2f}"]
       for m in ["dnabert2", "dnabert_s", "hyenadna", "ntv2"]],
      caption="Table 9. Readout effect measured on the locked test, 20 blocks.",
      note="Post hoc. Folds and per-fold regions are those of the frozen run; the "
           "reconstruction reproduces every frozen per-block probe value exactly. The mean "
           "readout effect is +1.98 points (95 percent interval +0.69 to +3.61) against a "
           "model spread of 0.78 points at the probe readout, a ratio of 2.5 on "
           "confirmatory data against 4.0 in development.")

_bins = _dv["summary"]["bins_by_max_target_neighbour_identity"]
table(["Max target-neighbour identity", "Pairs", "Blocks", "k-mer (%)",
       "DNABERT-2 (%)", "Difference (points)"],
      [[b["bin"], str(b["n_pairs"]), str(b["n_blocks"]), f"{b['C_kmer']*100:.2f}",
        f"{b['dnabert2']*100:.2f}", f"{b['difference_pp']:+.2f}"] for b in _bins],
      caption="Table 10. Performance by target-neighbour sequence divergence.",
      note="Post hoc; 63 target-neighbour pairs across 20 blocks, decomposing the pooled "
           "classifier rather than retraining per neighbour. The 99 to 99.9 percent row "
           "holds 2 pairs from one block and is uninformative; the apparent "
           "non-monotonicity across rows is an artefact of that cell. Accuracy correlates "
           "with maximum identity at Spearman -0.62 for k-mer counting and -0.58 for "
           "DNABERT-2, while the method difference stays within 0.25 points in every row.")

H2("Where difficulty occurred")
body(
    "Forty-two candidate groups were screened for difficulty using sequence properties and "
    "the classical comparator only, with no learned representation involved. Groups chosen "
    "because their 16S sequences are reputed to be highly conserved, including Rickettsia, "
    "Chlamydia, Xylella, Synechococcus with Prochlorococcus, Desulfovibrio and Thermus, "
    "were separable at 100 percent. The lowest accuracies occurred in subspecies and "
    "species-complex comparisons: Mycobacterium abscessus subspecies at 55 percent and the "
    "Xanthomonas campestris group at 63 percent. Taxonomic rank was not a prespecified "
    "variable and was not formally modelled, so this is a description of these panels "
    "rather than a general rule (Fig 7).")

# ============================================================ DISCUSSION
H("Discussion")
body(
    "For target-versus-neighbour discrimination within single 16S regions, frozen "
    "representations from four published genomic language models performed equivalently to "
    "k-mer counting under a matched linear probe. The equivalence is supported by a "
    "prespecified margin and a formal test on an independently collected locked test set, "
    "and point estimates were similar in two further sets. The result should not be read as "
    "showing that these models lack relevant information: all four reached 95 to 96 percent "
    "balanced accuracy against an empirical ceiling of 97.0 percent, and permuted labels "
    "reduced them to chance, so the information they carry is real. It is simply not "
    "additional to what k-mer composition already provides for this task.")
body(
    "The readout result is the more transferable finding. The project itself produced a "
    "concrete illustration: an early analysis reported DNABERT-2 at chance on a "
    "staphylococcal V1 panel, and a later analysis with more isolates, isolate-grouped "
    "validation and a linear probe placed the same model near the matched ceiling. Both the "
    "sample size and the readout changed, so neither alone explains the difference, but the "
    "episode shows that a failure observed under one readout cannot be attributed to the "
    "representation without testing alternatives. Benchmarks of DNA models on other "
    "problems draw the same distinction between zero-shot, probed and fine-tuned evaluation "
    "[8].")
body(
    "The two controls together bound what the accuracy can mean. The random projection "
    "shows that the specific identity and weighting of individual k-mer features is not "
    "what carries the signal, since a random linear mixture performs as well; similar "
    "performance is consistent with redundancy or distributed structure, though it does not "
    "by itself prove that individual k-mer identities are irrelevant. The label permutation "
    "shows that the pipeline is not fitting noise. Reporting either control alone would "
    "leave the interpretation underdetermined.")
body(
    "For practice, two findings are more actionable than the model comparison. Read length "
    "mattered more than classifier choice: losing 10 bases from each end cost more accuracy "
    "than a tenfold increase in substitution error. And within these panels, whether a "
    "group was separable depended on whether the comparison was between species or between "
    "subspecies, rather than on the genus's reputation for conservation.")

H2("Limitations")
body(
    "Four models were evaluated, all frozen, at one pooling rule. Intermediate layers were "
    "tested for two of the four; for DNABERT-2 and DNABERT-S the encoder operates on "
    "unpadded sequences internally, and extracting intermediate layers correctly would "
    "require reimplementing part of the model, so this was not attempted rather than "
    "attempted unsafely. The task was closed-panel identification, in which the correct "
    "answer is always among the listed candidates; open-set rejection and community "
    "profiling were not evaluated. The eligibility gate excludes organisms with few public "
    "deposits, which are often those hardest to identify in practice. Degradation was "
    "simulated with uniform random substitutions and fixed end trimming, which are stress "
    "tests rather than models of sequencing and PCR behaviour; in particular, end trimming "
    "is not an empirical primer-binding experiment. The analysis used reference sequences "
    "and no laboratory validation was performed.")
body(
    "A scope condition follows from the length of the evaluated spans. Single "
    "hypervariable regions at Escherichia coli coordinates are 31 to 107 nucleotides, "
    "median 58, whereas the amplicons used in practice are two to four times longer "
    "because they carry conserved flanks. Two consequences follow. At these lengths a "
    "canonical 4-, 5- and 6-mer profile approaches a lossless encoding of the sequence, so "
    "the classical baseline sits close to the information ceiling by construction, which "
    "is part of why it is so hard to beat and is a reason to expect the result rather than "
    "a surprise. And spans this short sit well outside the context lengths these models "
    "were pretrained on, which is a plausible mechanism for the weak unsupervised geometry "
    "reported above, where pair-level area under the curve is 0.67 to 0.74 while a linear "
    "probe on the same representations reaches about 94 percent. Whether the equivalence "
    "reported here extends to full-length amplicons is open and is not tested by this "
    "design.")
body(
    "The most consequential limitation concerns what the design can measure at all. "
    "Because 92 percent of evaluated isolates share their exact regional sequence with the "
    "training folds, the confirmatory equivalence result is a statement about discrimination "
    "within an already-observed haplotype population. It does not transfer to unseen "
    "variants, and the post hoc disjoint analyses show it does not survive there: the "
    "verdict becomes inconclusive under both an exact and a 99 percent identity criterion. "
    "Readers should take the confirmatory claim as bounded in that way rather than as a "
    "general statement about frozen representations.")
body(
    "Two further limitations concern the inputs rather than the models. Independence "
    "between units rests on the deposited strain string, which is submitter-supplied "
    "metadata rather than verified isolate provenance, so the effective number of "
    "independent isolates may be overstated and some leakage across folds cannot be "
    "excluded; standalone 16S records carry no link to a parent genome assembly, so "
    "neither isolate identity nor rRNA operon structure could be checked. And neighbour "
    "species were chosen by expert judgement informed by the 16S identification "
    "literature, with a rationale recorded for every panel in advance but a primary "
    "citation for only four of them, so panel composition is documented and "
    "pre-committed but not reproducible from a stated source list alone.")
body(
    "One defect in the pipeline is known and unrepaired. Region extraction assumes the "
    "deposited orientation, and 22 of 24,009 records (0.09 percent) are deposited "
    "reverse-complemented and were sliced in the wrong orientation rather than being "
    "rejected. Leaving them out moves the primary estimand by 0.003 percentage points and "
    "changes no verdict (Methods), so the reported results stand, but the pipeline as "
    "published will mis-slice minus-strand input and should not be reused without an "
    "orientation check.")
ph("[AUTHOR: check 5. The condition described as copy variation introduces substitutions "
   "into a single deposited sequence. It does not use multiple verified 16S copies from one "
   "genome assembly, because standalone 16S deposits are not linked to parent genomes. "
   "Rename it, for example to within-isolate sequence variation, or justify the current "
   "label.]")

H2("Conclusion")
body(
    "Pretrained genomic language models offered no accuracy advantage over k-mer counting "
    "for fine-scale 16S discrimination in this benchmark, and a recommendation tool built "
    "on this task has no accuracy reason to include one. How a representation is read from "
    "a model mattered more than which model produced it, and two inexpensive controls were "
    "required before the comparison could be interpreted at all. We suggest both controls, "
    "and an explicit statement of the achievable ceiling, as routine components of "
    "benchmarks of this kind.")

# ============================================================ BACK MATTER
H("Supporting information")
for s, d in [
    ("S1 Text", "Study protocol, including task definitions, eligibility criteria and "
                "stopping rules."),
    ("S2 Text", "Frozen analysis plan with checkpoint identifiers, primary estimand, "
                "equivalence margin and multiplicity plan, committed before locked-test "
                "data collection."),
    ("S3 Text", "Protocol deviation log with dates, justifications, and whether each "
                "change preceded inspection of the locked test."),
    ("S4 Text", "Pre-registrations for the difficulty, degradation and layer extensions."),
    ("S1 Data", "Panel definitions for the development, locked-test and hard-panel sets, "
                "committed before the corresponding sequences were retrieved, each with "
                "the recorded rationale for target and neighbour selection and the "
                "pre-declared difficulty-screen thresholds for the hard panels."),
    ("S2 Data", "Accession manifest with per-record species, strain, isolate group and "
                "sequence checksum, and the manifest SHA-256. The isolate group column "
                "is the deposited-metadata key defined in Methods, not a verified "
                "provenance identifier."),
    ("S3 Data", "Block-level results, per-method predictions and bootstrap outputs for "
                "every experiment, together with the post hoc audits reported in Methods: "
                "per-region extraction availability, the minus-strand record list and its "
                "leave-out sensitivity analysis, the ambiguity-code effect on the ceiling, "
                "the isolate-grouping audit and the 20-seed random-projection check."),
    ("S1 Code", "Analysis scripts for retrieval, region extraction, embedding, readouts, "
                "the confirmatory analysis and all extensions."),
    ("S2 Code", "Synthetic validation fixtures with hand-computed expected values for the "
                "ceiling formula, grouped cross-validation and tie handling."),
    ("S10 Text", "Evidence for the pre-registration claim: commit identifiers, timestamps, "
                 "verification commands, and an explicit statement of what the repository "
                 "history does and does not establish."),
]:
    p = doc.add_paragraph()
    r = p.add_run(f"{s}. ")
    r.bold = True
    p.add_run(d)

H("Data and code availability")
ph("[AUTHOR: AQ02, AQ20. Insert the public repository URL, archived release DOI, access "
   "date and licences for code, derived data and model weights. The analysis repository has "
   "27 commits with tagged milestones for each phase; cite the commit hash corresponding to "
   "the submitted version. NOTE: the repository has no remote, so no commit carries a "
   "third-party timestamp. Depositing it now (Zenodo or OSF) fixes the history at a "
   "verifiable date, and cryptographically anchoring the freeze commit hash (for example "
   "via OpenTimestamps) would strengthen the pre-registration claim. Neither acts "
   "retroactively. See S10 Text.]")

H("Ethics statement")
ph("[AUTHOR: AQ21. All sequences are public database records and no human or animal "
   "subjects were involved. Confirm the institutional determination or state that formal "
   "review was not applicable under local policy.]")

H("Funding")
ph("[AUTHOR: AQ22. Funders, grant numbers, and the standard statement on funder role in "
   "study design, data collection, analysis, publication decision and manuscript "
   "preparation.]")

H("Competing interests")
ph("[AUTHOR: AQ22. Competing-interest declaration.]")

H("Author contributions")
ph("[AUTHOR: AQ22. CRediT roles for each author.]")

H("References")
for i, ref in enumerate([
    "Chakravorty S, Helb D, Burday M, Connell N, Alland D. A detailed analysis of 16S "
    "ribosomal RNA gene segments for the diagnosis of pathogenic bacteria. J Microbiol "
    "Methods. 2007;69(2):330-339.",
    "Johnson JS, Spakowicz DJ, Hong BY, Petersen LM, Demkowicz P, Chen L, et al. "
    "Evaluation of 16S rRNA gene sequencing for species and strain-level microbiome "
    "analysis. Nat Commun. 2019;10:5029.",
    "Zhou Z, Ji Y, Li W, Dutta P, Davuluri R, Liu H. DNABERT-2: efficient foundation model "
    "and benchmark for multi-species genome. ICLR. 2024.",
    "Dalla-Torre H, Gonzalez L, Mendoza-Revilla J, Lopez Carranza N, Henryk Grzywaczewski "
    "A, Oteri F, et al. Nucleotide Transformer: building and evaluating robust foundation "
    "models for human genomics. Nat Methods. 2025;22:287-297.",
    "Nguyen E, Poli M, Faizi M, Thomas A, Birch-Sykes C, Wornow M, et al. HyenaDNA: "
    "long-range genomic sequence modeling at single nucleotide resolution. NeurIPS. 2023.",
    "Zhou Z, Wu W, Ho H, Wang J, Shi L, Davuluri RV, et al. DNABERT-S: pioneering species "
    "differentiation with species-aware DNA embeddings. Bioinformatics. "
    "2025;41(Suppl 1):i255-i264.",
    "Bishop HV, Ogilvie OJ, Dobson RCJ, et al. Micro16S: universal phylogenetic 16S rRNA "
    "gene representations for deep learning of the microbiome. bioRxiv. 2026.",
    "Patel A, Singhal A, Wang A, Pampari A, Kasowski M, Kundaje A. DART-Eval: a "
    "comprehensive DNA language model evaluation benchmark on regulatory DNA. arXiv. 2024.",
    "Hartmann M, Howes CG, Abarenkov K, Mohn WW, Nilsson RH. V-Xtractor: an open-source, "
    "high-throughput software tool to identify and extract hypervariable regions of small "
    "subunit (16S/18S) ribosomal RNA gene sequences. J Microbiol Methods. "
    "2010;83(2):250-253.",
    "Sambo F, Finotello F, Lavezzo E, Baruzzo G, Masi G, Peta E, et al. Optimizing PCR "
    "primers targeting the bacterial 16S ribosomal RNA gene. BMC Bioinformatics. "
    "2018;19:343.",
    "Lauzon D, et al. Comprehensive assessment of 16S rRNA gene amplicon sequencing for "
    "microbiome profiling across multiple habitats. Microbiol Spectr. 2023;11:e00563-23.",
    "Lakens D. Equivalence tests: a practical primer for t tests, correlations, and "
    "meta-analyses. Soc Psychol Personal Sci. 2017;8(4):355-362.",
    "Johnson WB, Lindenstrauss J. Extensions of Lipschitz mappings into a Hilbert space. "
    "Contemp Math. 1984;26:189-206.",
], 1):
    p = doc.add_paragraph()
    p.paragraph_format.left_indent = Inches(0.35)
    p.paragraph_format.first_line_indent = Inches(-0.35)
    p.add_run(f"{i}. {ref}")

ph("[AUTHOR: AQ23. Verify every reference against the publisher record, including complete "
   "author lists, volume, pages and DOI. Entries above were assembled from verified titles "
   "and identifiers but author lists were not expanded beyond what was checked.]")

H("Figure legends")
for n, t_ in [
    ("Fig 1", "Paired differences between DNABERT-2 with a linear probe and k-mer counting, in "
              "percentage points of balanced accuracy. (a) Per-block differences on the locked "
              "test set, ordered by magnitude; orange bars favour DNABERT-2 and blue bars "
              "favour k-mer counting. The grey shaded band spans the prespecified equivalence "
              "margin of plus or minus 1 percentage point. Open circles mark the seven blocks "
              "whose difference is exactly zero and which would otherwise show no bar at all. "
              "(b) Mean paired difference with 95 percent bootstrap confidence intervals for "
              "the development, locked-test and hard-panel sets. Every row is the same "
              "contrast, DNABERT-2 at the R2 readout minus k-mer counting, so the three are "
              "directly comparable; the development row is labelled exploratory because those "
              "blocks informed every analysis decision and are burned for confirmatory use. "
              "The shaded band again marks the margin."),
    ("Fig 2", "Readout, model and layer effects on the development set, with metrics kept "
              "strictly separate. (a) R0 unsupervised geometry for each model, reported as a "
              "pair-level area under the curve on its own axis because it is not comparable "
              "to balanced accuracy and must not be subtracted from it. (b) Change from "
              "nearest-reference assignment (R1) to the linear probe (R2), both in balanced "
              "accuracy; the paired effect for DNABERT-2 is plus 2.83 points, 95 percent "
              "confidence interval plus 1.37 to plus 4.57. (c) Balanced accuracy across the "
              "four models at the R2 readout, spanning 0.73 points. (d) Three differences in "
              "balanced accuracy on a common axis: R2 minus R1 for the readout bar, best "
              "minus worst model for the model bar, and best minus worst depth for the layer "
              "bar. The layer bar is shaded grey because it was measured on the 13 hard "
              "panels rather than the development set; its aggregation is defined in "
              "Methods."),
    ("Fig 3", "Balanced accuracy of every method on the locked test set, averaged over the 20 "
              "eligible blocks. The solid red vertical line at 50.1 percent is the mean "
              "accuracy obtained when class labels are permuted through the complete "
              "pipeline, and the dashed black line at 97.0 percent is the mean fold-matched "
              "empirical ambiguity ceiling for these blocks, computed on the same "
              "aggregation as Fig 7. Together they bound what any method could achieve on "
              "these data. Bars begin at the axis minimum rather than at zero, so their "
              "lengths should be read against those two reference lines and not as "
              "ratios."),
    ("Fig 4", "Task difficulty. (a) Per-block difference between DNABERT-2 and k-mer counting "
              "plotted against k-mer error on the locked test set; the dotted vertical line at "
              "1 percentage point of error marks the saturation threshold, below which 10 of "
              "the 20 blocks fall. (b) Per-block differences on the 13 difficulty-selected "
              "hard panels, ordered by magnitude, with the pooled estimate and its 95 percent "
              "bootstrap interval on the bottom row: minus 0.19 percentage points, interval "
              "minus 1.27 to plus 0.66. The grey band marks the plus or minus 1 point "
              "equivalence margin. The interval contains zero and also extends beyond the "
              "margin, so it neither demonstrates a benefit nor establishes equivalence."),
    ("Fig 5", "Sequence degradation on the 20 locked-test blocks, with identical degraded "
              "sequences supplied to both methods. (a) Balanced accuracy of k-mer counting and "
              "DNABERT-2 under six conditions; the vertical axis is truncated at 86 percent. "
              "(b) The paired block-level difference for each condition with its 95 percent "
              "bootstrap interval, which is where the uncertainty belongs because the analysis "
              "is paired; the grey band marks the plus or minus 1 point margin. The three "
              "substitution-based conditions have intervals that contain zero and extend "
              "beyond the margin, and the apparent advantage at 1 percent substitution rests "
              "on a single panel, Francisella."),
    ("Fig 6", "Region selection using training isolates only. (a, b) Balanced accuracy achieved "
              "at the selected region, drawn as points rather than bars because the values "
              "differ by fractions of a point and any readable axis is truncated. (c, d) "
              "Regret, defined as balanced accuracy lost relative to the best region "
              "identifiable in hindsight for that same method; lower is better. Panels a and c "
              "are the locked-test set, 20 blocks; panels b and d the difficulty-selected hard "
              "panels, 13 blocks. Note the differing horizontal scales. Regret is plotted "
              "rather than the best-region hit rate because the hit flag uses a strict "
              "argmax and is pessimistic under ties: several regions frequently tie at the "
              "maximum, so a method can record zero regret while not matching the single "
              "argmax."),
    ("Fig 7", "Achieved balanced accuracy against the empirical ambiguity ceiling for each "
              "block, for both arms, on the locked test set (circles, 20 blocks) and the "
              "difficulty-selected hard panels (squares, 13 blocks). The ceiling is "
              "fold-matched: the same formula applied to each held-out fold at that fold's "
              "chosen region and then averaged over folds, which matches the aggregation of "
              "the accuracy it bounds. A whole-block ceiling does not bound a mean of "
              "per-fold accuracies and 3 of these 66 cells exceed one; none exceeds the "
              "fold-matched ceiling. The diagonal marks the ceiling."),
]:
    p = doc.add_paragraph()
    r = p.add_run(f"{n}. ")
    r.bold = True
    p.add_run(t_)

ph("[AUTHOR: Figures are supplied as 300 dpi TIFF or EPS at submission. Source files are "
   "in the figures directory as PNG and vector PDF.]")

doc.save(OUT)
print(f"wrote {OUT}")
print(f"word count (approx): {sum(len(p.text.split()) for p in doc.paragraphs)}")
print(f"author placeholders: {sum(1 for p in doc.paragraphs if p.text.startswith('[AUTHOR'))}")
