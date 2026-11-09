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

ph("[AUTHOR: Author names, ORCID iDs, numbered affiliations, and corresponding-author "
   "email. PLOS requires these on the title page. AQ22.]")
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
    "tests against a prespecified 1-point margin, p = 0.0001). Point estimates were "
    "similar in a development set and in a set of groups selected specifically for "
    "difficulty, although the latter interval was too wide to establish equivalence at the "
    "same margin. Two controls bound the interpretation: a Gaussian random projection of "
    "k-mer counts matched every learned representation, whereas permuting labels through "
    "the complete pipeline reduced every method to chance. Among classification readouts "
    "measured on a common endpoint, moving from nearest-reference assignment to a linear "
    "probe changed balanced accuracy by about 2.9 points, roughly four times the 0.73-point "
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

H2("Isolate grouping")
body(
    "The unit of independence is the isolate group, defined as the pair (species, strain) "
    "where esummary reported a usable strain field, and otherwise as the individual "
    "accession. Bare culture-collection prefixes without a catalogue number (for example "
    "\"ATCC\" alone) were treated as missing rather than as a shared group, because merging "
    "them would combine unrelated deposits. Strain metadata was available for essentially "
    "all records (794 of 796 in an early staphylococcal panel; 469 of 469 in a Bacillus "
    "panel).")
ph("[AUTHOR: AQ04 and check 6. The isolate group is derived from the NCBI strain string. "
   "This is deposited metadata, not verified independent provenance, and standalone 16S "
   "deposits are not linked to a parent genome assembly. State this explicitly as a "
   "limitation and confirm whether any normalization of collection identifiers was "
   "performed beyond the prefix rule described above.]")

H2("Region extraction")
body(
    "Each sequence was aligned globally to the Escherichia coli J01695 reference using the "
    "Biopython PairwiseAligner (match +2, mismatch -1, gap open -5, gap extend -2), and "
    "Escherichia coli coordinates for V1 to V9 were transferred through the alignment. An "
    "extracted region was retained only if its length fell between 0.75 and 1.30 times the "
    "expected length for that region; otherwise it was recorded as unavailable for that "
    "sequence, and these failures are reported rather than dropped silently.")
ph("[AUTHOR: AQ06. Alignment-based transfer was used rather than a hidden Markov model "
   "tool such as V-Xtractor. State this choice and its justification, and confirm strand "
   "handling and the treatment of ambiguous bases within retained regions.]")

H2("Models and representations")
body(
    "Four checkpoints were used, frozen, with no fine-tuning. DNABERT-2 (zhihan1996/"
    "DNABERT-2-117M; weight file SHA-256 beginning 7ff39ec77a48, 468,354,983 bytes) and "
    "DNABERT-S (zhihan1996/DNABERT-S; SHA-256 beginning f3cfc3d05418, 468,320,053 bytes) "
    "were loaded from local patched copies; the patches address CPU and library "
    "compatibility only and make no numerical change. HyenaDNA (LongSafari/"
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
    "the number of k-mer features, the projection matrix generated from a fixed seed and "
    "applied identically to training and test folds.")
ph("[AUTHOR: AQ24 and check 3. C_random was evaluated at a single projection seed. If the "
   "manuscript interprets invariance to projection, repeat across seeds; otherwise state "
   "that only one seed was used.]")

H2("Cross-validation and region selection")
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
    "One residual mismatch should be stated. The ceiling is computed on the full block "
    "population, whereas the reported accuracy is the mean of per-fold balanced accuracies, "
    "each computed on a fold with its own class composition. The two are therefore close "
    "but not formally comparable, and a small exceedance is possible. Across the locked "
    "test this occurred in 2 of 40 block-by-method cells, both in the Veillonella block, "
    "with a maximum excess of 0.30 percentage points.")

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
    "return and is smaller than the fold-to-fold variation observed within single blocks.")
body(
    "The confirmatory family comprised five tests: the primary contrast, and four secondary "
    "contrasts (DNABERT-S, HyenaDNA, Nucleotide Transformer v2 and a best-of-four rule, "
    "each against C_kmer) corrected by the Holm procedure within that family. All other "
    "analyses, including the hard-panel, degradation, layer and region-selection "
    "experiments, are exploratory and are labelled as such.")
ph("[AUTHOR: AQ15 and check 2. The TOST reported here is a bootstrap implementation on the "
   "same paired block-level estimand as the confidence interval, not a t-based TOST; its "
   "component one-sided p-values are 0.0001 and <0.00005. Confirm this description and "
   "deposit the script and raw output, and state why bootstrap and t intervals differ "
   "slightly.]")

H2("Software and availability")
body(
    "Analyses ran on an Apple M4 laptop with 16 GB memory and no GPU, in two conda "
    "environments. The main environment used Python 3.11.18, NumPy 2.4.6, SciPy 1.17.1, "
    "scikit-learn 1.9.1, PyTorch 2.15.0, Transformers 5.22.1, Biopython 1.88 and joblib "
    "1.6.0. Nucleotide Transformer v2 required a separate environment with Transformers "
    "pinned to 4.34.1 and the same PyTorch version.")

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
    "1-point margin, and two one-sided tests against that margin gave p = 0.0001. "
    "Seventeen of the 20 blocks differed by less than one point. Point estimates in the "
    "development set (+0.14 points) and the difficulty-selected set (-0.19 points) were "
    "similar in magnitude, although the latter interval (-1.28 to +0.66) spans the margin "
    "and does not establish equivalence at the same criterion (Fig 1, Table 2).")

table(["Contrast", "Difference (points)", "95% CI", "p", "Holm p", "Interpretation"],
      [["DNABERT-2 vs k-mer (primary)", f"{pm:+.2f}", f"{plo:+.2f} to {phi:+.2f}",
        f"{P['p_two_sided']:.3f}", "not applicable", "equivalent at 1 point"]] +
      [[k.replace(":R2 - C_kmer", " vs k-mer").replace("dnabert_s", "DNABERT-S")
         .replace("hyenadna", "HyenaDNA").replace("ntv2", "NTv2")
         .replace("best-of-4 R2 - C_kmer", "best-of-four vs k-mer"),
        f"{v['mean']*100:+.2f}", f"{v['ci95'][0]*100:+.2f} to {v['ci95'][1]*100:+.2f}",
        f"{v['p_raw']:.3f}", f"{v['p_holm']:.3f}",
        ("equivalent at 1 point" if v["verdict"].startswith("EQUIV")
         else "inconclusive at 1 point")]
       for k, v in LT["secondary"].items()],
      caption="Table 2. Locked-test contrasts in percentage points of balanced accuracy.",
      note="Secondary contrasts are Holm-corrected within their family of four. HyenaDNA is "
           "reliably below k-mer counting but by less than the practical margin.")

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
    "accuracy, nominally above every learned representation. This does not indicate that "
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
    "fourfold, and only the readout effect had an interval excluding zero. The unsupervised "
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
    "attributable to a single small block; excluding it moved the estimate from +0.99 to "
    "-0.06 points.")
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
    "(-0.61 to +0.65). No method reliably identified the single best region, but because "
    "several regions frequently tie at the maximum, a method can incur zero regret while "
    "not matching the single highest-scoring region (Fig 6, Table 7).")

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
      note="Regret is balanced accuracy lost relative to the best region identifiable in "
           "hindsight. Because several regions frequently tie at the maximum, a method can "
           "record zero regret without matching the single highest-scoring region.")

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
                "committed before the corresponding sequences were retrieved."),
    ("S2 Data", "Accession manifest with per-record species, strain, isolate group and "
                "sequence checksum, and the manifest SHA-256."),
    ("S3 Data", "Block-level results, per-method predictions and bootstrap outputs for "
                "every experiment."),
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
    ("Fig 1", "Paired differences between DNABERT-2 with a linear probe and k-mer counting. "
              "(a) Per-block differences on the locked test, ordered by magnitude; the "
              "shaded band is the prespecified 1-point equivalence margin. (b) Mean "
              "differences with 95 percent bootstrap intervals for the development, locked "
              "test and hard-panel sets."),
    ("Fig 2", "Readout, model and layer effects, with metrics kept separate. (a) R0 "
              "unsupervised geometry, reported as a pair-level area under the curve on its "
              "own axis because it is not comparable to accuracy. (b) Change from "
              "nearest-reference assignment to the linear probe, both balanced accuracy. "
              "(c) Balanced accuracy across models at the probe readout. (d) The two "
              "same-metric spreads on one axis; the layer bar is shaded differently because "
              "it was measured on the hard panels."),
    ("Fig 3", "Balanced accuracy of all methods on the locked test, bounded below by the "
              "permuted-label value and above by the empirical ambiguity ceiling."),
    ("Fig 4", "Task difficulty. (a) Per-block difference against k-mer error on the locked "
              "test, showing that half the blocks are saturated. (b) Per-block differences "
              "on the 13 difficulty-selected blocks."),
    ("Fig 5", "Balanced accuracy of k-mer counting and DNABERT-2 under six sequence "
              "degradation conditions, with identical degraded sequences supplied to both."),
    ("Fig 6", "Regret in the region-selection task, defined as balanced accuracy lost "
              "relative to the best region identifiable in hindsight."),
    ("Fig 7", "Achieved balanced accuracy against the empirical ambiguity ceiling for each "
              "block; the diagonal marks the ceiling."),
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
