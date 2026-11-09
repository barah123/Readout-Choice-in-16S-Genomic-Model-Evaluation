"""Taxonomic block definitions for the Protocol 3 benchmark (protocol §6.2, §6.3).

DECLARED BEFORE ANY MODEL RESULT ON THESE PANELS IS COMPUTED. This file is committed prior
to fetching so the panel list cannot be retrofitted to whatever makes embeddings look good or
bad. Protocol §6.2: "Select panels using source diversity, taxonomy, and sequence criteria,
not whichever panels make embeddings look weak. A positive result must be as publishable as a
negative result."

SELECTION RULE (fixed, applied uniformly)
-----------------------------------------
1. Target species: a species of established diagnostic or ecological importance for which
   16S-based identification against close congeners is a recognised practical problem.
2. Neighbour species: the congeners most frequently reported as confusable with the target in
   routine 16S identification, or its closest relatives by published phylogeny. 3-6 per block.
3. A block is one GENUS (or one tightly-linked genus pair where taxonomy has split a clade,
   e.g. Escherichia/Shigella, Clostridium/Clostridioides). Protocol §6.3: panels sharing
   species or genomes must be connected into a single block. Genus-level blocking guarantees
   no species appears in two blocks.
4. Blocks are additionally tagged with `family` so a STRINGENT sensitivity analysis can group
   related genera (the Enterobacteriaceae blocks in particular) and report a reduced,
   conservative block count alongside the primary genus-level count.
5. Species-level and subspecies-level questions are SEPARATE STRATA (protocol §3). The single
   subspecies block below is flagged `rank="subspecies"` and must never be pooled with the
   species-level blocks in a headline estimate.
6. Eligibility (applied after fetching, never before): a block is CONFIRMATORY-ELIGIBLE only
   if the target has >=10 independent isolates and >=2 neighbour species each have >=5, in at
   least one region. Blocks that fail are retained and reported as DESCRIPTIVE-ONLY
   (protocol §6.2's "underpowered_panel_policy: descriptive_only"), never silently dropped.

EXPECTED PANEL CLASS
--------------------
`expected_class` records the a-priori expectation (A informative / B partial / C highly
ambiguous / D low-reference) with a one-line rationale, so that the measured ambiguity can be
compared against a prediction made in advance. Protocol §3 is explicit that a group must NOT
be declared null before measuring its reference-specific ambiguity - so `expected_class` is a
hypothesis recorded for calibration, NOT an input to any analysis or filtering step.
"""

# (block_id, family, rank, target, [neighbours], expected_class, rationale)
BLOCKS = [
    # ---------------------------------------------------------------- Bacillota (Firmicutes)
    dict(block_id="staphylococcus", family="Staphylococcaceae", rank="species",
         target="Staphylococcus aureus",
         neighbours=["Staphylococcus epidermidis", "Staphylococcus haemolyticus",
                     "Staphylococcus saprophyticus", "Staphylococcus hominis",
                     "Staphylococcus lugdunensis", "Staphylococcus capitis"],
         expected_class="A",
         rationale="V1 reported to separate S. aureus from CoNS (Chakravorty 2007)."),

    dict(block_id="bacillus_cereus_group", family="Bacillaceae", rank="species",
         target="Bacillus anthracis",
         neighbours=["Bacillus cereus", "Bacillus thuringiensis",
                     "Bacillus weihenstephanensis", "Bacillus mycoides"],
         expected_class="C",
         rationale="anthracis is a monophyletic subclade of cereus sensu lato; "
                   "near-identical 16S expected across much of the gene."),

    dict(block_id="streptococcus_mitis_group", family="Streptococcaceae", rank="species",
         target="Streptococcus pneumoniae",
         neighbours=["Streptococcus mitis", "Streptococcus oralis",
                     "Streptococcus pseudopneumoniae", "Streptococcus sanguinis"],
         expected_class="C",
         rationale="pneumoniae/mitis/oralis exceed 99% 16S identity; a classic "
                   "16S-unresolvable group."),

    dict(block_id="enterococcus", family="Enterococcaceae", rank="species",
         target="Enterococcus faecalis",
         neighbours=["Enterococcus faecium", "Enterococcus durans",
                     "Enterococcus hirae", "Enterococcus gallinarum"],
         expected_class="A",
         rationale="faecalis is well separated from the faecium/durans/hirae cluster."),

    dict(block_id="listeria", family="Listeriaceae", rank="species",
         target="Listeria monocytogenes",
         neighbours=["Listeria innocua", "Listeria welshimeri",
                     "Listeria seeligeri", "Listeria ivanovii"],
         expected_class="C",
         rationale="Listeria sensu stricto is highly conserved at 16S; "
                   "monocytogenes/innocua notoriously close."),

    dict(block_id="clostridium", family="Clostridiaceae", rank="species",
         target="Clostridioides difficile",
         neighbours=["Clostridium perfringens", "Clostridium butyricum",
                     "Clostridium sporogenes", "Clostridium septicum"],
         expected_class="A",
         rationale="C. difficile is phylogenetically distant from Clostridium "
                   "sensu stricto (hence its reassignment)."),

    dict(block_id="lactobacillus_vaginal", family="Lactobacillaceae", rank="species",
         target="Lactobacillus crispatus",
         neighbours=["Lactobacillus iners", "Lactobacillus gasseri",
                     "Lactobacillus jensenii", "Lactobacillus johnsonii"],
         expected_class="B",
         rationale="vaginal Lactobacillus group; crispatus/gasseri/johnsonii are close "
                   "but iners is distinct. Ties to the microbiome-relevant Panel B."),

    # ---------------------------------------------------------- Actinomycetota
    dict(block_id="mycobacterium", family="Mycobacteriaceae", rank="species",
         target="Mycobacterium tuberculosis",
         neighbours=["Mycobacterium avium", "Mycobacterium intracellulare",
                     "Mycobacterium kansasii", "Mycobacterium abscessus"],
         expected_class="B",
         rationale="V2 reported best for Mycobacterium (Chakravorty 2007); MTBC internally "
                   "clonal but distinct from NTM."),

    dict(block_id="corynebacterium", family="Corynebacteriaceae", rank="species",
         target="Corynebacterium diphtheriae",
         neighbours=["Corynebacterium ulcerans", "Corynebacterium pseudotuberculosis",
                     "Corynebacterium striatum", "Corynebacterium amycolatum"],
         expected_class="B",
         rationale="diphtheriae/ulcerans/pseudotuberculosis form a close cluster; "
                   "striatum/amycolatum more distant."),

    dict(block_id="bifidobacterium_species", family="Bifidobacteriaceae", rank="species",
         target="Bifidobacterium longum",
         neighbours=["Bifidobacterium breve", "Bifidobacterium bifidum",
                     "Bifidobacterium adolescentis", "Bifidobacterium animalis"],
         expected_class="A",
         rationale="distinct Bifidobacterium species; infant-gut relevance (links to the "
                   "systematic review motivating this project)."),

    dict(block_id="bifidobacterium_longum_subsp", family="Bifidobacteriaceae",
         rank="subspecies",
         target="Bifidobacterium longum subsp. infantis",
         neighbours=["Bifidobacterium longum subsp. longum",
                     "Bifidobacterium longum subsp. suis"],
         expected_class="C",
         rationale="SUBSPECIES STRATUM - do not pool with species-level blocks. "
                   "Within-species 16S variation is minimal by definition."),

    # ---------------------------------------------- Pseudomonadota: Enterobacteriaceae
    dict(block_id="escherichia_shigella", family="Enterobacteriaceae", rank="species",
         target="Escherichia coli",
         neighbours=["Shigella flexneri", "Shigella sonnei",
                     "Shigella dysenteriae", "Shigella boydii"],
         expected_class="C",
         rationale="Shigella is genomically E. coli; the canonical 16S-unresolvable pair."),

    dict(block_id="klebsiella", family="Enterobacteriaceae", rank="species",
         target="Klebsiella pneumoniae",
         neighbours=["Klebsiella oxytoca", "Klebsiella variicola",
                     "Klebsiella aerogenes", "Klebsiella quasipneumoniae"],
         expected_class="C",
         rationale="pneumoniae/variicola/quasipneumoniae are near-identical at 16S."),

    dict(block_id="enterobacter", family="Enterobacteriaceae", rank="species",
         target="Enterobacter cloacae",
         neighbours=["Enterobacter hormaechei", "Enterobacter asburiae",
                     "Enterobacter kobei", "Cronobacter sakazakii"],
         expected_class="C",
         rationale="the E. cloacae complex is explicitly defined by 16S "
                   "indistinguishability."),

    dict(block_id="salmonella_citrobacter", family="Enterobacteriaceae", rank="species",
         target="Salmonella enterica",
         neighbours=["Citrobacter freundii", "Citrobacter koseri",
                     "Citrobacter braakii"],
         expected_class="B",
         rationale="Salmonella/Citrobacter are close but generally separable; "
                   "partial separation expected."),

    # -------------------------------------------- Pseudomonadota: other Gammaproteobacteria
    dict(block_id="pseudomonas", family="Pseudomonadaceae", rank="species",
         target="Pseudomonas aeruginosa",
         neighbours=["Pseudomonas fluorescens", "Pseudomonas putida",
                     "Pseudomonas stutzeri", "Pseudomonas syringae"],
         expected_class="A",
         rationale="aeruginosa is well separated from the fluorescens/putida complexes."),

    dict(block_id="acinetobacter", family="Moraxellaceae", rank="species",
         target="Acinetobacter baumannii",
         neighbours=["Acinetobacter calcoaceticus", "Acinetobacter pittii",
                     "Acinetobacter nosocomialis", "Acinetobacter lwoffii"],
         expected_class="C",
         rationale="the Acb complex (baumannii/calcoaceticus/pittii/nosocomialis) is "
                   "defined by 16S near-identity."),

    dict(block_id="haemophilus", family="Pasteurellaceae", rank="species",
         target="Haemophilus influenzae",
         neighbours=["Haemophilus parainfluenzae", "Haemophilus haemolyticus",
                     "Haemophilus ducreyi", "Haemophilus parahaemolyticus"],
         expected_class="B",
         rationale="V3 reported best for Haemophilus (Chakravorty 2007); "
                   "influenzae/haemolyticus are close."),

    dict(block_id="vibrio", family="Vibrionaceae", rank="species",
         target="Vibrio cholerae",
         neighbours=["Vibrio parahaemolyticus", "Vibrio vulnificus",
                     "Vibrio alginolyticus", "Vibrio mimicus"],
         expected_class="B",
         rationale="cholerae/mimicus are very close; the others more distant."),

    dict(block_id="legionella", family="Legionellaceae", rank="species",
         target="Legionella pneumophila",
         neighbours=["Legionella longbeachae", "Legionella micdadei",
                     "Legionella bozemanae", "Legionella dumoffii"],
         expected_class="A",
         rationale="Legionella species are comparatively well resolved at 16S."),

    dict(block_id="yersinia", family="Yersiniaceae", rank="species",
         target="Yersinia pestis",
         neighbours=["Yersinia pseudotuberculosis", "Yersinia enterocolitica",
                     "Yersinia frederiksenii", "Yersinia intermedia"],
         expected_class="C",
         rationale="Y. pestis is a clone of Y. pseudotuberculosis; "
                   "16S-indistinguishable from it (Chakravorty 2007 Panel A)."),

    dict(block_id="aeromonas", family="Aeromonadaceae", rank="species",
         target="Aeromonas hydrophila",
         neighbours=["Aeromonas caviae", "Aeromonas veronii",
                     "Aeromonas salmonicida", "Aeromonas media"],
         expected_class="C",
         rationale="Aeromonas 16S is highly conserved; species assignment normally "
                   "requires housekeeping genes."),

    dict(block_id="burkholderia", family="Burkholderiaceae", rank="species",
         target="Burkholderia pseudomallei",
         neighbours=["Burkholderia thailandensis", "Burkholderia cepacia",
                     "Burkholderia multivorans", "Burkholderia cenocepacia"],
         expected_class="B",
         rationale="pseudomallei/thailandensis are very close; the cepacia complex "
                   "is a separate cluster."),

    # ------------------------------------------------- Pseudomonadota: Betaproteobacteria
    dict(block_id="neisseria", family="Neisseriaceae", rank="species",
         target="Neisseria meningitidis",
         neighbours=["Neisseria gonorrhoeae", "Neisseria lactamica",
                     "Neisseria sicca", "Neisseria mucosa"],
         expected_class="B",
         rationale="meningitidis/gonorrhoeae/lactamica are close; "
                   "sicca/mucosa commensals more distant."),

    dict(block_id="bordetella", family="Alcaligenaceae", rank="species",
         target="Bordetella pertussis",
         neighbours=["Bordetella parapertussis", "Bordetella bronchiseptica",
                     "Bordetella holmesii", "Bordetella avium"],
         expected_class="C",
         rationale="pertussis/parapertussis/bronchiseptica are a single 16S-identical "
                   "clade."),

    # ------------------------------------------------- Campylobacterota
    dict(block_id="campylobacter", family="Campylobacteraceae", rank="species",
         target="Campylobacter jejuni",
         neighbours=["Campylobacter coli", "Campylobacter lari",
                     "Campylobacter fetus", "Campylobacter upsaliensis"],
         expected_class="B",
         rationale="jejuni/coli are close but usually separable at 16S."),

    dict(block_id="helicobacter", family="Helicobacteraceae", rank="species",
         target="Helicobacter pylori",
         neighbours=["Helicobacter hepaticus", "Helicobacter bilis",
                     "Helicobacter cinaedi", "Helicobacter felis"],
         expected_class="A",
         rationale="H. pylori is well separated from enterohepatic Helicobacter spp."),

    # ------------------------------------------------- Bacteroidota
    dict(block_id="bacteroides", family="Bacteroidaceae", rank="species",
         target="Bacteroides fragilis",
         neighbours=["Bacteroides thetaiotaomicron", "Bacteroides ovatus",
                     "Bacteroides uniformis", "Bacteroides vulgatus"],
         expected_class="A",
         rationale="distinct gut Bacteroides species; microbiome relevance."),

    dict(block_id="prevotella", family="Prevotellaceae", rank="species",
         target="Prevotella melaninogenica",
         neighbours=["Prevotella intermedia", "Prevotella bivia",
                     "Prevotella denticola", "Prevotella nigrescens"],
         expected_class="B",
         rationale="Prevotella is diverse; intermedia/nigrescens are a close pair."),

    # ------------------------------------------------- Spirochaetota
    dict(block_id="borrelia", family="Borreliaceae", rank="species",
         target="Borrelia burgdorferi",
         neighbours=["Borrelia afzelii", "Borrelia garinii",
                     "Borrelia hermsii", "Borrelia miyamotoi"],
         expected_class="C",
         rationale="the B. burgdorferi sensu lato complex is 16S-conserved; "
                   "relapsing-fever spp. more distant."),

    dict(block_id="leptospira", family="Leptospiraceae", rank="species",
         target="Leptospira interrogans",
         neighbours=["Leptospira borgpetersenii", "Leptospira biflexa",
                     "Leptospira kirschneri", "Leptospira santarosai"],
         expected_class="B",
         rationale="pathogenic Leptospira cluster tightly; biflexa (saprophytic) "
                   "is an outgroup-like contrast."),
]


def all_species():
    """Every species string needed, deduplicated, preserving block association."""
    out = []
    for b in BLOCKS:
        out.append((b["block_id"], b["target"], True))
        for n in b["neighbours"]:
            out.append((b["block_id"], n, False))
    return out


def summary():
    fams = {}
    for b in BLOCKS:
        fams.setdefault(b["family"], []).append(b["block_id"])
    sp = all_species()
    return dict(
        n_blocks=len(BLOCKS),
        n_blocks_species_rank=sum(1 for b in BLOCKS if b["rank"] == "species"),
        n_blocks_subspecies_rank=sum(1 for b in BLOCKS if b["rank"] == "subspecies"),
        n_families=len(fams),
        families_with_multiple_blocks={k: v for k, v in fams.items() if len(v) > 1},
        n_species_queries=len(sp),
        expected_class_counts={
            c: sum(1 for b in BLOCKS if b["expected_class"] == c) for c in "ABCD"},
    )


if __name__ == "__main__":
    import json
    print(json.dumps(summary(), indent=2))
