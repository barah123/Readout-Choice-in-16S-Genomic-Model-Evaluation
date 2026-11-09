"""LOCKED TEST BLOCKS — declared and committed before any fetch, slice, embedding or score.

Governed by ../phase3_freeze/FREEZE.md. Single use. Scored once, in Phase 5.

INDEPENDENCE FROM DEVELOPMENT
-----------------------------
Every family below is DISJOINT from the 27 families used in the 31 development blocks.
Protocol §6.3 requires panels sharing species or genomes to be one block; going further and
excluding whole families means no test block shares a family-level ancestor with a
development block, so development-derived choices cannot have been tuned on a close relative.

Development families (excluded here): Staphylococcaceae, Bacillaceae, Streptococcaceae,
Enterococcaceae, Listeriaceae, Clostridiaceae, Lactobacillaceae, Mycobacteriaceae,
Corynebacteriaceae, Bifidobacteriaceae, Enterobacteriaceae, Pseudomonadaceae, Moraxellaceae,
Pasteurellaceae, Vibrionaceae, Legionellaceae, Yersiniaceae, Aeromonadaceae,
Burkholderiaceae, Neisseriaceae, Alcaligenaceae, Campylobacteraceae, Helicobacteraceae,
Bacteroidaceae, Prevotellaceae, Borreliaceae, Leptospiraceae.

SELECTION RULE — identical to blocks_definition.py, applied without reference to any result:
1. Target of established diagnostic/ecological importance where 16S identification against
   close congeners is a recognised practical problem.
2. Neighbours = the congeners most often reported as confusable, or closest relatives by
   published phylogeny. 3-5 per block.
3. One block = one genus (or a tightly-linked genus pair where taxonomy split a clade).
4. Families disjoint from development, as above.
5. Species rank only. No subspecies block in the locked test (the development subspecies
   stratum stays exploratory and is not generalised from).

SIZING: FREEZE.md §5 requires >=14 eligible blocks for the planned precision. Development
lost 16% of declared blocks to the eligibility gate; these taxa are rarer on average, so a
~35% loss is assumed and 22 blocks are declared to land ~14 eligible. If <10 prove eligible
the primary analysis is reported as underpowered and inconclusive.

`expected_class` is again an a-priori hypothesis recorded for calibration only. It is NOT an
input to any analysis, filter, or weighting.
"""

TEST_BLOCKS = [
    # ------------------------------------------------------------ Actinomycetota
    dict(block_id="streptomyces", family="Streptomycetaceae", rank="species",
         target="Streptomyces griseus",
         neighbours=["Streptomyces albus", "Streptomyces coelicolor",
                     "Streptomyces lividans", "Streptomyces venezuelae"],
         expected_class="C",
         rationale="Streptomyces 16S is highly conserved; species assignment normally "
                   "needs MLSA."),

    dict(block_id="nocardia", family="Nocardiaceae", rank="species",
         target="Nocardia farcinica",
         neighbours=["Nocardia asteroides", "Nocardia brasiliensis",
                     "Nocardia nova", "Nocardia cyriacigeorgica"],
         expected_class="B",
         rationale="clinically important; 16S partially resolves the complex."),

    dict(block_id="cutibacterium", family="Propionibacteriaceae", rank="species",
         target="Cutibacterium acnes",
         neighbours=["Cutibacterium granulosum", "Cutibacterium avidum",
                     "Propionibacterium freudenreichii"],
         expected_class="A",
         rationale="skin microbiome relevance; acnes well separated from congeners."),

    dict(block_id="rothia_micrococcaceae", family="Micrococcaceae", rank="species",
         target="Rothia mucilaginosa",
         neighbours=["Rothia dentocariosa", "Micrococcus luteus",
                     "Kocuria rhizophila", "Arthrobacter globiformis"],
         expected_class="A",
         rationale="oral/skin taxa spanning several genera within one family."),

    dict(block_id="actinomyces", family="Actinomycetaceae", rank="species",
         target="Actinomyces israelii",
         neighbours=["Actinomyces naeslundii", "Actinomyces odontolyticus",
                     "Actinomyces viscosus"],
         expected_class="B",
         rationale="oral Actinomyces; partial 16S resolution reported."),

    # ------------------------------------------------------------ Bacillota (new families)
    dict(block_id="paenibacillus", family="Paenibacillaceae", rank="species",
         target="Paenibacillus polymyxa",
         neighbours=["Paenibacillus larvae", "Paenibacillus alvei",
                     "Brevibacillus brevis"],
         expected_class="A",
         rationale="distinct species within a family not used in development."),

    dict(block_id="leuconostoc", family="Leuconostocaceae", rank="species",
         target="Leuconostoc mesenteroides",
         neighbours=["Leuconostoc citreum", "Leuconostoc lactis",
                     "Weissella confusa"],
         expected_class="B",
         rationale="food-fermentation LAB; mesenteroides/citreum are close."),

    dict(block_id="veillonella", family="Veillonellaceae", rank="species",
         target="Veillonella parvula",
         neighbours=["Veillonella dispar", "Veillonella atypica",
                     "Veillonella rogosae"],
         expected_class="C",
         rationale="oral Veillonella spp. are near-identical at 16S."),

    dict(block_id="granulicatella", family="Carnobacteriaceae", rank="species",
         target="Granulicatella adiacens",
         neighbours=["Granulicatella elegans", "Abiotrophia defectiva",
                     "Carnobacterium maltaromaticum"],
         expected_class="B",
         rationale="nutritionally variant streptococci; historically hard to identify."),

    # ------------------------------------------------------------ Pseudomonadota (new families)
    dict(block_id="proteus_morganellaceae", family="Morganellaceae", rank="species",
         target="Proteus mirabilis",
         neighbours=["Proteus vulgaris", "Providencia stuartii",
                     "Morganella morganii", "Providencia rettgeri"],
         expected_class="A",
         rationale="distinct genera within Morganellaceae; good separability expected."),

    dict(block_id="stenotrophomonas", family="Xanthomonadaceae", rank="species",
         target="Stenotrophomonas maltophilia",
         neighbours=["Stenotrophomonas rhizophila", "Xanthomonas campestris",
                     "Xanthomonas oryzae", "Lysobacter enzymogenes"],
         expected_class="B",
         rationale="maltophilia/rhizophila are close; Xanthomonas more distant."),

    dict(block_id="brucella_ochrobactrum", family="Brucellaceae", rank="species",
         target="Brucella melitensis",
         neighbours=["Brucella abortus", "Brucella suis",
                     "Brucella canis", "Ochrobactrum anthropi"],
         expected_class="C",
         rationale="classical Brucella spp. are essentially identical at 16S; "
                   "a strong expected-ambiguous case."),

    dict(block_id="francisella", family="Francisellaceae", rank="species",
         target="Francisella tularensis",
         neighbours=["Francisella novicida", "Francisella philomiragia",
                     "Francisella noatunensis"],
         expected_class="C",
         rationale="tularensis/novicida are 16S-indistinguishable."),

    dict(block_id="bartonella", family="Bartonellaceae", rank="species",
         target="Bartonella henselae",
         neighbours=["Bartonella quintana", "Bartonella bacilliformis",
                     "Bartonella clarridgeiae"],
         expected_class="B",
         rationale="vector-borne pathogens; moderate 16S divergence."),

    dict(block_id="sphingomonas", family="Sphingomonadaceae", rank="species",
         target="Sphingomonas paucimobilis",
         neighbours=["Sphingobium yanoikuyae", "Novosphingobium aromaticivorans",
                     "Sphingopyxis alaskensis"],
         expected_class="A",
         rationale="genus-level contrasts within Sphingomonadaceae."),

    dict(block_id="rhizobium_agrobacterium", family="Rhizobiaceae", rank="species",
         target="Agrobacterium tumefaciens",
         neighbours=["Rhizobium leguminosarum", "Rhizobium etli",
                     "Sinorhizobium meliloti"],
         expected_class="B",
         rationale="Agrobacterium/Rhizobium taxonomy is contested; partial separation."),

    dict(block_id="comamonadaceae", family="Comamonadaceae", rank="species",
         target="Delftia acidovorans",
         neighbours=["Comamonas testosteroni", "Acidovorax avenae",
                     "Variovorax paradoxus"],
         expected_class="A",
         rationale="distinct genera; environmental and occasional clinical isolates."),

    dict(block_id="methylobacterium", family="Methylobacteriaceae", rank="species",
         target="Methylobacterium extorquens",
         neighbours=["Methylobacterium radiotolerans", "Methylobacterium mesophilicum",
                     "Methylobacterium populi"],
         expected_class="C",
         rationale="Methylobacterium 16S is conserved within the genus."),

    dict(block_id="acetobacteraceae", family="Acetobacteraceae", rank="species",
         target="Acetobacter aceti",
         neighbours=["Gluconobacter oxydans", "Komagataeibacter xylinus",
                     "Acetobacter pasteurianus"],
         expected_class="B",
         rationale="acetic acid bacteria; genus-level contrasts expected to resolve."),

    # ------------------------------------------------------------ Bacteroidota (new families)
    dict(block_id="porphyromonas", family="Porphyromonadaceae", rank="species",
         target="Porphyromonas gingivalis",
         neighbours=["Porphyromonas asaccharolytica", "Porphyromonas endodontalis",
                     "Porphyromonas somerae"],
         expected_class="A",
         rationale="periodontal pathogen; well separated from congeners."),

    dict(block_id="flavobacteriaceae", family="Flavobacteriaceae", rank="species",
         target="Elizabethkingia meningoseptica",
         neighbours=["Chryseobacterium gleum", "Chryseobacterium indologenes",
                     "Flavobacterium psychrophilum"],
         expected_class="A",
         rationale="distinct genera within Flavobacteriaceae."),

    # ------------------------------------------------------------ Fusobacteriota
    dict(block_id="fusobacterium", family="Fusobacteriaceae", rank="species",
         target="Fusobacterium nucleatum",
         neighbours=["Fusobacterium necrophorum", "Fusobacterium varium",
                     "Fusobacterium mortiferum"],
         expected_class="B",
         rationale="oral/gut anaerobe; nucleatum subspecies complicate the genus."),
]


def all_species():
    out = []
    for b in TEST_BLOCKS:
        out.append((b["block_id"], b["target"], True))
        for n in b["neighbours"]:
            out.append((b["block_id"], n, False))
    return out


def summary():
    fams = {b["family"] for b in TEST_BLOCKS}
    dev_fams = {
        "Staphylococcaceae", "Bacillaceae", "Streptococcaceae", "Enterococcaceae",
        "Listeriaceae", "Clostridiaceae", "Lactobacillaceae", "Mycobacteriaceae",
        "Corynebacteriaceae", "Bifidobacteriaceae", "Enterobacteriaceae",
        "Pseudomonadaceae", "Moraxellaceae", "Pasteurellaceae", "Vibrionaceae",
        "Legionellaceae", "Yersiniaceae", "Aeromonadaceae", "Burkholderiaceae",
        "Neisseriaceae", "Alcaligenaceae", "Campylobacteraceae", "Helicobacteraceae",
        "Bacteroidaceae", "Prevotellaceae", "Borreliaceae", "Leptospiraceae"}
    return dict(
        n_blocks=len(TEST_BLOCKS),
        n_families=len(fams),
        overlap_with_development_families=sorted(fams & dev_fams),
        n_species_queries=len(all_species()),
        expected_class_counts={
            c: sum(1 for b in TEST_BLOCKS if b["expected_class"] == c) for c in "ABC"},
    )


if __name__ == "__main__":
    import json
    print(json.dumps(summary(), indent=2))
