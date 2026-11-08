# Implementation contracts, example configuration, and audit checklists

Version 1.0 | All paths/interfaces below are proposed, not pre-existing software.

## 1. Recommended repository organization

```text
regioneer/
  README.md
  LICENSE
  CITATION.cff
  pyproject.toml
  environment-lock.yml
  config/
    study.yaml
    models.yaml
    assays.yaml
  manifests/
    genomes.tsv
    taxonomy_crosswalk.tsv
    panels.tsv
    primers.tsv
    exclusions.tsv
    splits.tsv
  workflow/
    Snakefile
    rules/
  src/regioneer/
    references.py
    taxonomy.py
    extraction.py
    pcr.py
    observations.py
    features.py
    embeddings.py
    collisions.py
    splits.py
    scoring.py
    recommendation.py
    training.py
    reporting.py
  tests/
    fixtures/
    test_orientation.py
    test_primer_matching.py
    test_collisions.py
    test_splits.py
    test_no_query_fit.py
  analysis/
    01_qc.R
    02_primary_comparisons.R
    03_robustness.R
    report.qmd
  results/
    qc/
    predictions/
    metrics/
    figures/
    runtime/
  docs/
    model_card.md
    data_card.md
    deviations.md
```

Put raw bulk data and checkpoint files outside Git. Store checksums and acquisition instructions in Git. Select one environment/container strategy and record actual package versions after a successful pilot. Exact version pins cannot be meaningfully guessed before testing model compatibility.

## 2. Data schemas and integrity rules

### 2.1 `genomes.tsv`

```text
genome_id assembly_accession accession_version isolate_group source source_release retrieved_at biosample taxonomy_id taxonomy_system taxonomy_release assembly_level genome_checksum included exclusion_reason
```

One row per analyzed assembly. `isolate_group` links technical/reassembly duplicates. Missing BioSample is allowed if flagged; missing genome provenance is not silently replaced by a sequence identifier.

### 2.2 `sequences.parquet`

```text
sequence_id genome_id copy_id contig start end strand gene_length sequence_hash_rcanonical sequence species_id taxonomy_conflict extraction_tool extraction_version qc_status
```

Canonical hash is the hash of a deterministic choice between sequence and reverse complement. Preserve original orientation separately. Copy IDs are unique within a genome, not globally unique strings without a parent.

### 2.3 `panels.tsv`

```text
panel_id block_id target_species neighbor_species neighbor_selection_rule rank_relatedness application_role eligibility_status
```

A long table has one target-neighbor relationship per row. Species IDs must resolve to the frozen taxonomy. Verify that overlap of species/genomes across panels is reflected in block construction.

### 2.4 `primers.tsv`

```text
assay_id forward_name forward_sequence_5to3 reverse_name reverse_sequence_5to3 source_citation expected_region domain read_layout usable_read_length minimum_overlap mismatch_rule three_prime_rule
```

Document whether usable read length is before or after primer removal. A different oligo variant or read configuration is a different assay ID, even if its region label is unchanged.

### 2.5 `observations.parquet`

```text
observation_id sequence_id genome_id assay_id product_index template_start template_end product_length forward_mismatches reverse_mismatches terminal_mismatches compatibility_status observed_sequence observed_hash read_layout error_profile simulation_seed role
```

Keep failure rows. For paired unmerged observations, store each read and orientation explicitly rather than embedding an artificial contiguous string. If a separator is used, it must be supported/validated by the representation and matched baselines.

### 2.6 `splits.tsv`

```text
entity_id entity_type block_id global_role local_role fold_id seed split_manifest_hash
```

`global_role`: adaptation_train/development/locked_test.
`local_role`: support/query; optional support_train/support_validation for nested evaluation.

The same genome/isolate group cannot cross local support/query. The same block cannot cross global adaptation/development/test. Do not confuse region-specific identical-sequence overlap with a duplicate-genome leak; report both.

### 2.7 `predictions.parquet`

```text
run_id method_id checkpoint_hash panel_id block_id assay_id genome_id copy_id observation_id true_species true_target prediction predicted_probability abstained product_failed local_role global_role
```

Probabilities may be missing for a pure distance method. Do not serialize arbitrary distance as a probability. Include every expected query row or an explicit failure status.

### 2.8 `recommendations.parquet`

```text
run_id panel_id method_id assay_id selection_score score_definition eligible rank selected reference_count target_compatibility empirical_ambiguity uncertainty_definition status rejection_reason
```

Test utility belongs in a separate evaluation table, not in features consumed by recommendation generation.

## 3. Example prospective configuration

```yaml
study:
  version: '1.0'
  intended_use: research_assay_planning
  taxonomy_snapshot: REQUIRED
  accession_manifest_sha256: REQUIRED
  split_manifest_sha256: REQUIRED
  primary_task: target_vs_declared_neighbors
  primary_unit: taxonomic_block

panels:
  primary_neighbor_rule: genome_relatedness_plus_declared_application
  neighbor_count_target: 5
  min_distinct_genomes_per_species: 5
  underpowered_panel_policy: descriptive_only

observations:
  remove_primers: true
  retain_failure_rows: true
  primary_orientation: extracted_gene_strand
  truncated_input_policy: error
  paired_read_policy: REQUIRE_EXPLICIT_IMPLEMENTATION

features:
  kmer_k: [4, 6, 8]
  canonical_kmers: true
  primary_embedding_pooling: masked_mean_if_supported
  embedding_checkpoint_revisions: REQUIRED

splits:
  global_group: nonoverlapping_taxonomic_block
  local_group: isolate_group
  sequence_overlap_policy: report_primary_and_novelty_sensitivity
  protected_test_unseal_requires_approval: true

selection:
  primary_score: inner_support_species_balanced_binary_accuracy
  tie_tolerance: 0.01
  tie_break_order: [declared_cost, target_compatibility, assay_id]
  missing_assay_policy: explicit_failure

statistics:
  confidence_level: 0.95
  bootstrap_unit: taxonomic_block
  bootstrap_replicates: 2000
  primary_comparison: FREEZE_AFTER_DEVELOPMENT
  multiplicity: Holm_for_declared_confirmatory_family

adaptation:
  run_only_after_gate: true
  primary_loss: same_assay_supervised_contrastive
  impossible_cross_label_identity_constraints: mask
  update_test_encoder: false
  test_panel_local_support_classifier: allowed_fixed_recipe
```

All example numeric fields are proposals to confirm, not evidence-based universal defaults. A config validation step must reject `REQUIRED` placeholders before running.

## 4. Primary recommendation evaluation pseudocode

```python
for block in protected_test_blocks:
    for panel in block.panels:
        support, query = load_locked_local_roles(panel)
        # Method recipes and hyperparameters were frozen on development blocks.
        for method in frozen_method_recipes:
            scores = {}
            for assay in declared_assays:
                support_observed = observe(support, assay)
                scores[assay] = inner_grouped_support_score(
                    method, support_observed,
                    groups=support.isolate_group,
                    retain_assay_failures=True,
                )
            chosen = select_with_locked_constraints_and_ties(scores)
            # No query labels or query-derived assay oracle used above.
            fitted_local_scorer = method.fit_support_only(support, chosen)
            predictions = fitted_local_scorer.predict(observe(query, chosen))
            save_complete_predictions(predictions, including_failures=True)
            # Evaluating all assays is permitted only for post-selection metrics.
            save_evaluation_only_assay_utilities(method, support, query)
```

This is conceptual pseudocode, not a runnable library API. If a ranker is trained globally, it must be frozen before this loop. If a method has no local fitting stage, replace the stage with its fixed reference-matching rule.

## 5. Shared-sequence ambiguity audit pseudocode

```python
# weights implement the declared target/neighbor and genome weighting.
weights = observed_rows.groupby(['observed_hash', 'class_label']).weight.sum()
ceiling = weights.groupby('observed_hash').max().sum() / weights.sum()
label_sets = observed_rows.groupby('observed_hash').class_label.agg(set)
ambiguous_hashes = label_sets[label_sets.map(len) > 1]
```

Run on clean, actually observed inserts, not primer-containing templates. A single pooled binary ceiling can change with class weights. Also calculate target-versus-each-neighbor collisions. This is a reference-population audit, not a proof of universal irresolvability.

## 6. Pair-generation safeguards for adaptation

For each training assay:
1. Sample a taxonomic block.
2. Sample species with controlled weighting.
3. Sample independent genome groups where possible.
4. Create same-species positives and within-block cross-species negatives.
5. Mask exact cross-label observed identities and disputed labels.
6. Keep metadata sufficient to audit every sampled relationship.
7. Count skipped anchors and impossible constraints.
8. Do not inspect development/test queries during mining.

Assert that a species with only one genome is not presented as having multiple independent positives. Different regions are not automatic positive evidence for a same-assay task.

## 7. R analysis pattern for paired block effects

The following is an example analysis pattern, not a complete analysis of an existing dataset. Input must have one aggregate value per independent block/method after averaging nested panels and seeds according to the protocol.

```r
library(dplyr)
library(tidyr)

wide <- block_metrics |>
  select(block_id, method, utility) |>
  pivot_wider(names_from = method, values_from = utility)

stopifnot(!anyDuplicated(wide$block_id))
stopifnot(all(c('classical', 'candidate') %in% names(wide)))
stopifnot(!anyNA(wide[c('classical', 'candidate')]))

wide <- wide |> mutate(delta = candidate - classical)
set.seed(1729)
boot_delta <- replicate(2000, {
  ii <- sample(seq_len(nrow(wide)), replace = TRUE)
  mean(wide$delta[ii])
})

estimate <- mean(wide$delta)
interval <- quantile(boot_delta, c(0.025, 0.975))
```

Do not silently complete-case away unavailable assays; first implement the explicit failure policy and show missingness. If only a few independent blocks exist, interpret the bootstrap cautiously and emphasize the limited scope. A hierarchical bootstrap additionally resamples genomes within sampled blocks if the estimand includes that sampling uncertainty.

## 8. Pilot acceptance checklist

- [ ] Every sequence maps to a source accession and parent biological unit.
- [ ] Taxonomy crosswalk resolves names without outcome-driven relabeling.
- [ ] Both strands and IUPAC primers pass synthetic fixture tests.
- [ ] Primer clipping does not leak oligo sequence into classification.
- [ ] Coordinate mapping handles insertion/deletion variation.
- [ ] Failed amplification remains in the denominator where required.
- [ ] Multiple 16S copies remain linked to genomes.
- [ ] Identical cross-species inserts are visibly ambiguous.
- [ ] No input is silently truncated by a tokenizer.
- [ ] Same input/settings yield identical inference embeddings.
- [ ] Classical alignment and k-mer baselines are competitive and audited.
- [ ] Each learned model has an exact checkpoint revision and permitted license.
- [ ] Global taxonomy blocks and local genome roles are locked.
- [ ] Query labels cannot be loaded by ranking/training components.
- [ ] One end-to-end fixture runs from manifest to report.

## 9. Test-set unsealing checklist

- [ ] Primary method selection rule, comparator, endpoint, and inference unit declared.
- [ ] Hyperparameter budget and all tried model families logged.
- [ ] Split hash and accessions archived.
- [ ] Assay candidate universe fixed; ties and failures specified.
- [ ] Calibration and abstention thresholds frozen.
- [ ] Test is genuinely new to the current development cycle.
- [ ] Pretraining exposure uncertainty documented.
- [ ] Analysis scripts operate on complete predictions with failure checks.
- [ ] A negative result will be retained rather than triggering post hoc test tuning.

## 10. Compute/resource log template

```text
run_id; date; code_commit; config_hash; data_hash; checkpoint_hash;
CPU; RAM; GPU; driver; accelerator_library; numerical_precision;
observations; unique_sequences; token_distribution; batch_size;
trainable_parameters; total_parameters; optimizer_steps;
wall_time; peak_RAM; peak_GPU_memory; cache_state; failure_reason
```

Record cold-start and cached inference separately. Compute financial costs from measured runtime and the actual provider rate, not a guessed hourly estimate. Do not claim carbon savings without an explicitly documented measurement/estimation method.

## 11. Model/data card headings

Intended population; taxonomic coverage; supported assays; source releases; label provenance; known ambiguity; input lengths/strand handling; pretraining uncertainty; fine-tuning objective; held-out design; primary performance; subgroup failures; calibration; abstention behavior; runtime; licenses; unsupported diagnostic uses; maintenance policy.

## 12. Protocol deviation form

```text
Date:
Protocol/version:
Change:
Reason:
Before or after protected test inspection:
Who approved:
Affected hypotheses/endpoints:
Whether a new independent test is required:
Sensitivity analysis retaining original plan:
```

## 13. Final publication/implementation checks

- [ ] Qscore overlap discussed explicitly.
- [ ] DNABERT-S and Micro16S cited accurately, with publication status.
- [ ] Species classification, community profiling, and assay selection not conflated.
- [ ] Full-length and short-assay observations compared fairly.
- [ ] Test-derived oracle used only for evaluation.
- [ ] Paired distances/reads/seeds not used to inflate n.
- [ ] Support/reference availability distinguished from zero-shot label prediction.
- [ ] All source records have verified bibliographic exports before submission.
- [ ] Software interface examples labeled implemented only after they exist.
- [ ] Public release includes a minimal reproducible example and an explicit unresolved case.
