# S10 Text. Evidence for the pre-registration claim

The manuscript states that the analysis plan, panel definitions and eligibility criteria
were committed to version control before the confirmatory data were collected. This document
sets out the evidence, and states plainly what that evidence does and does not establish.

---

## 1. The claim, broken into checkable parts

| claim | evidence | verdict |
|---|---|---|
| A frozen analysis plan exists | commit `24ad4928b4e3`, 2026-11-08 07:20:42 −0500 | documented |
| It contained the plan, panel list and deviation log | that commit added exactly `FREEZE.md` (194 lines), `test_blocks_definition.py` (233 lines) and `DEVIATIONS.md` (75 lines), and nothing else | documented |
| It was not edited afterwards | `git log --follow FREEZE.md` returns one commit; `git diff HEAD` is empty | documented |
| Locked-test data came later | earliest raw record file mtime 2026-11-08 08:54:19; first commit containing locked-test artifacts 08:55:48 | documented, 94 min gap |
| Development panels were likewise pre-declared | commit `1918e6f` 2026-11-07 22:38:22 declared the 31 blocks; the fetch commit `69c5c19` followed at 22:53:24 | documented, 15 min gap |
| History was not rewritten | 27 commits, dates strictly non-decreasing, linear with no merges or amends | consistent |

## 2. Reproducing the check

```bash
# the freeze commit and its exact contents
git log --format='%H %ci %s' -1 24ad492
git show --stat --format="" 24ad492

# the plan has one commit in its entire history
git log --oneline --follow -- ".../phase3_freeze/FREEZE.md"

# earliest locked-test record on disk
find .../phase4_locked_test/raw_species_test -name '*.json' \
  -exec stat -f '%Sm %N' -t '%Y-%m-%d %H:%M:%S' {} \; | sort | head -1

# history is linear and chronological
git log --reverse --format='%ci' | sort -c
```

## 3. What this evidence does NOT establish

This is the part a reviewer is entitled to press, and it should not be glossed.

**Git timestamps are self-reported.** Both author and committer dates can be set to any
value through `GIT_AUTHOR_DATE` and `GIT_COMMITTER_DATE`. They are recorded by the person
making the commit, not by an independent party.

**File modification times are local metadata** and can be changed trivially.

**The repository has no remote.** Nothing was pushed to a server that would hold an
independent record of when each commit appeared. There is therefore no third-party
corroboration of any date in this history.

**In principle the entire history could have been constructed retrospectively.** Nothing in
the repository alone rules that out.

The correct characterization is therefore: *the pre-registration is thoroughly documented
and internally consistent, but it is not independently timestamped.* It rests on the
authors' integrity, in the same way that an unregistered analysis plan described in a methods
section does. It is stronger than an undocumented claim and weaker than a public
registration.

## 4. Circumstantial support

Two features of the record are consistent with genuine pre-declaration and would have been
awkward to fabricate.

**The pre-declared panels did not all work.** Of 22 locked-test blocks, 2 failed the
eligibility gate (`actinomyces`, 7 target isolates; `bartonella`, one qualifying neighbour
species). Of 31 development blocks, 5 failed. A panel list written after seeing the data
would be unlikely to contain blocks that fail the author's own criteria.

**Pre-registered predictions were wrong in recorded ways.** The a-priori difficulty classes
in the panel files are systematically mis-calibrated against the measured ceilings, and the
Extension A hypothesis, pre-registered with an explicit commitment to report
non-replication, did not replicate. Both are recorded in the committed files rather than
quietly revised.

**Deviations are logged against themselves.** `DEVIATIONS.md` records three departures from
the original protocol, each stating whether it preceded inspection of the locked test.

## 5. Recommended actions before submission

1. **Deposit the repository publicly and obtain a DOI** (Zenodo, OSF or equivalent). This
   does not retroactively timestamp the commits, but it fixes the history at a verifiable
   date and lets reviewers inspect it.
2. **Consider cryptographic timestamping** of the freeze commit hash, for example via
   OpenTimestamps, which anchors a hash to a public blockchain. Again this certifies from
   the date of anchoring onward, not retroactively.
3. **State the evidentiary basis precisely in the manuscript** rather than implying external
   verification. Suggested wording:

   > The analysis plan, panel definitions and eligibility criteria were committed to a
   > version-controlled repository before the corresponding data were retrieved; the freeze
   > commit is dated 94 minutes before the first locked-test record. These are
   > author-recorded timestamps in a repository without a remote, so they document the
   > intended sequence rather than provide independent certification. The full history is
   > available at [repository DOI].

4. **For future work, register the protocol with a third party before collecting data.** OSF
   registration or an equivalent service provides exactly the external timestamp this study
   lacks, at no cost.

## 6. Bottom line

The sequence described in the manuscript is documented in detail, internally consistent, and
supported by circumstantial features that would be difficult to arrange after the fact. It
is not independently certified, and the manuscript should say so rather than leave a reader
to assume otherwise.
