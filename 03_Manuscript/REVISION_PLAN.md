# CHIK Benefit-Risk Manuscript Revision Strategy

## Context

The manuscript ("Benefit-risk of Ixchiq for travel vaccination and age-specific outbreak
response immunisation") was reviewed by 5 reviewers for PLOS Medicine. Reviewer #2
recommended **reject**, centred on one issue: the baseline immunity / susceptible pool
(S0) is fixed *before* fitting the 52-week, per-state transmission rate, and those free
weekly parameters are flexible enough to absorb any error in S0 — so a visually good fit
is not evidence that S0, and therefore the herd-effect benefit attributed to
infection-blocking vaccination in older adults, is correct. Reviewer #5 raised a related,
independent concern about the catalytic model's smooth age-immunity assumption vs.
chikungunya's real step-like, outbreak-driven immunity. Reviewers #2 and #4 also flagged
the manuscript's density as a major barrier: 14 dense probabilistic-cloud panels across
Figure 1/Figure 2, plus two ~70-row tables, "function as model diagnostics rather than
figures suited to a main text article."

The user has since made the critical structural fix in the companion simulation repo
(`CHIK_vaccine_impact/CHIK_ORV_impact`): the forward vaccination-scenario simulation was
silently drawing baseline immunity from an unrelated long-term-average FOI distribution,
while the transmission-rate *fitting* stage for the "finite-history" model was already
using short-term, state-specific FOI (`bra_state_short_term_foi$H_median`). This was a
bug, not just an unvalidated assumption — the two stages weren't even internally
consistent. The fix makes both stages draw from the same short-term FOI source
(`make_region_foi_draws_shortterm()`, `foi_draws_list_shortterm`), and the corrected
pipeline now feeds `postsim_vc_ixchiq_model_finite.RData` into this repo's
`03_brazil_all_draws_ori_v3.R` (`DATA_MODE <- "flat"`), which has been run through to a
regenerated `ft_brr` table.

Before those new numbers go into a revision, they need a fast sanity check (not a full
formal parameter-recovery exercise — that's out of scope for this pass). Separately, the
figure/table structure needs a genuine redesign, not just a re-shuffle, aimed at a
compact, high-quality "one clear look" storyline in the spirit of Reviewer #5's own
suggestion: *"a summary figure that integrates the key results and clearly highlights the
combinations of transmission setting, vaccination strategy, and age group for which
vaccination provides a net benefit."*

This plan covers: (A) what to check to confirm the new finite-history results are sound,
(B) a proposed compact figure/table storyline, and (C) the concrete code changes each
piece needs, mapped to files already found in this repo. **Table 1/Table 2 unification
(travel vs. outbreak) is deliberately left as an open decision** — the user wants to
think through that separately before committing.

---

## Part A — Fast consistency check (before trusting the new `ft_brr`)

Goal: catch anything that looks *wrong*, not prove the model is right. Everything below
reuses existing objects/toggles — no new simulation code needed, only a short diagnostic
script and a couple of comparison runs.

1. **FOI self-consistency (the actual bug that was fixed).**
   Confirm the forward-simulation FOI draws now agree with the fitting-stage FOI by
   construction: pull `foi_draws_list_shortterm[[region]]` (built in
   `age_struc_fitting_region_func_updated.R`) and compare its median/95% range per state
   against `bra_state_short_term_foi$H_median` / `H_lo` / `H_hi` directly. They should
   match closely (it's a `truncnorm` fit around those exact values) — a mismatch means
   the wiring is still off somewhere.

2. **Old vs. new S0, by state.**
   In this repo, `DATA_MODE` already toggles between `"original"` (long-term average
   FOI) and `"flat"` (finite/short-term FOI). Run the age-seroprevalence calculation
   (`R0_vec` / `sero_vec` equivalent) for both modes per region and plot the two
   age-immunity curves overlaid, one panel per state (11 small multiples). Sanity check:
   states with a large, recent (2022) outbreak should show *higher* short-term immunity
   than the long-term-average assumption implied, and vice versa for states with little
   recent transmission — i.e. the direction of the shift should track `H_median`
   magnitude, not be random noise.

3. **Headline number check: does the reviewer-challenged finding survive?**
   The single most important number to check is the one Reviewer #2 explicitly
   questioned: *older adults (65+) only show BRR > 1 under the infection-blocking
   mechanism.* Filter the new `brr_table_final_long` / `ft_brr` to `Age group == "65+"`
   and compare `BRR_base`/`BRR_adj` and `prob_base`/`prob_adj` for both mechanisms,
   across settings, against the same slice from a `DATA_MODE = "original"` run. Flag
   explicitly if the qualitative conclusion (disease-blocking-only insufficient for 65+)
   changes sign anywhere — that's the number the revision response letter will need to
   speak to directly.

4. **Burden-level plausibility, not just BRR.**
   Compare `all_draws_ix_true` / `all_draws_hosp_true` / `all_draws_fatal_true` /
   `all_draws_daly_true` pre-vaccination totals (flat vs. original mode) per region.
   Expect a level shift (new S0 changes remaining susceptible pool → changes projected
   future burden) but not an implausible jump (e.g. >10x) and not a sign flip in
   cross-state ranking (a state that had the highest projected burden under the old
   assumption shouldn't suddenly rank lowest without a clear reason tied to its
   `H_median`).

5. **Data hygiene re-check.**
   Quick grep confirmation that no other object in the pipeline still uses the stale
   18-age-bin scheme (the `hosp`/`fatal`/`nh_fatal` bug already found and fixed) — check
   `length()` of `hosp`, `fatal`, `nh_fatal`, and any `age_groups = 1:18` default
   argument left in `age_struc_fitting_region_func_updated.R` or
   `03_brazil_all_draws_ori_v3.R`. Also check `draw_level_xy_serostatus` and
   `brr_table_final_long` for `NA`/`Inf`/negative values in `brr_base`, `brr_adj`,
   `Benefit`, `Risk_base`, `Risk_adj` before those numbers go in the manuscript.

6. **Travel scenario — confirm it's out of scope for this fix.**
   The travel-vaccination engine (`compute_ar()`/`compute_outcome()` in `01_setup.R`,
   driven by `06_brazil_travel_final.R`/`07_brazil_travel_final_update.R`) uses the
   fitted time-varying daily FOI from the actual 2022 outbreak directly, not the
   catalytic-model S0. Confirm this is still true (it should be unaffected by the S0
   fix) so the plan doesn't waste time re-validating numbers that didn't change — but
   also confirm Reviewer #1's flagged inconsistency (line 315: text says "risks
   outweigh benefits" for 90-day/18–64/low-transmission while BRR > 1, 91.5% probability
   of benefit) is a **wording bug, not a numbers bug**, so it can be fixed in text
   without touching the travel pipeline.

**Output of Part A:** a short diagnostic script (`02_Scripts/09_verify_finite_history.R`)
producing (a) the 11-state old-vs-new S0 overlay/comparison, (b) a delta table for the
65+ mechanism-sensitivity headline numbers, and (c) a data-hygiene check printed to
console. This script doubles as raw material for a **supplementary robustness figure**
answering Reviewers #2/#5 directly in the revision.

---

## Part B — Storyline & figure/table redesign (Nature Medicine–compact)

### The one-sentence message to design around

Vaccination benefit vs. risk is *conditional* — on setting (transmission intensity), age,
and (for outbreak response) protection mechanism — and the paper's job is to let a reader
see that whole conditional structure in one or two looks, with the quantitative detail
available but not forced on them. That's a decision-matrix problem, not a distribution
problem, which is why the current cloud-scatter-per-outcome approach (14 panels) fights
the story rather than telling it.

### Proposed structure (3 main figures + 1 compact table; full PSA output to supplement)

**Figure 1 — Benefit–risk decision map** *(new figure, the headline)*
- Two panels: **1a Travel vaccination**, **1b Outbreak response immunisation**.
- Grid layout: rows = age group (18–64, 65+), columns = setting (Low/Moderate/High) ×
  the secondary driver for each scenario (travel duration for 1a; protection mechanism
  for 1b — could be a within-cell split or paired columns).
- Fill = probability(BRR > 1) for **DALY** (the outcome the paper's own Discussion
  already treats as most policy-relevant, since it captures chronic morbidity); tile
  label = median BRR. This is a direct, higher-quality realization of what
  `p_brr_heatmap` (`05_brr_outbreak_heatmap.R`) already does for one outcome/one
  scenario — extended to cover both scenarios and answering Reviewer #5's minor comment
  almost verbatim.
- This single figure replaces the current Fig 1a + Fig 2a and gives the "so what" in one
  glance, before any distributional detail.

**Figure 2 — Benefit–risk ratios with uncertainty** *(new figure, point-interval style)*
- Forest/point-interval plot: one row per Outcome × Setting × Age × (Mechanism or
  Duration), median BRR (log x-axis) with 95% UI, vertical reference line at BRR = 1.
  Small-multiple columns for DALY / Death / SAE (3 side-by-side panels, or 3 stacked).
- This is exactly what Reviewer #2 asked for by name ("point and interval type standard
  in benefit risk reporting") and is the direct replacement for the 6 dense cloud-scatter
  panels (1b–d, 2b–d) that reviewers singled out as unreadable and as 65MB files.
- Needs an explicit, stated convention for the zero-risk cells Reviewer #4 flagged
  (`brr_base`/`brr_adj` = `NA` when the caused-outcome denominator is 0) — e.g. show an
  open-ended arrow / right-censored point rather than silently relabeling the cell
  "beneficial." Worth deciding this convention once, in the shared plotting helper.

**Figure 3 — Acceptability curves, compacted** *(existing figure, restructured)*
- Keep CEAC-style Pr(BRR > threshold) curves but as **one** multi-panel figure (3 columns
  for DALY/Death/SAE) instead of 6 separate files, combining travel and outbreak by
  colour/linetype rather than duplicating the whole figure. Matches Reviewer #4's
  concrete suggestion ("CEAC plots could make a 3-panel Figure 2"). Candidate for
  supplement if Figure 2's point-intervals are judged sufficient on their own — flag as
  a call to make once Figures 1–2 are drafted and reviewed for redundancy.

**Table 1 — Compact BRR summary (DALY-focused, main text)**
- Reviewer #2's explicit ask: one focused table restricted to DALY, full breakdown
  (all outcomes × all settings) moved to supplement in its entirety.
- Mechanically: filter `brr_table_final_long` (outbreak) / `ar_long_table` (travel) to
  `Outcome == "DALY"` before building the `flextable` — cuts each from ~70 rows to
  roughly 12–24 depending on how travel/outbreak are laid out.
- **Left open, per the author's request**: whether Table 1 stays as two parallel tables
  (current `ft_brr`/`ft_ar` structure, just DALY-filtered) or gets unified into one
  table with a `Scenario` column. Revisit after Figure 1/2 drafts make it clearer how
  much table detail readers still need alongside the decision-map figure.

### Supplement mapping (nothing gets deleted, just relocated)
Full outcome × setting × age × mechanism/duration tables (current Table 1/2 in full);
the 6 cloud-scatter figures (`brr_ori_*_cloud.pdf`, `brr_travel_*_cloud.pdf`), ideally
hexbinned/density-plotted per Reviewer #4 rather than 65MB point clouds if kept at all;
coverage-sensitivity tables (`BRR_table_ori_setting_cov10/cov90.docx`); the Part A
robustness diagnostics (old-vs-new S0 overlay, headline-number delta table); CEAC raw
data (`ceac_travel.xlsx` etc.); full LHS parameter tables.

### Other reviewer items worth tracking (not figure/table, quick checklist for later)
Text-only or scope-limited fixes to not lose track of while focused on figures: define
"Moderate"/acronyms self-contained; disambiguate "state" (Reviewer #1 #7); add
model-type sentence (deterministic ODE/SEIR) and a model-calibration paragraph to main
text (Reviewer #1 #1–2); add key-parameter table + DW/duration references (Reviewer #1
#3–4); fix the travel/low-transmission wording inconsistency (Reviewer #1 #13); add a
funding-source statement (Reviewer #1 #10); add VAERS-limitation caveats and the 3
missing references Reviewer #3 supplied; add a model-validation paragraph (Reviewer #2
recommends showing fit isn't just absorbing S0 error — ties back to Part A); a README
for the pipeline and removal of duplicate/dead script versions (Reviewer #4 #2, separate
from this figure-focused pass — the `old codes/` folder and multiple `brazil_all_draws*`
versions).

---

## Part C — Concrete code changes

| Deliverable | New/reused code |
|---|---|
| Part A diagnostics | New `02_Scripts/09_verify_finite_history.R`: reruns S0 calc for both `DATA_MODE`s, builds the 11-panel overlay plot, builds the 65+ mechanism-sensitivity delta table, runs the NA/Inf/negative hygiene check. |
| Figure 1 (decision map) | Extend `prep_brr_heatmap_data()`/`plot_brr_outbreak_heatmap()` (`05_brr_outbreak_heatmap.R`) to loop over all 3 outcomes and both scenarios instead of one hardcoded `OUTCOME`; build an equivalent `prep_brr_heatmap_data_travel()` from `ar_summary_all`/`pr_gt1_wide` (`07_brazil_travel_final_update.R`) — no travel heatmap currently exists, this is new. |
| Figure 2 (forest/point-interval) | New `plot_brr_forest()`, built on the already-available `brr_draw_summary_true_filtered` (outbreak, `03_brazil_all_draws_ori_v3.R`); needs an equivalent travel-side summary table built the same way from `ar_summary_all`. Pure plotting code — the median/95%CI stats already exist. |
| Figure 3 (compact CEAC) | Reuse `plot_brr_ceac_outbreak_ve()` / `plot_brr_ceac()`; combine outputs via `patchwork` (already loaded in `01_setup.R`) instead of separate `ggsave` calls per outcome. |
| Table 1 (compact) | Filter `brr_table_final_long` / `ar_long_table` to `Outcome == "DALY"` before the existing `flextable::flextable()` calls (`03_brazil_all_draws_ori_v3.R:1967`, `07_brazil_travel_final_update.R:436`). |
| **Blocker to flag** | `06_brazil_travel_final.R` currently calls undefined functions (`fn_br_space_*`, `fn_panel_range`, `fn_br_grid`, `fn_br_summ`, a differently-signatured `plot_brr_outcome`) — it will not run standalone. These need restoring (a matching definition exists in `02_Scripts/old codes/brazil_travel.R:808` as a starting point) before the new travel heatmap/forest figures can be built from its outputs. |

---

## Verification

1. `02_Scripts/09_verify_finite_history.R` runs end-to-end and produces the S0
   comparison table + delta table.
2. The Figure 1 decision-map prototype renders for at least the outbreak-response side
   (existing `p_brr_heatmap` code path) before extending to travel.
3. The compact Table 1 (DALY-filtered) row count drops from ~70 to the expected
   12–24 range, confirming the filter is wired correctly against `brr_table_final_long`.
4. Once Figures 1–2 exist, do a visual gut-check against the cover letter's two headline
   conclusions (travel: 18–64y + ≥30 days + moderate/high only; outbreak: 18–64y all
   settings, 65+ only moderate/high + infection-blocking) — the decision map in Figure 1
   should make both conclusions readable without reading the text.

---

## Progress log

- **2026-08-11**: Part A checks 1, 2, 4 (partial), 5 run. See PR/commit history and
  session notes for the finding that median short-term FOI (`H_median`) is running
  ~5x higher than the long-term-average FOI across all 11 states, which is a much
  larger baseline-immunity shift than a routine sensitivity change — flagged for
  author review before proceeding to Figure/Table build-out.
