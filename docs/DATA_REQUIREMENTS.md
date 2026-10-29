# Data requirements

This document lists every input the publication pipeline reads, where it is
expected, what it is for, where it comes from, and whether it is available to a
reviewer. File and object names were taken from the scripts (static scan in
[`SCRIPT_INVENTORY.md`](SCRIPT_INVENTORY.md)) and checked by loading each file
(`Rscript tools/run_reproduction.R --mode preflight --deep`).

## Where inputs live

| Location | Default path | Override | In git? |
|---|---|---|---|
| Project data | `01_Data/` | `CHIK_DATA_DIR` | **No** (`01_Data/` is in `.gitignore`; ~56 GB on the author's machine, much of it historical) |
| Companion project | `../CHIK_vaccine_impact/CHIK_ORV_impact/` | `CHIK_ORV_ROOT` | **No** (separate, private repository) |
| Reproduced outputs | `reproduced_outputs/` | `CHIK_OUTPUT_DIR` | No (ignored; regenerated) |

The scientific scripts contain the author's absolute paths
(`c:/Users/user/OneDrive/...`). They are **not** edited; the reproduction runner
translates those paths to the locations above (see
[`tools/repro/sandbox.R`](../tools/repro/sandbox.R)). A missing input stops that
item with a message naming the file and its role; nothing is substituted.

Data types used below: **raw** (observed data), **processed** (derived tables),
**fit** (transmission-model posterior), **sim** (cached simulations),
**draws** (PSA/posterior draw-level benefit-risk outputs), **summary**
(summaries used directly by figure/table scripts), **param** (parameter
samples).

## A. Inputs for reproducing manuscript figures and tables (Mode B)

These are sufficient to regenerate every figure and table in the submission
from archived model outputs. Total size ~1.4 GB, dominated by the posterior.

| File | Location | Type | R object(s) and key columns | Purpose | Produced by | Used by (manuscript item) | Available | Redistributable? |
|---|---|---|---|---|---|---|---|---|
| `psa_grid_bra_travel_finite_setting.RData` (1 MB) | `01_Data` | summary | `psa_grid_travel_setting` (8,736 rows): `age_group, duration, entry_day, setting, AR_*, brr_{sae,death,daly}_{median,lo,hi,prob_gt1}, averted_*, excess_*, daly_*` | Setting-level travel grid (entry day x duration), posterior summaries pooled across states and draws | `06_brazil_travel_final_finite.R` | Fig. 1, Supp. Fig. 3, Supp. Data 2 | Yes (author) | Likely yes (model output); author to confirm |
| `psa_df_bra_travel_finite.RData` (12 MB) | `01_Data` | draws | `psa_df` (88,000 rows): `draw, state, days, age_group, averted_10k_*, excess_10k_*, daly_averted, daly_sae, brr_*` | Travel draw-level benefits/risks for 7/14/30/90-day trips (entry timing sampled within draw) | `06_brazil_travel_final_finite.R` | Table 1, Supp. Tables 2-3 | Yes | Likely yes; author to confirm |
| `draw_level_xy_serostatus_finite.RData` (46 MB) | `01_Data` | draws | `draw_level_xy_serostatus` (1,056,000 rows): `draw_id, Region, Coverage, outcome, AgeCat, VE_label, setting, RR_seropos, averted_10k, q_seroneg_vacc, x_10k_base, x_10k_adj, brr_base, brr_adj` | Outbreak-response (ORI) draw-level benefits and risks, base and serostatus-adjusted | `03_brazil_all_draws_ori_v3.R` | Fig. 2, Table 2, Supp. Tables 4-5, Supp. Fig. 4A-B | Yes | Likely yes; author to confirm |
| `setting_key.RData` (<1 kB) | `01_Data` | processed | `setting_key`: named character, state -> Low/Moderate/High | Transmission-setting classification of the 11 states | `08_attack_rate_classification_finite.R` | Tables 1-2 (via travel), Supp. Table 1, Supp. Fig. 2 | Yes | Yes |
| `pop_by_state.RData` (<1 kB) | `01_Data` | processed | `pop_by_state`: `region, tot_pop` | State population | author (source to be documented) | Supp. Fig. 2 | Yes | Likely yes (public census); author to confirm source |
| `state_ar_hazard_summary_finite.RData` (<1 kB) | `01_Data` | summary | `ar_summary`: `region, ar_med, ar_lo, ar_hi` | Posterior cumulative attack rate by state | author (diagnostic script; producer not identified in this repository) | Supp. Table 1 | Yes | Yes |
| `posterior_finite_all.RData` (1.3 GB) | `01_Data` **and** `CHIK_ORV_impact/00_Data/0_2_Processed/` (byte-identical copies) | fit | `posterior_*` per state (rho, base_beta, I0, ...) | Transmission-model posterior draws (finite-history immunity) | Stan fit in `CHIK_ORV_impact` | Supp. Table 1, Supp. Figs 1-2 | Yes | Model output; size may require a data repository rather than GitHub |
| `observed_2022.RData` (<1 MB) | `CHIK_ORV_impact/00_Data/0_2_Processed/` | raw/processed | `observed_all` (572 rows: `Week, Observed, Type, region`), per-state `observed_*`, `bra_all_sum` | Weekly reported symptomatic chikungunya cases, 2022, 11 states | `CHIK_ORV_impact` data processing (Brazilian national surveillance) | Supp. Figs 1-2, Supp. Table 1 | Yes (author) | **Author to confirm** (surveillance-derived) |
| `brr_weeksweep_summary_finite.RData`, `brr_weeksweep_summary_finite_ve0.RData` (<1 MB each) | `CHIK_ORV_impact/00_Data/0_2_Processed/` | summary | `brr_weeksweep_summary` (4,032 rows): `outcome, AgeCat, VE_label, setting, Coverage, week, brr_{base,adj}_{med,lo,hi,prob_gt1}` | ORI benefit-risk across campaign start week x coverage (infection+disease blocking / disease blocking only) | `CHIK_ORV_impact/01_Script/1_4_Analysis/finite_weeksweep_brr.R` | Supp. Fig. 4C-E | Yes | Likely yes; author to confirm |
| `finite_weeksweep_heatmap_probgt1.R` (code) | `CHIK_ORV_impact/01_Script/1_4_Analysis/` | code | n/a | Plotting functions sourced by `46_daly_weeksweep_probgt1_panelE.R` | `CHIK_ORV_impact` | Supp. Fig. 4C-E | Yes | Author's code; can be copied into this repository with author approval. **Note:** the copy in `02_Scripts/` differs from the one actually sourced. |

Supplementary Data 1 (parameter table,
`Manuscript/Supplementary Data table 1.xlsx`) is a compiled table with no
generating script in this repository; its values correspond to the parameter
distributions coded in `02_Scripts/01_setup.R` (LHS block) and
`02_setup_age_props.R`, but this correspondence has **not** been checked
cell by cell.

## B. Additional inputs for the upstream pipeline in this repository (Mode C)

Needed only to regenerate the draw-level outputs in section A.

| File | Type | Purpose | Read by | Available | Notes |
|---|---|---|---|---|---|
| `combined_nnv_national_age_ixchiq.RData`, `mortality_chikv.RData`, `rho_df.RData`, `all_risk_four_age.csv` | processed | Age-specific risk inputs, reporting rates | `01_setup.R` (sourced by 03 and 06) | Yes | |
| `lhs_sample.RData` | param | Cached LHS parameter sample (1,000 draws) | `03_brazil_all_draws_ori_v3.R` | Yes | `01_setup.R` also regenerates the sample with `set.seed(123)`; `regenerate_lhs_sample_cache.R` rebuilds the cache |
| `chikv_fatal_hosp_rate.RData` | processed | Infection fatality / hospitalisation rates | 03 | Yes | |
| `le_sample.RDS`, `lhs_sample_young.RDS`, `lhs_old.RDS` | param | Life expectancy and age-specific parameter samples | 03 | Yes, **but stored as lower-case `.rds`** | Scripts request `.RDS`. Works on Windows and default macOS (case-insensitive), **fails on Linux**. Fix by renaming the files or the references (author decision). |
| `lhs_orv.RData` | param | ORI parameter sample | 03 | Yes | |
| `preui_flat_all.RData` (335 MB) | sim | Pre-vaccination simulations | 03 | Yes | From `CHIK_ORV_impact` |
| `sim_results_vc_ixchiq_model_finite.RData` (1.4 GB), `postsim_vc_ixchiq_model_finite.RData` (2.7 GB), `combined_nnv_df_region_coverage_model_finite.RData` | sim | Outbreak-response vaccination simulations | 03 (and 29-31 diagnostics) | Yes | From `CHIK_ORV_impact` |
| `observed_2022.RData` | raw/processed | as above | `08_attack_rate_classification_finite.R` | Yes (ORV) | |

## C. Upstream of this repository (companion project `CHIK_ORV_impact`)

The age-structured transmission model fit (Stan), the outbreak-response
simulations and the campaign week-sweep are produced in `CHIK_ORV_impact`
(`01_Script/1_1_Functions/sim_functions_final.R`,
`01_Script/1_4_Analysis/finite_weeksweep_{fulldraws,nnv,brr}.R`, and the
fitting code). That project is not part of this repository, its raw inputs
(municipality-week case data, population estimates) are not audited here, and
re-running it is outside the scope of the reproduction tooling.
**Status: BLOCKED for independent reviewers** until the author decides which
parts of `CHIK_ORV_impact` to release.

Copies of four week-sweep scripts exist in `02_Scripts/`
(`finite_weeksweep*.R`). Only `finite_weeksweep_nnv.R` is identical to the
companion-project version; `finite_weeksweep_brr.R`,
`finite_weeksweep_fulldraws.R`, `finite_weeksweep_heatmap.R` and
`finite_weeksweep_heatmap_probgt1.R` differ. The companion-project versions are
the ones the publication outputs depend on.

## How a reviewer obtains the data

Not yet decided by the authors. Options, in order of preference:

1. Deposit the section-A files (~1.4 GB) in a data repository (e.g. Zenodo)
   and point `CHIK_DATA_DIR` / `CHIK_ORV_ROOT` at the downloaded folder; this
   is enough for Mode B.
2. Provide them to reviewers through the journal's confidential channel.
3. For full regeneration, additionally release the section-B files and the
   relevant parts of `CHIK_ORV_impact`.

Any surveillance-derived file (`observed_2022.RData` and its raw sources)
needs the authors' confirmation that redistribution is permitted.
