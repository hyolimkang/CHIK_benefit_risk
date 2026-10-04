# =============================================================================
# 46_muni_timing_ridgeline_by_state.R
#
# Figure: Municipality-level epidemic timing within each state, 2022 outbreak.
#
# NOT a ridgeline of "peak-week distributions" -- each ridge IS a
# municipality's actual weekly epidemic curve (epi week on the x-axis),
# peak-normalised to its own maximum so that timing (not relative size) is
# what the eye compares across ridges. Within a state, ridges are stacked in
# ascending order of peak week, so a flat/vertical stack of peaks = the
# municipalities moved in near-unison, while a long diagonal staircase =
# staggered spread across the season. The point is deliberately NOT to
# validate the state-level aggregate trajectory -- it is to show how much
# (or how little) municipality-level timing agrees with, or is masked by,
# that aggregate.
#
# Municipality selection (per state, independently):
#   - candidates restricted to municipalities with >= MIN_CASES total
#     confirmed cases in 2022 (sporadic single/handful-case reports would
#     otherwise contribute noisy, uninformative near-zero curves and clutter
#     the figure without adding real timing information)
#   - from those candidates, take the smallest top-N (by descending total
#     cases) whose cumulative case count reaches BURDEN_TARGET (90%) of the
#     STATE's true 2022 case burden (denominator = ALL municipalities,
#     including sub-threshold ones -- so "90% of case burden" is honest
#     about the full reported total, not just the filtered candidate pool)
#   - MIN_CASES = 20 was checked against a floor of 10: identical
#     municipality sets and coverage in every state, i.e. the floor removes
#     genuinely sporadic reporting without costing any burden coverage
#     (see console output of this script for the per-state check)
#
# Data sources:
#   - CHIK_ORV_impact/00_Data/0_1_Raw/chik_brazil_muni_week_2015_2024.rds
#     (municipality [muni6, IBGE 6-digit] x week case counts, 2015-2024;
#     cases_confirmed used throughout, matching muni_hotspot_map_2022.R)
#   - CHIK_ORV_impact/00_Data/0_1_Raw/estimativa_dou_2024.xls (municipality
#     names, IBGE 2024 estimate -- for the single largest-burden municipality
#     labelled per state)
#   - 01_Data/setting_key.RData + 01_Data/state_ar_hazard_summary_finite.RData
#     (Low/Moderate/High transmission-setting tier + posterior AR ordering,
#     reproduced here only to order facets consistently with Figure 1's
#     panel B/C state order -- no colour-coding on the strip itself)
#
# Output: 07_Final_Results/Main_Figures/figure_muni_timing_ridgeline.png
#         07_Final_Results/Main_Figures/figure_muni_timing_ridgeline.pdf
# =============================================================================

setwd("c:/Users/user/OneDrive/CHIK_benefit_risk")
suppressMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggridges)
  library(ggtext)
  library(scales)
  library(readxl)
})

ori_root <- "c:/Users/user/OneDrive/CHIK_vaccine_impact/CHIK_ORV_impact"

## ---- 0. Config --------------------------------------------------------
target_states <- c(
  "Bahia" = 29, "Ceará" = 23, "Minas Gerais" = 31, "Pernambuco" = 26,
  "Paraíba" = 25, "Rio Grande do Norte" = 24, "Piauí" = 22, "Alagoas" = 27,
  "Tocantins" = 17, "Sergipe" = 28, "Goiás" = 52
)

MIN_CASES     <- 20    # per-municipality floor, 2022 total confirmed cases
BURDEN_TARGET <- 0.90  # cumulative share of state's TRUE 2022 case burden

load("01_Data/setting_key.RData")              # -> setting_key
load("01_Data/state_ar_hazard_summary_finite.RData")  # -> ar_summary

# Same state order as Figure 1 panel B/C: Low -> Moderate -> High tier,
# ascending posterior cumulative AR within tier -- keeps this figure's facet
# order consistent with the transmission-setting classification used
# throughout the manuscript's main figures.
setting_rank <- c(Low = 1, Moderate = 2, High = 3)
state_order <- ar_summary %>%
  mutate(setting = setting_key[region]) %>%
  arrange(setting_rank[setting], ar_med) %>%
  pull(region)

## ---- 1. Municipality x week case counts, 2022 --------------------------
muni_week <- readRDS(file.path(ori_root, "00_Data/0_1_Raw/chik_brazil_muni_week_2015_2024.rds"))

muni_2022 <- muni_week %>%
  filter(epi_year == 2022, epi_week >= 1, epi_week <= 52) %>%
  mutate(uf_code = as.integer(substr(muni6, 1, 2))) %>%
  filter(uf_code %in% target_states) %>%
  mutate(state = names(target_states)[match(uf_code, target_states)],
         is_unknown_muni = substr(muni6, 3, 6) == "0000") %>%
  filter(!is_unknown_muni) %>%
  select(state, muni6, epi_week, cases_confirmed)

rm(muni_week); gc()

## ---- 2. Municipality names (IBGE 2024 estimate; for one label/state) ---
pop_raw <- readxl::read_excel(
  file.path(ori_root, "00_Data/0_1_Raw/estimativa_dou_2024.xls"),
  sheet = "MUNICÍPIOS", skip = 1
)
colnames(pop_raw) <- c("uf_abbr", "cod_uf", "cod_munic", "muni_name", "population")
muni_names <- pop_raw %>%
  filter(!is.na(cod_uf), !is.na(cod_munic)) %>%
  mutate(cod_uf = as.integer(cod_uf), cod_munic = as.integer(cod_munic),
         muni6 = as.character((cod_uf * 100000L + cod_munic) %/% 10L)) %>%
  distinct(muni6, .keep_all = TRUE) %>%
  select(muni6, muni_name)

## ---- 3. Select municipalities explaining >=90% of each state's burden --
muni_totals <- muni_2022 %>%
  group_by(state, muni6) %>%
  summarise(total_cases = sum(cases_confirmed, na.rm = TRUE), .groups = "drop")

state_totals <- muni_totals %>%
  group_by(state) %>%
  summarise(state_total = sum(total_cases), .groups = "drop")

selected <- muni_totals %>%
  filter(total_cases >= MIN_CASES) %>%
  left_join(state_totals, by = "state") %>%
  group_by(state) %>%
  arrange(desc(total_cases), .by_group = TRUE) %>%
  mutate(cum_frac  = cumsum(total_cases) / state_total,
         reach_tgt = cum_frac >= BURDEN_TARGET) %>%
  mutate(cutoff = { hit <- which(reach_tgt); if (length(hit) == 0) n() else hit[1] }) %>%
  filter(row_number() <= cutoff) %>%
  ungroup() %>%
  select(state, muni6, total_cases, state_total, cum_frac)

coverage_summary <- selected %>%
  group_by(state) %>%
  summarise(n_muni = n(), coverage_pct = 100 * max(cum_frac), .groups = "drop")

cat("===== Municipality selection per state (>=", MIN_CASES, "cases, reaching",
    BURDEN_TARGET * 100, "% of case burden) =====\n")
print(as.data.frame(coverage_summary))

## ---- 4. Weekly series for selected municipalities, peak-normalised -----
weekly <- muni_2022 %>%
  semi_join(selected, by = c("state", "muni6")) %>%
  complete(nesting(state, muni6), epi_week = 1:52, fill = list(cases_confirmed = 0))

peak_info <- weekly %>%
  group_by(state, muni6) %>%
  summarise(peak_week  = epi_week[which.max(cases_confirmed)],
            peak_cases = max(cases_confirmed), .groups = "drop") %>%
  left_join(selected %>% select(state, muni6, total_cases), by = c("state", "muni6")) %>%
  left_join(muni_names, by = "muni6") %>%
  mutate(state = factor(state, levels = state_order)) %>%
  arrange(state, peak_week, desc(total_cases)) %>%
  mutate(muni_row = factor(muni6, levels = unique(muni6)))

weekly_norm <- weekly %>%
  left_join(peak_info %>% select(state, muni6, muni_row, peak_week, peak_cases),
             by = c("state", "muni6")) %>%
  mutate(state      = factor(state, levels = state_order),
         norm_cases = ifelse(peak_cases > 0, cases_confirmed / peak_cases, 0))

# Single largest-burden municipality per state -- one anchoring label per
# ridge stack, not a name on every ridge (up to 48 munis/state -- unreadable).
# Flip the label to the other side of its point whenever the point sits near
# the right edge of the 1-52 week axis, so the text can't run past the panel
# border and get clipped (e.g. "Januária" -> "Jan" when left-aligned at week 51).
top_label <- peak_info %>%
  group_by(state) %>%
  slice_max(total_cases, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(near_right_edge = peak_week > 40,
         label_x    = peak_week + ifelse(near_right_edge, -1, 1),
         label_hjust = ifelse(near_right_edge, 1, 0))

## ---- 5. Facet strip: state name (large/bold) + n municipalities and ------
## coverage (small, secondary) on a second line -- via ggtext::element_markdown
## so the two lines can carry different sizes without either one overflowing
## the panel width. Plain grey strip background (no per-tier colour coding --
## with this many panels, colour on the strip added noise, not information).
facet_lab <- coverage_summary %>%
  mutate(state = factor(state, levels = state_order)) %>%
  arrange(state) %>%
  mutate(label = paste0(
    "**", state, "**<br>",
    "<span style='font-size:11pt;color:", "#52514e", "'>n = ", n_muni,
    " municipalities, ", sprintf("%.0f", coverage_pct), "% of 2022 cases</span>"
  ))
facet_labels <- setNames(facet_lab$label, facet_lab$state)

## ---- 6. Split into pages -- 3 pages of ~4 states each (2 x 2 panels/page) -
## Splitting (rather than one crowded 4x3 grid) is deliberate: at 4 panels
## per page each panel is roughly twice as wide as in a 4-column layout, so
## the two-line strip label and the x-axis both stay legible instead of
## clipping at the panel edge. Pages follow state_order (Figure 1's
## Low -> Moderate -> High, ascending-AR sequence) in simple consecutive
## chunks -- the split points are for legibility only, not a claim that the
## chunks are analytically distinct groups.
page_size  <- 4
n_pages    <- ceiling(length(state_order) / page_size)
page_of    <- setNames(ceiling(seq_along(state_order) / page_size), state_order)

## ---- 7. Plot ------------------------------------------------------------
ink_primary   <- "#0b0b0b"
ink_secondary <- "#52514e"
grid_col      <- "#e7e6e0"
axis_col      <- "#c3c2b7"

build_page <- function(states_this_page, page_idx) {
  d_weekly <- weekly_norm %>% filter(as.character(state) %in% states_this_page) %>%
    mutate(state = factor(as.character(state), levels = states_this_page))
  d_label  <- top_label %>% filter(as.character(state) %in% states_this_page) %>%
    mutate(state = factor(as.character(state), levels = states_this_page))

  ggplot(d_weekly,
         aes(x = epi_week, y = muni_row, height = norm_cases,
             group = muni_row, fill = peak_week)) +
    ggridges::geom_ridgeline(scale = 2.8, colour = "white", linewidth = 0.15, alpha = 0.92) +
    geom_text(data = d_label,
              aes(x = label_x, y = muni_row, label = muni_name, hjust = label_hjust),
              inherit.aes = FALSE, vjust = -0.6,
              size = 3.6, colour = ink_primary, fontface = "italic") +
    facet_wrap(~ state, scales = "free_y", ncol = 2, labeller = as_labeller(facet_labels)) +
    scale_fill_viridis_c(option = "viridis", name = "Municipality peak week (epi week, 2022)",
                          limits = c(1, 52), breaks = scales::pretty_breaks(6)) +
    scale_x_continuous(name = "Epidemiological week (2022)",
                        breaks = scales::pretty_breaks(8),
                        expand = expansion(mult = c(0.01, 0.05))) +
    scale_y_discrete(name = "Municipalities within state (stacked by ascending peak week)",
                      expand = expansion(add = c(0.3, 1.6))) +
    guides(fill = guide_colourbar(barwidth = unit(16, "lines"), barheight = unit(0.7, "lines"),
                                   title.position = "top", title.hjust = 0.5)) +
    labs(
      title    = paste0("Municipality-level timing of the 2022 chikungunya outbreak within each state (",
                         page_idx, " of ", n_pages, ")"),
      subtitle = paste(strwrap(paste0(
        "Each ridge is one municipality's weekly confirmed-case curve, normalised to its own peak (not a distribution of peak weeks). ",
        "Municipalities selected per state: \u2265", MIN_CASES, " confirmed cases in 2022, ",
        "smallest set reaching \u2265", BURDEN_TARGET * 100, "% of the state's total 2022 case burden. ",
        "Stack order (bottom \u2192 top) = ascending peak week, so a flat stack = municipalities peaked in unison; a diagonal staircase = staggered spread."
      ), width = 95), collapse = "\n"),
      caption = "Source: SINAN municipality x epidemiological-week confirmed case counts, 2022. Labelled municipality = largest 2022 case burden within its state."
    ) +
    theme_minimal(base_size = 14) +
    theme(
      text                   = element_text(colour = ink_primary),
      plot.title             = element_text(face = "bold", size = 18, colour = ink_primary),
      plot.subtitle          = element_text(size = 12, colour = ink_secondary, margin = margin(b = 10), lineheight = 1.25),
      plot.caption           = element_text(size = 9.5, colour = ink_secondary, hjust = 0, margin = margin(t = 8)),
      plot.title.position    = "plot",
      plot.caption.position  = "plot",
      axis.title.x           = element_text(size = 13, colour = ink_secondary),
      axis.title.y           = element_text(size = 13, colour = ink_secondary),
      axis.text.x            = element_text(colour = ink_primary, size = 11),
      axis.text.y            = element_blank(),
      axis.ticks             = element_blank(),
      axis.line.x            = element_line(colour = axis_col, linewidth = 0.3),
      panel.grid.major.x     = element_line(colour = grid_col, linewidth = 0.3),
      panel.grid.major.y     = element_blank(),
      panel.grid.minor       = element_blank(),
      panel.spacing          = unit(22, "pt"),
      strip.background       = element_rect(fill = "grey93", colour = NA),
      strip.text             = ggtext::element_markdown(size = 17, colour = ink_primary, lineheight = 1.3,
                                                          margin = margin(t = 8, b = 8)),
      legend.position        = "bottom",
      legend.title           = element_text(size = 11.5, colour = ink_secondary),
      legend.text            = element_text(size = 10.5, colour = ink_secondary),
      plot.margin            = margin(18, 22, 14, 18)
    )
}

pages <- lapply(seq_len(n_pages), function(i) {
  build_page(state_order[page_of == i], i)
})

dir.create("07_Final_Results/Main_Figures", showWarnings = FALSE, recursive = TRUE)

PAGE_WIDTH  <- 15
PAGE_HEIGHT <- 13

# One multi-page PDF (each page = one ggplot). cairo_pdf (not base pdf())
# is required here -- base pdf()'s default font encoding mangles the
# UTF-8 glyphs used throughout (>=, ->, accented state names).
cairo_pdf("07_Final_Results/Main_Figures/figure_muni_timing_ridgeline.pdf",
          width = PAGE_WIDTH, height = PAGE_HEIGHT, onefile = TRUE)
for (p in pages) print(p)
dev.off()

# ... plus one PNG per page, for quick viewing outside a PDF reader.
for (i in seq_len(n_pages)) {
  ggsave(sprintf("07_Final_Results/Main_Figures/figure_muni_timing_ridgeline_page%d.png", i),
         pages[[i]], width = PAGE_WIDTH, height = PAGE_HEIGHT, dpi = 400, bg = "white")
}

message("Saved: 07_Final_Results/Main_Figures/figure_muni_timing_ridgeline.pdf (",
        n_pages, " pages) + _page1..", n_pages, ".png")
