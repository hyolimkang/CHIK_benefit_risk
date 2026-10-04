suppressMessages({
  library(dplyr)
})

rdata <- "01_Data/draw_level_xy_serostatus_finite_weeksweep.RData"
if (!file.exists(rdata)) {
  stop("RData not found: ", rdata)
}
loaded <- load(rdata)
cat("Loaded objects:", paste(loaded, collapse=", "), "\n")

obj <- NULL
if (exists("draw_level_xy_serostatus")) obj <- draw_level_xy_serostatus
if (is.null(obj) && exists("all_weeks_brr")) obj <- all_weeks_brr
if (is.null(obj)) stop("No expected object (draw_level_xy_serostatus or all_weeks_brr) found in RData")

cat("Rows,cols:", nrow(obj), ncol(obj), "\n")
cat("Column names:\n")
print(names(obj))

# detect infection / symptomatic / cumulative columns
cols <- names(obj)
prefer <- c("symp_nv_10k","symp_10k","inf_10k","inf_10000","AR_travel","symp_per10k","x_10k_base","cum_inf","cum_attack","attack_cum")
chosen <- intersect(prefer, cols)[1]
if (is.na(chosen) || is.null(chosen)) {
  chosen <- cols[grepl("symp|inf|attack|ar|AR|cum", cols, ignore.case=TRUE)][1]
}
if (is.na(chosen) || is.null(chosen)) stop("No infection-like column found")
cat("Chosen column for remaining-risk:", chosen, "\n")

# ensure week and Region exist
if (!('week' %in% cols)) stop('no week column')
if (!('Region' %in% cols)) stop('no Region column')

# compute median final (week 52) per region
finals <- obj %>% filter(week == 52) %>% group_by(Region) %>% summarise(final_med = median(.data[[chosen]], na.rm=TRUE)) %>% arrange(desc(final_med))
print(head(finals, 6))
if (nrow(finals) == 0) stop('no final rows')

top_region <- finals$Region[1]
cat('Top region selected:', top_region, '\n')

series <- obj %>% filter(Region == top_region) %>% group_by(week) %>% summarise(med = median(.data[[chosen]], na.rm=TRUE)) %>% arrange(week)

# convert per10k to percent if magnitude suggests per10k (>1)
if (max(series$med, na.rm=TRUE) > 20) {
  series <- series %>% mutate(cum_AR_pct = med / 100)
} else {
  series <- series %>% mutate(cum_AR_pct = med)
}
final_pct <- series %>% filter(week == max(week, na.rm=TRUE)) %>% pull(cum_AR_pct)
series <- series %>% mutate(remaining_pct = pmax(0, final_pct - cum_AR_pct))

cat('\nWeek \t cumulative_AR_pct \t remaining_pct (percent)\n')
print(series %>% select(week, cum_AR_pct, remaining_pct))

# also write to CSV for inspection
if (!dir.exists('06_Results')) dir.create('06_Results')
write.csv(series, file = '06_Results/tmp_remaining_series.csv', row.names = FALSE)
cat('Wrote 06_Results/tmp_remaining_series.csv\n')
