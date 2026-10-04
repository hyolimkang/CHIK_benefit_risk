suppressMessages({
  library(dplyr)
})

rdata <- "01_Data/draw_level_xy_serostatus_finite.RData"
if (!file.exists(rdata)) stop("RData not found: ", rdata)
loaded <- load(rdata)
cat("Loaded objects:", paste(loaded, collapse=", "), "\n")

if (!exists("draw_level_xy_serostatus")) stop("draw_level_xy_serostatus not found in RData")
obj <- draw_level_xy_serostatus
cat("Rows, cols:", nrow(obj), ncol(obj), "\n")
cat("Columns:\n")
print(names(obj))

# look for inf or AR columns
cols <- names(obj)
candidates <- cols[grepl('inf|AR|attack|symp', cols, ignore.case = TRUE)]
cat('Candidates:', paste(candidates, collapse=', '), '\n')

# prefer inf_10k or symp_10k
prefer <- c('inf_10k','symp_10k','symp_nv_10k','AR_travel')
chosen <- intersect(prefer, cols)[1]
if (is.na(chosen) || is.null(chosen)) chosen <- candidates[1]
cat('Chosen:', chosen, '\n')

# ensure Region & week
if (!('Region' %in% cols)) stop('no Region')
if (!('week' %in% cols)) stop('no week')

# compute final median per region
finals <- obj %>% filter(week == max(week, na.rm=TRUE)) %>% group_by(Region) %>% summarise(final_med = median(.data[[chosen]], na.rm=TRUE)) %>% arrange(desc(final_med))
print(head(finals,6))

top_region <- finals$Region[1]
cat('Top region:', top_region, '\n')

series <- obj %>% filter(Region == top_region) %>% group_by(week) %>% summarise(med = median(.data[[chosen]], na.rm=TRUE)) %>% arrange(week)

# If chosen looks like per10k (values > 20), convert to percent by /100; if values in [0,1] maybe it's fraction
if (max(series$med, na.rm=TRUE) > 20) {
  series <- series %>% mutate(cum_AR_pct = med / 100)
} else if (max(series$med, na.rm=TRUE) <= 1) {
  series <- series %>% mutate(cum_AR_pct = med * 100)
} else {
  series <- series %>% mutate(cum_AR_pct = med)
}
final_pct <- series$cum_AR_pct[which.max(series$week)]
series <- series %>% mutate(remaining_pct = pmax(0, final_pct - cum_AR_pct))

cat('\nWeek \t cumulative_AR_pct \t remaining_pct\n')
print(series %>% select(week, cum_AR_pct, remaining_pct))
write.csv(series, '06_Results/tmp_remaining_series2.csv', row.names = FALSE)
cat('Wrote 06_Results/tmp_remaining_series2.csv\n')
