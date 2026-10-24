###############################################################################
## 22_use_cases.R
## The paper's practical claim, made visible: which species suit a continuous
## temporal reconstruction, and which only a paired before/after comparison.
##
##   x = historical-modern overlap   (is there a spatial baseline?)
##   y = temporal coverage           (is there a continuous time series?)
###############################################################################
suppressPackageStartupMessages({library(data.table); library(ggplot2)})

MET <- "HDI_output/metrics_complete.csv"
OUT <- "figures"
dir.create(OUT, showWarnings=FALSE, recursive=TRUE)

cat("=== 22_use_cases.R ===\n"); print(Sys.time())
d <- fread(MET)
if ("in_analysis" %in% names(d)) d <- d[in_analysis == TRUE]
d <- d[!is.na(HDI) & !is.na(overlap_pct)]
d[, coverage_pct := 100 * pmin(n_bins_occupied / n_bins_possible, 1)]
cat("species:", nrow(d), "\n")

CUT_OV <- median(d$overlap_pct); CUT_CV <- median(d$coverage_pct)
cat(sprintf("\nquadrant cuts at the medians: overlap %.1f%%, coverage %.1f%%\n",
            CUT_OV, CUT_CV))

d[, quad := fcase(
  coverage_pct >= CUT_CV & overlap_pct >= CUT_OV, "Continuous reconstruction",
  coverage_pct <  CUT_CV & overlap_pct >= CUT_OV, "Paired before / after",
  coverage_pct >= CUT_CV & overlap_pct <  CUT_OV, "Timeline, no spatial anchor",
  default = "Neither")]
LEV <- c("Continuous reconstruction","Paired before / after",
         "Timeline, no spatial anchor","Neither")
d[, quad := factor(quad, levels=LEV)]
COLS <- c("Continuous reconstruction"="#2C7FB8","Paired before / after"="#6BAE75",
          "Timeline, no spatial anchor"="#E8A33A","Neither"="#C2544D")

cat("\n=== how the flora divides ===\n")
tb <- d[, .(species=.N, median_HDI=round(median(HDI),1),
            median_sheets=median(n_sheets)), by=quad][order(quad)]
tb[, pct := round(100*species/sum(species),1)]
print(tb)
fwrite(tb, file.path(OUT,"use_case_counts.csv"))

cat("\n=== sensitivity of the counts to where the cuts are placed ===\n")
for (q in c(.33,.5,.67)) {
  co <- quantile(d$overlap_pct,q); cc <- quantile(d$coverage_pct,q)
  cat(sprintf("  cuts at the %.0fth percentile (overlap %.1f%%, coverage %.1f%%): continuous %d, paired %d\n",
      100*q, co, cc, d[coverage_pct>=cc & overlap_pct>=co, .N],
      d[coverage_pct<cc & overlap_pct>=co, .N]))
}
cat("  The split is a presentational device, not a result. Say so in the caption.\n")

cat("\n=== three examples per quadrant (highest HDI in each) ===\n")
print(d[order(-HDI), head(.SD,3), by=quad,
        .SDcols=c("species","HDI","n_sheets","coverage_pct","overlap_pct")][
        , .(quad, species, HDI=round(HDI,1), n_sheets,
            coverage=round(coverage_pct), overlap=round(overlap_pct,1))])

lab <- d[, .(x=c(90,90,8,8), y=c(97,8,97,8),
             t=LEV)][, .(x,y,t=factor(t,levels=LEV))]
p <- ggplot(d, aes(overlap_pct, coverage_pct)) +
  geom_point(aes(colour=quad), alpha=.30, size=.8) +
  geom_hline(yintercept=CUT_CV, linetype=2, colour="grey35") +
  geom_vline(xintercept=CUT_OV, linetype=2, colour="grey35") +
  geom_label(data=lab, aes(x,y,label=t,colour=t), fill="white",
             label.size=.3, size=3.3, fontface="bold", show.legend=FALSE) +
  scale_colour_manual(values=COLS, guide="none") +
  scale_x_continuous(limits=c(0,100)) + scale_y_continuous(limits=c(0,100)) +
  labs(x="Historical-modern overlap (% of modern cells with a specimen)",
       y="Temporal coverage (% of decades occupied)",
       title="What each species can actually support",
       subtitle=sprintf("%s species; dashed lines at the medians", format(nrow(d), big.mark=","))) +
  theme_minimal(base_size=11) +
  theme(plot.title=element_text(face="bold", colour="#1F3864"),
        plot.subtitle=element_text(colour="grey40"))
ggsave(file.path(OUT,"fig_use_cases.pdf"), p, width=8, height=7.5)
ggsave(file.path(OUT,"fig_use_cases.png"), p, width=8, height=7.5, dpi=300)

fwrite(d[, .(species, HDI, HDI_rank, n_sheets, coverage_pct, overlap_pct, quad)],
       file.path(OUT,"use_case_assignment.csv"))
cat("\nwritten to", OUT, "\n"); print(Sys.time())
