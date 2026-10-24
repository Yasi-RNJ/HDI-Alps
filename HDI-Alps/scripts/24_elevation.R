###############################################################################
## 24_elevation.R
## Is documentation worse at altitude? The abstract claims high-elevation gaps;
## this is the figure behind that half of the sentence.
## Also writes species_elevation.csv, reused by 23_retention.R.
###############################################################################
suppressPackageStartupMessages({library(data.table); library(ggplot2)})

HERB <- "herbaria_occurrences_filtered.csv"
MET  <- "HDI_output/metrics_complete.csv"
OUT  <- "figures"
dir.create(OUT, showWarnings=FALSE, recursive=TRUE)
MIN_EL <- 5   # specimens with an elevation needed before a species median means anything

cat("=== 24_elevation.R ===\n"); print(Sys.time())
h <- fread(HERB, select=c("species","elevation"))
h[, species := trimws(species)]
cat("specimens:", nrow(h),
    "  with elevation:", h[!is.na(elevation), .N],
    sprintf("(%.0f%%)\n", 100*mean(!is.na(h$elevation))))
he <- h[!is.na(elevation) & elevation > -100 & elevation < 5000]

el <- he[, .(n_elev=.N, elev_median=median(elevation),
             elev_q25=quantile(elevation,.25), elev_q75=quantile(elevation,.75)),
         by=species]
fwrite(el, file.path(OUT,"species_elevation.csv"))
cat("species with any elevation:", nrow(el),
    "  with >=", MIN_EL, ":", el[n_elev >= MIN_EL, .N], "\n")

d <- fread(MET)
if ("in_analysis" %in% names(d)) d <- d[in_analysis == TRUE]
d <- merge(d[!is.na(HDI)], el, by="species", all.x=TRUE)
e <- d[!is.na(elev_median) & n_elev >= MIN_EL]
cat("species in the elevation analysis:", nrow(e), "of", nrow(d), "scorable\n\n")

rho <- cor(e$elev_median, e$HDI, method="spearman")
cat("=== elevation vs HDI ===\n")
cat("Spearman rho:", round(rho,4), "\n")
e[, band := cut(elev_median, c(-100,500,1000,1500,2000,2500,5000),
   labels=c("<500","500-1000","1000-1500","1500-2000","2000-2500",">2500"))]
print(e[, .(species=.N, median_HDI=round(median(HDI),1),
            median_sheets=median(n_sheets),
            median_overlap=round(median(overlap_pct, na.rm=TRUE),1)), by=band][order(band)])
k <- kruskal.test(HDI ~ band, data=e)
cat("\nKruskal-Wallis across bands: chi-sq", round(k$statistic,1),
    " df", k$parameter, " p", format.pval(k$p.value, digits=3), "\n")

th <- theme_minimal(base_size=11) +
  theme(plot.title=element_text(face="bold", colour="#1F3864"),
        plot.subtitle=element_text(colour="grey40"))

p1 <- ggplot(e, aes(elev_median, HDI)) +
  geom_point(alpha=.22, size=.8, colour="#2C7FB8") +
  geom_smooth(method="gam", formula=y~s(x, bs="cs"), colour="#C2544D", fill="#C2544D", alpha=.2) +
  labs(x="Median specimen elevation (m)", y="HDI",
       title="Documentation against altitude",
       subtitle=sprintf("%s species with at least %d specimens carrying an elevation; Spearman rho = %.3f",
                format(nrow(e), big.mark=","), MIN_EL, rho)) + th
ggsave(file.path(OUT,"fig_elevation_hdi.pdf"), p1, width=8, height=5.5)
ggsave(file.path(OUT,"fig_elevation_hdi.png"), p1, width=8, height=5.5, dpi=300)

p2 <- ggplot(e[!is.na(band)], aes(band, HDI, fill=band)) +
  geom_violin(colour="white", linewidth=.3, alpha=.85) +
  geom_boxplot(width=.13, outlier.shape=NA, fill="white", linewidth=.35) +
  scale_fill_viridis_d(option="turbo", guide="none") +
  labs(x="Median specimen elevation (m)", y="HDI",
       title="Documentation by elevation band",
       subtitle=sprintf("Kruskal-Wallis p = %s", format.pval(k$p.value, digits=3))) + th
ggsave(file.path(OUT,"fig_elevation_bands.pdf"), p2, width=8, height=5.5)
ggsave(file.path(OUT,"fig_elevation_bands.png"), p2, width=8, height=5.5, dpi=300)

cat("\nCAUTION FOR THE CAPTION\n")
cat("  Elevation is the median of the specimens that happen to carry one, so\n")
cat("  species with few georeferenced elevations have a noisy value, and any\n")
cat("  bias in WHICH specimens got an elevation carries straight into this\n")
cat("  figure. Report the percentage of specimens with an elevation.\n")
cat("\nwritten to", OUT, "\n"); print(Sys.time())
