###############################################################################
## 28_elevation_panels.R
##   B  elevation distribution of the herbarium specimens
##   C  species richness per 100 m elevation band
## Written as standalone panels and as one combined figure.
###############################################################################
suppressPackageStartupMessages({
  library(data.table); library(ggplot2)
  if (!requireNamespace("patchwork", quietly=TRUE))
    install.packages("patchwork", repos="https://cloud.r-project.org")
  library(patchwork)
})

HERB <- "herbaria_occurrences_filtered.csv"
OUT  <- "figures"
dir.create(OUT, showWarnings=FALSE, recursive=TRUE)
BAND <- 100

cat("=== 28_elevation_panels.R ===\n"); print(Sys.time())

h <- fread(HERB, select=c("species","elevation"))
h[, species := trimws(species)]
cat("specimens:", nrow(h), "\n")
cat("with an elevation value:", h[!is.na(elevation), .N],
    sprintf("(%.1f%%)\n", 100*mean(!is.na(h$elevation))))

he <- h[!is.na(elevation) & elevation > -100 & elevation < 5000]
cat("usable (-100 to 5000 m):", nrow(he),
    sprintf("(%.1f%% of all specimens)\n", 100*nrow(he)/nrow(h)))
cat("species represented:", uniqueN(he$species), "\n\n")

print(summary(he$elevation))

he[, eband := cut(elevation, c(-100, 1000, 2000, 5000),
                  labels=c("< 1000 m", "1000-2000 m", "> 2000 m"))]
cat("\n=== specimens per elevation band ===\n")
tb <- he[, .(specimens=.N, species=uniqueN(species)), by=eband][order(eband)]
tb[, pct := round(100*specimens/sum(specimens), 1)]
print(tb)

he[, ebin := floor(elevation/BAND)*BAND]
rich <- he[, .(species=uniqueN(species), specimens=.N), by=ebin][order(ebin)]
cat("\n=== richness peak ===\n")
cat("peaks at", rich[which.max(species), ebin], "-",
    rich[which.max(species), ebin]+BAND, "m with", rich[, max(species)], "species\n")
cat("median specimen elevation:", median(he$elevation), "m\n")
cat("quartiles:", paste(quantile(he$elevation, c(.25,.75)), collapse=" - "), "m\n")

th <- theme_minimal(base_size=11) +
  theme(plot.title=element_text(face="bold", colour="#1F3864", size=12),
        plot.subtitle=element_text(colour="grey40", size=9),
        panel.grid.minor=element_blank(),
        legend.position="top", legend.title=element_blank(),
        legend.key.size=unit(.4,"cm"), legend.text=element_text(size=9))

EC <- c("< 1000 m"="#E8857A", "1000-2000 m"="#7FBF7F", "> 2000 m"="#7FA8D9")

pB <- ggplot(he, aes(elevation, fill=eband)) +
  geom_histogram(binwidth=BAND, colour="white", linewidth=.15) +
  geom_vline(xintercept=median(he$elevation), colour="grey15", linewidth=.6) +
  geom_vline(xintercept=quantile(he$elevation, c(.25,.75)),
             linetype=2, colour="grey30", linewidth=.5) +
  scale_fill_manual(values=EC) +
  scale_x_continuous(limits=c(0, 4000), breaks=seq(0, 4000, 1000)) +
  labs(x="Elevation (m)", y="Number of specimens",
       title="Elevation distribution of the herbarium records",
       subtitle=sprintf("%s specimens carrying an elevation (%.0f%% of all); solid line = median, dashed = quartiles",
                format(nrow(he), big.mark=","), 100*nrow(he)/nrow(h))) + th

pC <- ggplot(rich, aes(ebin, species)) +
  geom_col(fill="#8EB8E0", width=BAND*.92) +
  scale_x_continuous(limits=c(0, 4000), breaks=seq(0, 4000, 1000)) +
  labs(x="Elevation (m)", y="Number of species (unique)",
       title="Species richness per elevation",
       subtitle=sprintf("unique species per %d m band; peak %d species at %d-%d m",
                BAND, rich[, max(species)], rich[which.max(species), ebin],
                rich[which.max(species), ebin]+BAND)) +
  th + theme(legend.position="none")

ggsave(file.path(OUT,"fig_elevation_specimens.pdf"), pB, width=7, height=5)
ggsave(file.path(OUT,"fig_elevation_specimens.png"), pB, width=7, height=5, dpi=300)
ggsave(file.path(OUT,"fig_elevation_richness.pdf"), pC, width=7, height=5)
ggsave(file.path(OUT,"fig_elevation_richness.png"), pC, width=7, height=5, dpi=300)

fig <- (pB | pC)
ggsave(file.path(OUT,"fig_elevation_panels.pdf"), fig, width=13, height=5)
ggsave(file.path(OUT,"fig_elevation_panels.png"), fig, width=13, height=5, dpi=300)

fwrite(rich, file.path(OUT,"richness_by_elevation.csv"))
fwrite(tb,   file.path(OUT,"specimens_by_elevation_band.csv"))

cat("\n---------------------------------------------------------------\n")
cat("WHY THE SECOND PANEL MATTERS\n")
cat("---------------------------------------------------------------\n")
cat("  Richness falls steeply above the montane belt. So 'fewer well-documented\n")
cat("  species at altitude' and 'fewer species at altitude' look identical\n")
cat("  unless this panel is shown. Any claim about a high-elevation\n")
cat("  documentation gap needs it as the denominator.\n\n")
cat("  Caveat for the caption: elevation is present on only",
    sprintf("%.0f%%", 100*nrow(he)/nrow(h)), "of specimens,\n")
cat("  and if it was recorded more often for some kinds of collection than\n")
cat("  others, that bias carries into both panels. State the percentage.\n")
cat("\nwritten to", OUT, "\n"); print(Sys.time())
