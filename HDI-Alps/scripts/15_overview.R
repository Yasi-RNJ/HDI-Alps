###############################################################################
## 25_overview.R  --  the overview panel figure
##   A  HDI distribution, with data-driven bands
##   B  elevation of the herbarium specimens themselves
##   C  species richness per elevation band   <- the denominator panel
##   D  documentation against richness, per band  <- the synthesis
##
## C and D are the point: a documentation gap at altitude only means something
## once you know how many species and specimens are up there to begin with.
###############################################################################
suppressPackageStartupMessages({
  library(data.table); library(ggplot2)
  if (!requireNamespace("patchwork", quietly=TRUE))
    install.packages("patchwork", repos="https://cloud.r-project.org")
  library(patchwork)
})

HERB <- "herbaria_occurrences_filtered.csv"
MET  <- "HDI_output/metrics_complete.csv"
OUT  <- "figures"
dir.create(OUT, showWarnings=FALSE, recursive=TRUE)
BAND <- 100   # m, for the richness curve

cat("=== 25_overview.R ===\n"); print(Sys.time())

d <- fread(MET)
if ("in_analysis" %in% names(d)) d <- d[in_analysis == TRUE]
sc <- d[!is.na(HDI)]
cat("scorable species:", nrow(sc), "\n")

## ---- bands from the data, not from the old index ------------------------
nz <- sc[HDI > 0]
CUTS <- quantile(nz$HDI, c(1/3, 2/3))
cat(sprintf("\nHDI terciles among non-zero species: %.1f and %.1f\n", CUTS[1], CUTS[2]))
cat("  (the old figure used 10 and 30, from the previous index - the new\n")
cat("   median is", round(median(nz$HDI),1), "so those breaks no longer separate anything)\n")
sc[, band := fcase(HDI == 0, "Unscorable (0)",
                   HDI <  CUTS[1], "Lower third",
                   HDI <  CUTS[2], "Middle third",
                   default = "Upper third")]
BL <- c("Unscorable (0)","Lower third","Middle third","Upper third")
sc[, band := factor(band, levels=BL)]
BC <- c("Unscorable (0)"="grey70","Lower third"="#C2544D",
        "Middle third"="#E8A33A","Upper third"="#2C7FB8")
print(sc[, .(species=.N, median_HDI=round(median(HDI),1)), by=band][order(band)])

th <- theme_minimal(base_size=10) +
  theme(plot.title=element_text(face="bold", colour="#1F3864", size=11),
        plot.subtitle=element_text(colour="grey40", size=9),
        legend.position="bottom", legend.title=element_blank(),
        legend.key.size=unit(.4,"cm"), legend.text=element_text(size=8))

pA <- ggplot(sc, aes(HDI, fill=band)) +
  geom_histogram(bins=55, colour="white", linewidth=.15) +
  geom_vline(xintercept=CUTS, linetype=2, colour="grey25") +
  scale_fill_manual(values=BC) +
  labs(x="HDI", y="Number of species", title="A   Documentation across the flora",
       subtitle=sprintf("%s species; dashed lines at the terciles (%.0f, %.0f)",
                format(nrow(sc), big.mark=","), CUTS[1], CUTS[2])) + th

## ---- elevation ----------------------------------------------------------
h <- fread(HERB, select=c("species","elevation"))
h[, species := trimws(species)]
he <- h[!is.na(elevation) & elevation > -100 & elevation < 5000]
cat("\nspecimens:", nrow(h), "  with a usable elevation:", nrow(he),
    sprintf("(%.0f%%)\n", 100*nrow(he)/nrow(h)))

he[, eband := cut(elevation, c(-100,1000,2000,5000),
                  labels=c("<1000 m","1000-2000 m",">2000 m"))]
EC <- c("<1000 m"="#C2544D","1000-2000 m"="#6BAE75",">2000 m"="#2C7FB8")
cat("\nspecimens per elevation band:\n")
print(he[, .(specimens=.N, species=uniqueN(species)), by=eband][order(eband)])

pB <- ggplot(he, aes(elevation, fill=eband)) +
  geom_histogram(binwidth=BAND, colour="white", linewidth=.1) +
  geom_vline(xintercept=median(he$elevation), colour="grey20") +
  geom_vline(xintercept=quantile(he$elevation, c(.25,.75)), linetype=2, colour="grey35") +
  scale_fill_manual(values=EC) +
  labs(x="Elevation (m)", y="Number of specimens",
       title="B   Where the specimens came from",
       subtitle=sprintf("%s specimens carrying an elevation (%.0f%% of all); solid line = median",
                format(nrow(he), big.mark=","), 100*nrow(he)/nrow(h))) + th

## ---- richness per band --------------------------------------------------
he[, ebin := floor(elevation/BAND)*BAND]
rich <- he[, .(species=uniqueN(species), specimens=.N), by=ebin][order(ebin)]
cat("\nrichness peaks at", rich[which.max(species), ebin], "m with",
    rich[, max(species)], "species\n")
pC <- ggplot(rich, aes(ebin, species)) +
  geom_col(fill="#8EC3E0", width=BAND*.9) +
  labs(x="Elevation (m)", y="Species recorded",
       title="C   How many species are up there at all",
       subtitle=sprintf("Unique species per %d m band - the denominator for any gap claim", BAND)) + th

## ---- D: documentation against richness ---------------------------------
el <- he[, .(n_elev=.N, elev_median=median(elevation)), by=species][n_elev >= 5]
md <- merge(sc, el, by="species")
md[, ebin := floor(elev_median/500)*500]
agg <- md[, .(species=.N, median_HDI=median(HDI)), by=ebin][species >= 10][order(ebin)]
cat("\n=== documentation by 500 m band (>= 10 species) ===\n"); print(agg)
pD <- ggplot(agg, aes(ebin, median_HDI)) +
  geom_line(colour="#C2544D", linewidth=.8) +
  geom_point(aes(size=species), shape=21, fill="#C2544D", colour="grey20") +
  scale_size_continuous(name="species", range=c(2,6)) +
  labs(x="Median specimen elevation (m)", y="Median HDI",
       title="D   Documentation against altitude",
       subtitle="Point size is how many species sit in each band") +
  th + theme(legend.position="right", legend.title=element_text(size=8))

fig <- (pA | pB) / (pC | pD) +
  plot_annotation(title="Herbarium documentation of the Alpine flora",
    theme=theme(plot.title=element_text(face="bold", size=14, colour="#1F3864")))
ggsave(file.path(OUT,"fig_overview.pdf"), fig, width=13, height=10)
ggsave(file.path(OUT,"fig_overview.png"), fig, width=13, height=10, dpi=300)
for (nm in c("A_hdi","B_specimen_elevation","C_richness","D_hdi_by_elevation")) {
  p <- list(pA,pB,pC,pD)[[match(nm, c("A_hdi","B_specimen_elevation","C_richness","D_hdi_by_elevation"))]]
  ggsave(file.path(OUT, paste0("fig_overview_", nm, ".pdf")), p, width=6.5, height=4.5)
}
fwrite(rich, file.path(OUT,"richness_by_elevation.csv"))
fwrite(agg,  file.path(OUT,"hdi_by_elevation_band.csv"))

cat("\n---------------------------------------------------------------\n")
cat("WHY PANEL C MATTERS\n---------------------------------------------------------------\n")
cat("  Species richness falls steeply above the montane belt. So 'fewer\n")
cat("  well-documented species at altitude' and 'fewer species at altitude'\n")
cat("  look the same unless C is shown. Panel D divides one by the other:\n")
cat("  median HDI per band, with the number of species behind each point.\n")
cat("  If D is flat while C declines, the gap is a richness effect, not a\n")
cat("  documentation effect - and that is worth knowing either way.\n")
cat("\nwritten to", OUT, "\n"); print(Sys.time())
