# ----------------------------------------------------------------------------
#  Figure 3  -  Documentation of the Alpine flora   (four panels)
#  v2 fixes:
#    1. panels B-D computed on the ANALYSIS set only (checklist species,
#       records with coordinates), not the raw 487,487-row download
#    2. panel A legend reads "Scored zero", not "Unscorable (0)"
#    3. elevation axes stop at 4,000 m
#    4. terciles printed to one decimal (50.1 / 60.9)
#    5. colourblind-safe palette (Okabe-Ito)  -> answers Referee 1
# ----------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork)
})

## ---- paths -----------------------------------------------------------------
BASE     <- "/data/fdalgrande_data/HDI-YASI"
METRICS  <- file.path(BASE, "01-build-metrics/baseline_grid10km_bin10yr/metrics_complete.csv")
HERB     <- file.path(BASE, "herbaria_occurrences_filtered.csv")
OUTDIR   <- file.path(BASE, "04-figures")
ELEV_MAX <- 4000          # axis ceiling, metres
BAND     <- 100           # richness band width, metres

## ---- colourblind-safe palette (Okabe & Ito) --------------------------------
OI <- c(blue = "#0072B2", orange = "#E69F00", green = "#009E73",
        grey = "#999999", vermillion = "#D55E00")

## ---- 1. metrics ------------------------------------------------------------
m <- fread(METRICS)
checklist <- unique(m$species)
s <- m[in_analysis == TRUE & !is.na(HDI)]
cat("scored species            :", nrow(s), "\n")

nz  <- s[HDI > 0]
t1  <- quantile(nz$HDI, 1/3, names = FALSE)
t2  <- quantile(nz$HDI, 2/3, names = FALSE)
cat(sprintf("terciles of the non-zero scores: %.2f and %.2f\n", t1, t2))
cat("species scoring exactly zero:", sum(s$HDI == 0), "\n\n")

## ---- 2. herbarium records, filtered to the analysis set --------------------
h <- fread(HERB, fill = TRUE, showProgress = FALSE,
           select = c("species", "decimalLatitude", "decimalLongitude", "elevation"))
cat("rows in the raw download  :", nrow(h), "\n")
h <- h[species %in% checklist &
       !is.na(decimalLatitude) & !is.na(decimalLongitude)]
cat("after checklist + coords  :", nrow(h),
    "   <-- should equal 426,830\n")
if (nrow(h) != 426830)
  cat("  !! does not match the Methods figure - check the filter before using\n")

he <- h[!is.na(elevation) & elevation >= 0 & elevation <= ELEV_MAX]
cat("with a usable elevation   :", nrow(he),
    sprintf(" (%.1f%% of the analysis set)\n", 100 * nrow(he) / nrow(h)))
cat("dropped for elevation > ", ELEV_MAX, " m : ",
    sum(!is.na(h$elevation) & h$elevation > ELEV_MAX), "\n", sep = "")

## ---- 3. numbers the Results text needs -------------------------------------
he[, band3 := fifelse(elevation < 1000, "<1,000 m",
              fifelse(elevation < 2000, "1,000-2,000 m", ">=2,000 m"))]
he[, band3 := factor(band3, levels = c("<1,000 m", "1,000-2,000 m", ">=2,000 m"))]
cat("\nSPECIMENS BY ELEVATION BAND (for the Results text)\n")
print(he[, .(specimens = .N, pct = round(100 * .N / nrow(he), 1)), by = band3][order(band3)])

he[, bandN := floor(elevation / BAND) * BAND]
rich <- he[, .(n_species = uniqueN(species)), by = bandN][order(bandN)]
cat("\nrichness peaks at", rich[which.max(n_species), bandN], "m with",
    rich[which.max(n_species), n_species], "species\n")
cat("(maximum possible is", length(checklist), "- if this is close to that, re-check the filter)\n")

sp_elev <- he[, .(med_elev = median(elevation), n_spec = .N), by = species]
d <- merge(sp_elev, s[, .(species, HDI)], by = "species")
d[, eband := cut(med_elev, breaks = seq(0, ELEV_MAX, by = 500), include.lowest = TRUE)]
dd <- d[, .(med_HDI = median(HDI), n_species = .N,
            mid = (as.numeric(eband) - 0.5) * 500), by = eband][!is.na(eband)][order(mid)]
cat("\nMEDIAN HDI BY 500 m BAND (for the Results text)\n"); print(dd)

## ---- 4. panels -------------------------------------------------------------
th <- theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank(),
        plot.title    = element_text(face = "bold", size = 11),
        plot.subtitle = element_text(size = 8.5, colour = "grey30"),
        legend.position = "bottom", legend.title = element_blank())

s[, tercile := fifelse(HDI == 0, "Scored zero",
                fifelse(HDI <= t1, "Lower third",
                 fifelse(HDI <= t2, "Middle third", "Upper third")))]
s[, tercile := factor(tercile, levels = c("Scored zero", "Lower third",
                                          "Middle third", "Upper third"))]

pA <- ggplot(s, aes(HDI, fill = tercile)) +
  geom_histogram(binwidth = 2, boundary = 0, colour = NA) +
  geom_vline(xintercept = c(t1, t2), linetype = "dashed", linewidth = 0.4) +
  scale_fill_manual(values = c("Scored zero"  = OI[["grey"]],
                               "Lower third"  = OI[["vermillion"]],
                               "Middle third" = OI[["orange"]],
                               "Upper third"  = OI[["blue"]])) +
  labs(title = "A  Documentation across the flora",
       subtitle = sprintf("%s species; dashed lines at the terciles of the non-zero scores (%.1f and %.1f)",
                          format(nrow(s), big.mark = ","), t1, t2),
       x = "HDI", y = "Number of species") + th

pB <- ggplot(he, aes(elevation, fill = band3)) +
  geom_histogram(binwidth = 100, boundary = 0, colour = NA) +
  geom_vline(xintercept = median(he$elevation), linewidth = 0.4) +
  scale_fill_manual(values = c(OI[["vermillion"]], OI[["green"]], OI[["blue"]])) +
  scale_x_continuous(limits = c(0, ELEV_MAX)) +
  labs(title = "B  Where the specimens came from",
       subtitle = sprintf("%s specimens of analysed species carrying an elevation; line = median (%d m)",
                          format(nrow(he), big.mark = ","), round(median(he$elevation))),
       x = "Elevation (m)", y = "Number of specimens") + th

pC <- ggplot(rich, aes(bandN, n_species)) +
  geom_col(fill = OI[["blue"]], width = BAND * 0.9) +
  scale_x_continuous(limits = c(0, ELEV_MAX)) +
  labs(title = "C  How many species occur at each elevation",
       subtitle = sprintf("Unique species per %d m band, of %s on the checklist",
                          BAND, format(length(checklist), big.mark = ",")),
       x = "Elevation (m)", y = "Species recorded") + th

pD <- ggplot(dd, aes(mid, med_HDI)) +
  geom_line(colour = OI[["vermillion"]], linewidth = 0.6) +
  geom_point(aes(size = n_species), colour = OI[["vermillion"]]) +
  scale_x_continuous(limits = c(0, ELEV_MAX)) +
  scale_size_continuous(name = "species", range = c(2, 7)) +
  labs(title = "D  Documentation against elevation",
       subtitle = "Median HDI by the median specimen elevation of a species",
       x = "Median specimen elevation (m)", y = "Median HDI") +
  th + theme(legend.position = "right", legend.title = element_text(size = 8))

fig <- (pA | pB) / (pC | pD) +
  plot_annotation(title = "Herbarium documentation of the Alpine flora",
                  theme = theme(plot.title = element_text(face = "bold", size = 13)))

dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
ggsave(file.path(OUTDIR, "fig_overview_v2.pdf"), fig, width = 13, height = 10,
       device = cairo_pdf, bg = "white")
ggsave(file.path(OUTDIR, "fig_overview_v2.png"), fig, width = 13, height = 10,
       dpi = 300, bg = "white")
cat("\nwritten: fig_overview_v2.pdf and .png\n")
