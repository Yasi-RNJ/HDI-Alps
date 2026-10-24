###############################################################################
## 21_maps.R  --  three maps on the same 10 km grid used by the metrics
##   1. collecting effort: specimens per cell
##   2. the same, split by era
##   3. per-species: historical only / both / modern only   <- the overlap term
###############################################################################
suppressPackageStartupMessages({library(data.table); library(sf); library(ggplot2)})

HERB <- "herbaria_occurrences_filtered.csv"
MOD  <- "current_occurrences_enriched.csv"
MET  <- "HDI_output/metrics_complete.csv"
OUT  <- "figures"
dir.create(OUT, showWarnings=FALSE, recursive=TRUE)

GRID_M <- 10000; EPSG <- 3035
HERB_WIN <- c(1740, 2019); MOD_WIN <- c(2000, 2019)

cat("=== 21_maps.R ===\n"); print(Sys.time())

## ---- projection, chunked, with an axis-order guard ---------------------
to_xy <- function(lon, lat) {
  keep <- !is.na(lon) & !is.na(lat)
  X <- rep(NA_real_, length(lon)); Y <- X
  idx <- which(keep); ch <- 2e6L
  for (s in seq(1L, length(idx), by=ch)) {
    i <- idx[s:min(s+ch-1L, length(idx))]
    p <- sf::sf_project("EPSG:4326", paste0("EPSG:", EPSG), cbind(lon[i], lat[i]))
    X[i] <- p[,1]; Y[i] <- p[,2]
  }
  list(cx = floor(X/GRID_M)*GRID_M, cy = floor(Y/GRID_M)*GRID_M)
}
chk <- sf::sf_project("EPSG:4326", paste0("EPSG:", EPSG), cbind(11, 46))
cat("axis check (11E,46N) -> X", round(chk[1]), " Y", round(chk[2]), "\n")
if (abs(chk[1]-4398558) > 5000 || abs(chk[2]-2543530) > 5000)
  stop("projection axis order wrong - every map would be plausible nonsense")

cat("\nreading herbarium records...\n")
h <- fread(HERB, select=c("species","year","decimalLatitude","decimalLongitude"))
h[, species := trimws(species)]
h <- h[!is.na(decimalLatitude) & !is.na(decimalLongitude)]
hc <- to_xy(h$decimalLongitude, h$decimalLatitude)
h[, `:=`(cx=hc$cx, cy=hc$cy)]; h <- h[!is.na(cx)]
cat("  mapped:", nrow(h), " cells:", uniqueN(h[, .(cx,cy)]), "\n")

cat("reading modern records (the big file)...\n")
mo <- fread(MOD, select=c("species","year","decimalLatitude","decimalLongitude"))
mo[, species := trimws(species)]
mo <- mo[!is.na(decimalLatitude) & !is.na(decimalLongitude) &
         year >= MOD_WIN[1] & year <= MOD_WIN[2]]
mc <- to_xy(mo$decimalLongitude, mo$decimalLatitude)
mo[, `:=`(cx=mc$cx, cy=mc$cy)]; mo <- mo[!is.na(cx)]
cat("  mapped:", nrow(mo), " cells:", uniqueN(mo[, .(cx,cy)]), "\n\n")

base_th <- theme_minimal(base_size=11) +
  theme(panel.grid=element_blank(), axis.text=element_blank(),
        axis.title=element_blank(),
        plot.title=element_text(face="bold", colour="#1F3864"),
        plot.subtitle=element_text(colour="grey40"))

## ---- 1. collecting effort ----------------------------------------------
e1 <- h[, .(sheets=.N), by=.(cx,cy)]
cat("=== map 1: collecting effort ===\n"); print(summary(e1$sheets))
p1 <- ggplot(e1, aes(cx, cy, fill=sheets)) + geom_raster() + coord_equal() +
  scale_fill_viridis_c(option="magma", trans="log10", direction=-1,
                       name="Specimens\n(log scale)") +
  labs(title="Herbarium collection effort across the Alps",
       subtitle=sprintf("%s specimens in %s cells of 10 x 10 km (ETRS89-LAEA)",
                format(nrow(h), big.mark=","), format(nrow(e1), big.mark=","))) + base_th
ggsave(file.path(OUT,"fig_map_effort.pdf"), p1, width=10, height=6.5)
ggsave(file.path(OUT,"fig_map_effort.png"), p1, width=10, height=6.5, dpi=300)

## ---- 2. effort by era ---------------------------------------------------
he <- h[!is.na(year) & year >= HERB_WIN[1] & year <= HERB_WIN[2]]
he[, era := cut(year, c(1739,1849,1899,1949,1999,2019),
                labels=c("1740-1849","1850-1899","1900-1949","1950-1999","2000-2019"))]
e2 <- he[!is.na(era), .(sheets=.N), by=.(cx,cy,era)]
cat("\n=== map 2: specimens per era ===\n")
print(he[!is.na(era), .(sheets=.N, cells=uniqueN(paste(cx,cy))), by=era][order(era)])
p2 <- ggplot(e2, aes(cx, cy, fill=sheets)) + geom_raster() + coord_equal() +
  facet_wrap(~era, nrow=2) +
  scale_fill_viridis_c(option="magma", trans="log10", direction=-1, name="Specimens") +
  labs(title="Where collecting happened, by era",
       subtitle="The same grid in every panel") + base_th
ggsave(file.path(OUT,"fig_map_effort_by_era.pdf"), p2, width=13, height=7.5)
ggsave(file.path(OUT,"fig_map_effort_by_era.png"), p2, width=13, height=7.5, dpi=300)

## ---- 3. overall cell coverage ------------------------------------------
hcell <- unique(h[, .(cx,cy)])[, hist := TRUE]
mcell <- unique(mo[, .(cx,cy)])[, modn := TRUE]
ov <- merge(hcell, mcell, by=c("cx","cy"), all=TRUE)
ov[is.na(hist), hist := FALSE][is.na(modn), modn := FALSE]
ov[, status := fcase(hist & modn, "Both periods",
                     hist & !modn, "Historical only",
                     !hist & modn, "Modern only")]
ov[, status := factor(status, levels=c("Historical only","Both periods","Modern only"))]
tb <- ov[, .(cells=.N), by=status][order(status)][, pct := round(100*cells/sum(cells),1)][]
cat("\n=== map 3: cell coverage by era ===\n"); print(tb)
COLS <- c("Historical only"="#D6407F","Both periods"="#8CC63E","Modern only"="#2C8FC9")
p3 <- ggplot(ov, aes(cx, cy, fill=status)) + geom_raster() + coord_equal() +
  scale_fill_manual(values=COLS, name=NULL) +
  labs(title="Which cells have historical records, modern records, or both",
       subtitle=sprintf("%s cells; %s%% hold records from both periods",
                format(nrow(ov), big.mark=","), tb[status=="Both periods", pct])) +
  base_th + theme(legend.position="bottom")
ggsave(file.path(OUT,"fig_map_overlap.pdf"), p3, width=10, height=7)
ggsave(file.path(OUT,"fig_map_overlap.png"), p3, width=10, height=7, dpi=300)
fwrite(ov, file.path(OUT,"cell_era_status.csv"))

## ---- 4. the same, per species -------------------------------------------
m <- fread(MET)
if ("in_analysis" %in% names(m)) m <- m[in_analysis == TRUE]
sc <- m[!is.na(HDI) & HDI > 0][order(HDI_rank)]
FOCAL <- unique(c(sc[1:2, species],
                  sc[floor(.N/2), species],
                  sc[.N, species],
                  "Larix decidua"))
FOCAL <- FOCAL[FOCAL %in% sc$species][1:4]
cat("\n=== per-species overlap panels ===\n")
print(sc[species %in% FOCAL, .(HDI_rank, species, HDI=round(HDI,1),
                               overlap=round(overlap_pct,1), n_sheets)])

sp <- rbindlist(lapply(FOCAL, function(s) {
  a <- unique(h [species==s, .(cx,cy)])[, hist := TRUE]
  b <- unique(mo[species==s, .(cx,cy)])[, modn := TRUE]
  if (!nrow(a) && !nrow(b)) return(NULL)
  z <- merge(a, b, by=c("cx","cy"), all=TRUE)
  z[is.na(hist), hist := FALSE][is.na(modn), modn := FALSE]
  z[, status := fcase(hist & modn, "Both periods",
                      hist & !modn, "Historical only",
                      !hist & modn, "Modern only")]
  info <- sc[species==s]
  z[, panel := sprintf("%s\nHDI %.1f  |  overlap %.1f%%", s, info$HDI, info$overlap_pct)]
  z
}))
sp[, status := factor(status, levels=c("Historical only","Both periods","Modern only"))]
p4 <- ggplot(sp, aes(cx, cy, fill=status)) + geom_raster() + coord_equal() +
  facet_wrap(~panel, ncol=2) +
  scale_fill_manual(values=COLS, name=NULL) +
  labs(title="Historical and contemporary occupancy",
       subtitle=paste("Green cells hold records from both periods and are what the overlap term counts;",
                      "blue cells are places the species is recorded today but was never collected", sep="\n")) +
  base_th + theme(legend.position="bottom",
                  strip.text=element_text(face="bold", size=9, lineheight=1.2))
ggsave(file.path(OUT,"fig_map_overlap_species.pdf"), p4, width=13, height=10)
ggsave(file.path(OUT,"fig_map_overlap_species.png"), p4, width=13, height=10, dpi=300)

cat("\n---------------------------------------------------------------\n")
cat("The species panels are the direct answer to referee 1's ninth point.\n")
cat("Pink cells were collected historically and have no modern record; blue\n")
cat("cells are newly surveyed ground with no historical baseline. Only the\n")
cat("green cells can support a before-and-after comparison at all.\n")
cat("\nwritten to", OUT, "\n"); print(Sys.time())
