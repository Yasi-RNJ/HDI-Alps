###############################################################################
## 20_ridges.R  --  collection-history ridges at three points on the index
##
## THREE GROUPS, chosen BEFORE the records are filtered:
##   top 10     - highest HDI
##   middle 10  - straddling the median HDI rank
##   bottom 10  - lowest HDI that still has enough dated records to draw
##
## Produces, for each group and for all three together:
##   year ridges (one colour per species)
##   year ridges coloured by period
##   elevation ridges
###############################################################################
suppressPackageStartupMessages({
  library(data.table); library(ggplot2)
  if (!requireNamespace("ggridges", quietly=TRUE))
    install.packages("ggridges", repos="https://cloud.r-project.org")
  library(ggridges)
})

HERB <- "herbaria_occurrences_filtered.csv"
MET  <- "HDI_output/metrics_complete.csv"
OUT  <- "figures"
dir.create(OUT, showWarnings=FALSE, recursive=TRUE)

WIN     <- c(1740, 2019)
N_EACH  <- 10
MIN_REC <- 20     # a density needs enough points to mean anything

cat("=== 20_ridges.R ===\n"); print(Sys.time())

m <- fread(MET)
if ("in_analysis" %in% names(m)) m <- m[in_analysis == TRUE]
sc <- m[!is.na(HDI) & HDI > 0]
setorder(sc, HDI_rank)
cat("scorable species:", nrow(sc), "\n")

## ---- how many DATED records does each species have? --------------------
h <- fread(HERB, select=c("species","year","elevation"))
h[, species := trimws(species)]
hd <- h[!is.na(year) & year >= WIN[1] & year <= WIN[2]]
cnt <- hd[, .(n_dated = .N), by=species]
sc <- merge(sc, cnt, by="species", all.x=TRUE)
sc[is.na(n_dated), n_dated := 0L]
setorder(sc, HDI_rank)

## ---- pick the three groups FIRST ---------------------------------------
ok  <- sc[n_dated >= MIN_REC]
cat("species with >=", MIN_REC, "dated records:", nrow(ok), "\n\n")

top <- ok[order(HDI_rank)][1:N_EACH]
bot <- ok[order(-HDI_rank)][1:N_EACH]
mid_i <- floor(nrow(ok)/2) - floor(N_EACH/2) + seq_len(N_EACH)
mid <- ok[order(HDI_rank)][mid_i]

stopifnot(nrow(top)==N_EACH, nrow(mid)==N_EACH, nrow(bot)==N_EACH,
          length(unique(c(top$species, mid$species, bot$species))) == 3*N_EACH)

GR <- rbindlist(list(
  top[, .(species, HDI, HDI_rank, n_sheets, n_dated, grp="Best documented")],
  mid[, .(species, HDI, HDI_rank, n_sheets, n_dated, grp="Middle of the index")],
  bot[, .(species, HDI, HDI_rank, n_sheets, n_dated, grp="Least documented")]))
GR[, grp := factor(grp, levels=c("Best documented","Middle of the index","Least documented"))]
GR[, lab := sprintf("%s   HDI %.1f", species, HDI)]
cat("=== the three groups ===\n")
print(GR[, .(grp, HDI_rank, species, HDI=round(HDI,1), n_sheets, n_dated)])
fwrite(GR, file.path(OUT,"ridge_groups.csv"))

## ---- NOW filter the records, once --------------------------------------
hd <- hd[species %in% GR$species]
hd <- merge(hd, GR[, .(species, grp, lab, HDI_rank)], by="species")
cat("\ndated records retained:", nrow(hd), "for", uniqueN(hd$species), "species\n")
stopifnot(uniqueN(hd$species) == 3*N_EACH)

ord <- GR[order(grp, HDI_rank), lab]
hd[, lab := factor(lab, levels=rev(ord))]

th <- theme_ridges(font_size=11, grid=TRUE) +
  theme(axis.text.y = element_text(face="italic", size=9),
        plot.title = element_text(face="bold", colour="#1F3864"),
        plot.subtitle = element_text(colour="grey40"),
        strip.text.y = element_text(angle=0, face="bold", size=9))

## =========================================================================
## 1. one figure per group: year
## =========================================================================
pal <- c("#E3572C","#F2A33A","#F5D33A","#8CC63E","#2FB49C",
         "#2C8FC9","#5A5FC7","#9B4FC7","#D6407F","#E8833A")
for (g in levels(GR$grp)) {
  dd <- hd[grp == g]
  dd[, lab := droplevels(lab)]
  p <- ggplot(dd, aes(year, lab, fill=lab)) +
    geom_density_ridges(scale=2.4, alpha=.9, colour="white", linewidth=.3,
                        rel_min_height=.008) +
    scale_fill_manual(values=rev(pal), guide="none") +
    scale_x_continuous(limits=WIN, breaks=seq(1750,2000,50)) +
    labs(x="Year of collection", y=NULL,
         title=sprintf("Collection dates: %s", tolower(g)),
         subtitle="Each ridge is the density of dated specimens through time") + th
  f <- gsub("[^a-z]+","_", tolower(g))
  ggsave(file.path(OUT, paste0("fig_ridge_year_", f, ".pdf")), p, width=9, height=7)
  ggsave(file.path(OUT, paste0("fig_ridge_year_", f, ".png")), p, width=9, height=7, dpi=300)
  cat("  wrote fig_ridge_year_", f, "\n", sep="")
}

## =========================================================================
## 2. all three groups together - the contrast figure
## =========================================================================
p_all <- ggplot(hd, aes(year, lab, fill=grp)) +
  geom_density_ridges(scale=1.7, alpha=.85, colour="white", linewidth=.3,
                      rel_min_height=.008) +
  scale_fill_manual(values=c("Best documented"="#2C7FB8",
                             "Middle of the index"="#6BAE75",
                             "Least documented"="#E3752C"), name=NULL) +
  scale_x_continuous(limits=WIN, breaks=seq(1750,2000,50)) +
  facet_grid(grp ~ ., scales="free_y", space="free_y") +
  labs(x="Year of collection", y=NULL,
       title="Collection histories across the index",
       subtitle=sprintf("Ten species at each of three positions; all have at least %d dated specimens", MIN_REC)) +
  th + theme(legend.position="none")
ggsave(file.path(OUT,"fig_ridge_year_three_groups.pdf"), p_all, width=10, height=12)
ggsave(file.path(OUT,"fig_ridge_year_three_groups.png"), p_all, width=10, height=12, dpi=300)

## =========================================================================
## 3. coloured by period, all three groups
## =========================================================================
p_grad <- ggplot(hd, aes(year, lab, fill=after_stat(x))) +
  geom_density_ridges_gradient(scale=1.7, colour="white", linewidth=.3,
                               rel_min_height=.008) +
  scale_fill_viridis_c(option="turbo", name="Year") +
  scale_x_continuous(limits=WIN, breaks=seq(1750,2000,50)) +
  facet_grid(grp ~ ., scales="free_y", space="free_y") +
  labs(x="Year of collection", y=NULL,
       title="Collection dates, coloured by period",
       subtitle="Colour tracks the year rather than the species") + th
ggsave(file.path(OUT,"fig_ridge_gradient_three_groups.pdf"), p_grad, width=10, height=12)
ggsave(file.path(OUT,"fig_ridge_gradient_three_groups.png"), p_grad, width=10, height=12, dpi=300)

## =========================================================================
## 4. elevation
## =========================================================================
if ("elevation" %in% names(h)) {
  he <- h[species %in% GR$species & !is.na(elevation) &
          elevation > -100 & elevation < 5000]
  he <- merge(he, GR[, .(species, grp, lab, HDI_rank)], by="species")
  he[, lab := factor(lab, levels=rev(ord))]
  pct <- round(100 * nrow(he) / nrow(h[species %in% GR$species]), 0)
  cat("\nspecimens with a usable elevation:", nrow(he),
      sprintf("(%d%% of this subset)\n", pct))
  drop <- setdiff(GR$species, unique(he$species))
  if (length(drop)) cat("  no elevation at all for:", paste(drop, collapse=", "), "\n")

  p_el <- ggplot(he, aes(elevation, lab, fill=after_stat(x))) +
    geom_density_ridges_gradient(scale=1.7, colour="white", linewidth=.3,
                                 rel_min_height=.008) +
    scale_fill_viridis_c(option="turbo", name="m a.s.l.") +
    facet_grid(grp ~ ., scales="free_y", space="free_y") +
    labs(x="Elevation (m)", y=NULL,
         title="Elevational distribution of herbarium specimens",
         subtitle=sprintf("Only %d%% of specimens for these species carry an elevation", pct)) + th
  ggsave(file.path(OUT,"fig_ridge_elevation_three_groups.pdf"), p_el, width=10, height=12)
  ggsave(file.path(OUT,"fig_ridge_elevation_three_groups.png"), p_el, width=10, height=12, dpi=300)
} else cat("\nNO elevation column in the herbarium file - elevation figures skipped\n")

cat("\n---------------------------------------------------------------\n")
cat("WHAT THESE SHOW\n---------------------------------------------------------------\n")
cat("  The three groups are chosen before any record filtering, so the\n")
cat("  bottom group actually appears - the earlier version filtered to the\n")
cat("  top species first and the contrast panel was left with only one end.\n\n")
cat("  All species shown have at least", MIN_REC, "dated specimens. The least\n")
cat("  documented species overall have one or two records and would draw as a\n")
cat("  single spike, which shows nothing. Restricting to species that CAN be\n")
cat("  drawn means the comparison is about the SHAPE of the collection history,\n")
cat("  not about how many records there are. State that in the caption.\n")
cat("\nwritten to", OUT, "\n"); print(Sys.time())
