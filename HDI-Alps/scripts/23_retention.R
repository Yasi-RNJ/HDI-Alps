###############################################################################
## 23_retention.R
## If you require HDI above a threshold, what are you left with - and what are
## you systematically losing? One curve per group.
## Needs species_elevation.csv from 24_elevation.R (falls back to computing it).
###############################################################################
suppressPackageStartupMessages({library(data.table); library(ggplot2)})

MET  <- "HDI_output/metrics_complete.csv"
HERB <- "herbaria_occurrences_filtered.csv"
OUT  <- "figures"
dir.create(OUT, showWarnings=FALSE, recursive=TRUE)
HIGH_EL <- 2000

cat("=== 23_retention.R ===\n"); print(Sys.time())
d <- fread(MET)
if ("in_analysis" %in% names(d)) d <- d[in_analysis == TRUE]
d <- d[!is.na(HDI)]
cat("scorable species:", nrow(d), "\n")

elf <- file.path(OUT,"species_elevation.csv")
if (file.exists(elf)) { el <- fread(elf); cat("using", elf, "\n") } else {
  cat("computing species elevations from the herbarium file\n")
  h <- fread(HERB, select=c("species","elevation")); h[, species := trimws(species)]
  el <- h[!is.na(elevation) & elevation > -100 & elevation < 5000,
          .(n_elev=.N, elev_median=median(elevation)), by=species]
  fwrite(el, elf)
}
d <- merge(d, el, by="species", all.x=TRUE)

## ---- define the groups --------------------------------------------------
grps <- list(`All species` = rep(TRUE, nrow(d)))
if ("wcvp_endemic" %in% names(d)) {
  isend <- tolower(as.character(d$wcvp_endemic)) %in% c("endemic","true","yes","1")
  usable <- !is.na(d$wcvp_endemic) & d$wcvp_endemic != "" &
            tolower(d$wcvp_endemic) != "missing data"
  grps[["Alpine endemics"]] <- isend & usable
  grps[["Non-endemics"]]    <- !isend & usable
}
if ("elev_median" %in% names(d)) {
  grps[[paste0("Above ", HIGH_EL, " m")]] <- !is.na(d$elev_median) & d$elev_median >= HIGH_EL & d$n_elev >= 5
  grps[[paste0("Below ", HIGH_EL, " m")]] <- !is.na(d$elev_median) & d$elev_median <  HIGH_EL & d$n_elev >= 5
}
if ("iucn_assessment" %in% names(d)) {
  thr <- grepl("vulnerable|endangered|critically", tolower(as.character(d$iucn_assessment)))
  if (sum(thr) >= 10) grps[["IUCN threatened"]] <- thr
  else cat("NOTE: only", sum(thr), "threatened species - group omitted\n")
}
if ("establishment" %in% names(d)) {
  intro <- tolower(as.character(d$establishment)) == "introduced"
  if (sum(intro) >= 10) grps[["Introduced"]] <- intro
}
cat("\ngroup sizes:\n")
for (g in names(grps)) cat(sprintf("  %-20s %d\n", g, sum(grps[[g]], na.rm=TRUE)))

## ---- the curves ---------------------------------------------------------
THR <- seq(0, max(d$HDI, na.rm=TRUE), length.out = 120)
ret <- rbindlist(lapply(names(grps), function(g) {
  v <- d$HDI[which(grps[[g]])]
  if (!length(v)) return(NULL)
  data.table(group=g, n=length(v), threshold=THR,
             retained=vapply(THR, function(t) 100*mean(v >= t), numeric(1)))
}))
ret[, group := factor(group, levels=names(grps))]
fwrite(ret, file.path(OUT,"retention_curves.csv"))

cat("\n=== percentage retained at selected thresholds ===\n")
key <- c(20,30,40,50,60,70)
tab <- dcast(ret[threshold %in% sapply(key, function(k) THR[which.min(abs(THR-k))])],
             group ~ round(threshold), value.var="retained")
print(tab[, lapply(.SD, function(x) if (is.numeric(x)) round(x,1) else x)])
fwrite(tab, file.path(OUT,"retention_table.csv"))

if (all(c("Alpine endemics","Non-endemics") %in% names(grps))) {
  cat("\n=== the endemic gap, in usable terms ===\n")
  for (k in c(40,50,60)) {
    t0 <- THR[which.min(abs(THR-k))]
    a <- ret[group=="Alpine endemics" & threshold==t0, retained]
    b <- ret[group=="Non-endemics"    & threshold==t0, retained]
    cat(sprintf("  at HDI >= %d : %.1f%% of endemics retained vs %.1f%% of non-endemics (gap %.1f points)\n",
                k, a, b, b-a))
  }
}

PAL <- c("All species"="grey35","Alpine endemics"="#E3752C","Non-endemics"="#2C7FB8",
         "IUCN threatened"="#C2544D","Introduced"="#9B59B6")
PAL[paste0("Above ", HIGH_EL, " m")] <- "#1A8F7A"
PAL[paste0("Below ", HIGH_EL, " m")] <- "#8CC63E"
p <- ggplot(ret, aes(threshold, retained, colour=group)) +
  geom_line(linewidth=.9) +
  scale_colour_manual(values=PAL[levels(ret$group)], name=NULL) +
  scale_y_continuous(limits=c(0,100)) +
  labs(x="HDI threshold", y="% of the group retained",
       title="What survives each HDI threshold",
       subtitle="Where lines diverge, the index is selecting unevenly across groups") +
  theme_minimal(base_size=11) +
  theme(plot.title=element_text(face="bold", colour="#1F3864"),
        plot.subtitle=element_text(colour="grey40"), legend.position="right")
ggsave(file.path(OUT,"fig_retention.pdf"), p, width=8.5, height=5.5)
ggsave(file.path(OUT,"fig_retention.png"), p, width=8.5, height=5.5, dpi=300)

cat("\n---------------------------------------------------------------\n")
cat("HOW TO USE THIS\n---------------------------------------------------------------\n")
cat("  This is the figure that makes the index practical: it tells a reader\n")
cat("  what a given cut-off costs them. Diverging lines quantify the bias -\n")
cat("  the endemic gap stops being an assertion and becomes a number.\n")
cat("\nwritten to", OUT, "\n"); print(Sys.time())
