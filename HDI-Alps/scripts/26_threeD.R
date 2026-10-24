###############################################################################
## 26_threeD.R  --  3D scatter, two versions
##
## VERSION 1 ("as drawn"): HDI on one axis, two of its own inputs on the others.
##   Honest warning: temporal evenness and overlap are INPUTS to HDI, so the
##   shape of this cloud is partly the arithmetic of the formula rather than a
##   finding. Usable as an illustration; do not describe it as a result.
##
## VERSION 2 (preferred): three inputs as the axes, HDI as the colour. Position
##   and colour are then independent, and the figure shows how the index maps
##   onto the input space - which is a finding.
##
## Generous margins are left for the species photographs.
###############################################################################
suppressPackageStartupMessages({library(data.table)})
if (!requireNamespace("plot3D", quietly=TRUE))
  install.packages("plot3D", repos="https://cloud.r-project.org")
library(plot3D)

MET <- "HDI_output/metrics_complete.csv"
OUT <- "figures"
dir.create(OUT, showWarnings=FALSE, recursive=TRUE)
THETA <- 40; PHI <- 18

cat("=== 26_threeD.R ===\n"); print(Sys.time())
d <- fread(MET)
if ("in_analysis" %in% names(d)) d <- d[in_analysis == TRUE]
d <- d[!is.na(HDI) & !is.na(overlap_pct) & !is.na(temporal_evenness)]
d[, overlap_prop := overlap_pct/100]
cat("species:", nrow(d), "\n")

## ---- classes recomputed for THIS index ---------------------------------
nz <- d[HDI > 0]; CUTS <- round(quantile(nz$HDI, c(1/3, 2/3)), 1)
cat(sprintf("\nHDI terciles among non-zero species: %.1f and %.1f\n", CUTS[1], CUTS[2]))
cat("  The earlier figure used <10 / 10-30 / >30, which came from the previous\n")
cat("  index. This one has a median of", round(median(nz$HDI),1),
    "and a maximum of", round(max(nz$HDI),1), "- under those\n  old breaks almost every species would fall in the top class.\n")
d[, cls := fcase(HDI == 0, "0 (unscorable)",
                 HDI <  CUTS[1], sprintf("< %.0f", CUTS[1]),
                 HDI <  CUTS[2], sprintf("%.0f-%.0f", CUTS[1], CUTS[2]),
                 default = sprintf("> %.0f", CUTS[2]))]
LEV <- c("0 (unscorable)", sprintf("< %.0f", CUTS[1]),
         sprintf("%.0f-%.0f", CUTS[1], CUTS[2]), sprintf("> %.0f", CUTS[2]))
d[, cls := factor(cls, levels=LEV)]
COL <- c("grey75","#C2544D","#E8A33A","#2C7FB8"); names(COL) <- LEV
print(d[, .(species=.N), by=cls][order(cls)])

## =========================================================================
## VERSION 1 - the layout as drawn
## =========================================================================
pdf(file.path(OUT,"fig_3d_as_drawn.pdf"), width=13, height=9)
par(mar=c(3,3,4,10), xpd=TRUE)
scatter3D(d$HDI, d$temporal_evenness, d$overlap_prop,
          colvar=NULL, col=COL[as.character(d$cls)], pch=16, cex=.45,
          theta=THETA, phi=PHI, bty="g", ticktype="detailed",
          xlab="HDI", ylab="Temporal evenness", zlab="Overlap (proportion)",
          main="Documentation in three dimensions",
          cex.lab=.85, cex.axis=.7)
legend("right", inset=c(-.16,0), legend=LEV, col=COL, pch=16, pt.cex=1.3,
       bty="n", title="HDI class", cex=.9)
mtext("Axes 2 and 3 are inputs to axis 1; the structure is partly arithmetic",
      side=1, line=1, cex=.75, col="grey40")
dev.off()
cat("\nwrote fig_3d_as_drawn.pdf\n")

## =========================================================================
## VERSION 2 - three inputs, HDI as colour
## =========================================================================
pdf(file.path(OUT,"fig_3d_inputs.pdf"), width=13, height=9)
par(mar=c(3,3,4,10), xpd=TRUE)
scatter3D(d$intensity, d$temporal_evenness, d$overlap_prop,
          colvar=d$HDI, col=ramp.col(c("#C2544D","#E8A33A","#8CC63E","#2C7FB8")),
          pch=16, cex=.45, theta=THETA, phi=PHI, bty="g", ticktype="detailed",
          xlab="Specimen intensity", ylab="Temporal evenness", zlab="Overlap (proportion)",
          main="How the index maps onto its inputs",
          clab="HDI", cex.lab=.85, cex.axis=.7)
mtext("Position and colour are independent: the axes are inputs, the colour is the score",
      side=1, line=1, cex=.75, col="grey40")
dev.off()
cat("wrote fig_3d_inputs.pdf\n")

## ---- a few viewing angles, so you can pick ------------------------------
pdf(file.path(OUT,"fig_3d_angles.pdf"), width=13, height=10)
par(mfrow=c(2,2), mar=c(2,2,3,2))
for (th in c(25, 45, 135, 215)) {
  scatter3D(d$intensity, d$temporal_evenness, d$overlap_prop,
            colvar=d$HDI, col=ramp.col(c("#C2544D","#E8A33A","#8CC63E","#2C7FB8")),
            pch=16, cex=.3, theta=th, phi=PHI, bty="g",
            xlab="Intensity", ylab="Temporal evenness", zlab="Overlap",
            main=paste("theta =", th), colkey=FALSE, cex.lab=.7, cex.axis=.6)
}
dev.off()
cat("wrote fig_3d_angles.pdf - pick an angle and set THETA\n")

## ---- species to annotate ------------------------------------------------
cat("\n=== candidates for photo annotation ===\n")
cat("highest HDI:\n")
print(d[order(-HDI)][1:5, .(species, HDI=round(HDI,1), n_sheets,
       temporal_evenness=round(temporal_evenness,2), overlap=round(overlap_pct,1))])
cat("\nhigh overlap, low intensity (the paired before/after cases):\n")
print(d[overlap_prop > .6 & intensity < median(intensity)][order(-overlap_prop)][1:5,
       .(species, HDI=round(HDI,1), n_sheets, overlap=round(overlap_pct,1))])
cat("\nmany specimens, low overlap (the published index's failure mode):\n")
print(d[order(-n_sheets)][overlap_prop < .15][1:5,
       .(species, HDI=round(HDI,1), n_sheets, overlap=round(overlap_pct,1))])
fwrite(d[, .(species, HDI, HDI_rank, cls, intensity, temporal_evenness,
             overlap_prop, n_sheets)], file.path(OUT,"threeD_points.csv"))

cat("\n---------------------------------------------------------------\n")
cat("WHICH TO USE\n---------------------------------------------------------------\n")
cat("  fig_3d_inputs is the defensible one. The axes are three inputs and the\n")
cat("  colour is the score, so nothing in the picture is forced by the formula.\n")
cat("  fig_3d_as_drawn reproduces your layout, but two of its axes feed the\n")
cat("  third - present it as an illustration, never as evidence.\n")
cat("\n  Both leave a wide right margin for the photographs.\n")
cat("\nwritten to", OUT, "\n"); print(Sys.time())
