###############################################################################
## 16_results.R  --  figures for the headline results
##   1. the HDI distribution
##   2. the four dimensions, on a common scale
##   3. endemics vs non-endemics
##   4. taxonomic order as a predictor
###############################################################################
suppressPackageStartupMessages({library(data.table); library(ggplot2)})

MET <- "HDI_output/metrics_complete.csv"
OUT <- "figures"
dir.create(OUT, showWarnings=FALSE, recursive=TRUE)

cat("=== 16_results.R ===\n"); print(Sys.time())
d <- fread(MET)
if ("in_analysis" %in% names(d)) d <- d[in_analysis == TRUE]
cat("species:", nrow(d), "\n")
th <- theme_minimal(base_size=11) +
  theme(plot.title=element_text(face="bold", colour="#1F3864"),
        plot.subtitle=element_text(colour="grey40"))

sc <- d[!is.na(HDI)]
cat("\n=== HDI distribution ===\n"); print(round(summary(sc$HDI),2))
cat("species scoring exactly 0:", sc[HDI==0,.N], sprintf("(%.1f%%)\n", 100*sc[HDI==0,.N]/nrow(sc)))
p1 <- ggplot(sc, aes(HDI)) +
  geom_histogram(bins=60, fill="#2C7FB8", colour="white", linewidth=.2) +
  geom_vline(xintercept=median(sc$HDI), colour="#C2544D", linetype=2) +
  labs(x="HDI", y="Species", title="How the index is distributed",
       subtitle=sprintf("%s species; median %.1f (dashed); %d score exactly zero",
                format(nrow(sc), big.mark=","), median(sc$HDI), sc[HDI==0,.N])) + th
ggsave(file.path(OUT,"fig_hdi_distribution.pdf"), p1, width=8, height=4.8)
ggsave(file.path(OUT,"fig_hdi_distribution.png"), p1, width=8, height=4.8, dpi=300)

DIMS <- c("intensity","spatial","temporal","overlap")
if (all(DIMS %in% names(d))) {
  L <- melt(d[, c("species", DIMS), with=FALSE], id.vars="species",
            variable.name="dimension", value.name="value")[!is.na(value)]
  cat("\n=== the four dimensions (all should sit on 0-1 with similar medians) ===\n")
  print(L[, .(min=round(min(value),3), q1=round(quantile(value,.25),3),
              median=round(median(value),3), q3=round(quantile(value,.75),3),
              max=round(max(value),3)), by=dimension])
  p2 <- ggplot(L, aes(dimension, value, fill=dimension)) +
    geom_violin(colour="white", linewidth=.3, alpha=.85) +
    geom_boxplot(width=.12, outlier.shape=NA, fill="white", linewidth=.35) +
    scale_fill_brewer(palette="Set2", guide="none") + coord_cartesian(ylim=c(0,1)) +
    labs(x=NULL, y="Normalised value", title="The four dimensions are on a common scale",
         subtitle="Each divided by its own attainable maximum") + th
  ggsave(file.path(OUT,"fig_dimensions.pdf"), p2, width=8, height=4.8)
  ggsave(file.path(OUT,"fig_dimensions.png"), p2, width=8, height=4.8, dpi=300)
}

if ("wcvp_endemic" %in% names(d)) {
  e <- sc[!is.na(wcvp_endemic) & wcvp_endemic != "" &
          tolower(wcvp_endemic) != "missing data"]
  e[, endemic := fifelse(tolower(as.character(wcvp_endemic)) %in%
                         c("endemic","true","yes","1"), "Alpine endemic", "Not endemic")]
  cat("\n=== endemic vs non-endemic ===\n")
  print(e[, .(species=.N, median_HDI=round(median(HDI),2), zero=sum(HDI==0)), by=endemic])
  if (uniqueN(e$endemic) == 2) {
    w <- wilcox.test(HDI ~ endemic, data=e)
    cat("\nWilcoxon rank-sum: W =", format(w$statistic, big.mark=","),
        " p =", format.pval(w$p.value, digits=3), "\n")
    cat("(a rank test: HDI is not normal and has a spike at zero)\n")
    p3 <- ggplot(e, aes(endemic, HDI, fill=endemic)) +
      geom_violin(colour="white", linewidth=.3, alpha=.85) +
      geom_boxplot(width=.13, outlier.shape=NA, fill="white", linewidth=.35) +
      scale_fill_manual(values=c("Alpine endemic"="#E3752C","Not endemic"="#2C7FB8"), guide="none") +
      labs(x=NULL, y="HDI", title="Alpine endemics are less well documented",
           subtitle=sprintf("Medians %.1f vs %.1f; Wilcoxon p = %s",
                     e[endemic=="Alpine endemic", median(HDI)],
                     e[endemic=="Not endemic", median(HDI)],
                     format.pval(w$p.value, digits=3))) + th
    ggsave(file.path(OUT,"fig_endemic.pdf"), p3, width=6.5, height=5)
    ggsave(file.path(OUT,"fig_endemic.png"), p3, width=6.5, height=5, dpi=300)
  }
  cat("\nspecies with no usable endemism value:",
      sc[is.na(wcvp_endemic) | wcvp_endemic=="" | tolower(wcvp_endemic)=="missing data", .N],
      " <- resolve before the Results are final\n")
}

if ("order" %in% names(d)) {
  o <- sc[!is.na(order) & order != ""]
  keep <- o[, .N, by=order][N >= 20, order]; o <- o[order %in% keep]
  cat("\n=== taxonomic order (orders with >= 20 species) ===\n")
  cat("orders kept:", length(keep), "  species:", nrow(o), "\n")
  k <- kruskal.test(HDI ~ factor(order), data=o)
  cat("Kruskal-Wallis: chi-sq", round(k$statistic,1), " df", k$parameter,
      " p", format.pval(k$p.value, digits=3), "\n")
  ss <- summary(aov(HDI ~ factor(order), data=o))
  r2 <- ss[[1]][1,"Sum Sq"] / sum(ss[[1]][,"Sum Sq"])
  cat("variance in HDI explained by order:", sprintf("%.1f%%", 100*r2),
      "\n  -> order is a weak predictor; species-level screening stays necessary\n")
  lev <- o[, .(med=median(HDI)), by=order][order(med), order]
  o[, order := factor(order, levels=lev)]
  p4 <- ggplot(o, aes(order, HDI)) +
    geom_boxplot(fill="#8EC3E0", outlier.size=.4, linewidth=.35) + coord_flip() +
    labs(x=NULL, y="HDI", title="Taxonomic order explains little of the variation",
         subtitle=sprintf("%d orders with 20+ species; order explains %.1f%% of the variance",
                   length(keep), 100*r2)) + th
  hh <- max(5, .22*length(keep) + 2)
  ggsave(file.path(OUT,"fig_order.pdf"), p4, width=8, height=hh)
  ggsave(file.path(OUT,"fig_order.png"), p4, width=8, height=hh, dpi=300)
}
cat("\nwritten to", OUT, "\n"); print(Sys.time())
