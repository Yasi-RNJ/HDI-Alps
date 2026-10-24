###############################################################################
## 27_collecting_curve.R
## How collecting effort has run through time - the context panel the paper
## does not yet have.
##
##   A  specimens per year
##   B  distinct species collected per year
##   C  cumulative species first collected  (when was the flora discovered?)
##   D  specimens per species per year      (effort per taxon, not just total)
##
## Also the figure that justifies the 1740 window start visually.
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
WIN <- c(1740, 2019)

cat("=== 27_collecting_curve.R ===\n"); print(Sys.time())

h <- fread(HERB, select=c("species","year"))
h[, species := trimws(species)]
cat("specimens:", nrow(h), "   dated:", h[!is.na(year), .N],
    sprintf("(%.0f%%)\n", 100*mean(!is.na(h$year))))

## ---- what is outside the window, and why it is excluded -----------------
cat("\n=== dated specimens outside the analysis window ===\n")
cat("before", WIN[1], ":", h[!is.na(year) & year < WIN[1], .N], "\n")
cat("after ", WIN[2], ":", h[!is.na(year) & year > WIN[2], .N], "\n")
pre <- h[!is.na(year) & year < WIN[1], .(specimens=.N), by=.(decade=floor(year/10)*10)][order(decade)]
if (nrow(pre)) { cat("\nby decade before the window:\n"); print(pre) }

hd <- h[!is.na(year) & year >= WIN[1] & year <= WIN[2]]
cat("\nin window:", nrow(hd), "specimens,", uniqueN(hd$species), "species\n")

yr <- hd[, .(specimens=.N, species=uniqueN(species)), by=year][order(year)]
yr[, per_species := round(specimens/species, 2)]

cat("\n=== peaks and troughs ===\n")
cat("peak specimen year :", yr[which.max(specimens), year], "with",
    yr[, max(specimens)], "specimens\n")
cat("peak species year  :", yr[which.max(species), year], "with",
    yr[, max(species)], "species\n")
cat("median specimens per year:", median(yr$specimens), "\n")
print(hd[, .(specimens=.N, species=uniqueN(species)),
         by=.(era=cut(year, c(1739,1849,1899,1949,1999,2019),
              labels=c("1740-1849","1850-1899","1900-1949","1950-1999","2000-2019")))][order(era)])

## ---- cumulative first collection ---------------------------------------
first <- hd[, .(first_year=min(year)), by=species][order(first_year)]
cum <- first[, .(new_species=.N), by=first_year][order(first_year)]
cum[, cumulative := cumsum(new_species)]
tot <- cum[, max(cumulative)]
for (p in c(.5,.9,.95,.99)) {
  y <- cum[cumulative >= p*tot, min(first_year)]
  cat(sprintf("%.0f%% of species first collected by %d\n", 100*p, y))
}

th <- theme_minimal(base_size=10) +
  theme(plot.title=element_text(face="bold", colour="#1F3864", size=11),
        plot.subtitle=element_text(colour="grey40", size=9),
        panel.grid.minor=element_blank())

pA <- ggplot(yr, aes(year, specimens)) +
  geom_col(fill="#2C7FB8", width=.9) +
  geom_smooth(method="loess", span=.15, se=FALSE, colour="#C2544D", linewidth=.8) +
  scale_x_continuous(limits=WIN, breaks=seq(1750,2000,50)) +
  labs(x="Year", y="Specimens",
       title="A   Collecting effort through time",
       subtitle=sprintf("%s dated specimens, %d-%d",
                format(nrow(hd), big.mark=","), WIN[1], WIN[2])) + th

pB <- ggplot(yr, aes(year, species)) +
  geom_col(fill="#6BAE75", width=.9) +
  geom_smooth(method="loess", span=.15, se=FALSE, colour="#C2544D", linewidth=.8) +
  scale_x_continuous(limits=WIN, breaks=seq(1750,2000,50)) +
  labs(x="Year", y="Distinct species",
       title="B   How much of the flora was touched each year",
       subtitle=sprintf("peak %d species in %d", yr[,max(species)], yr[which.max(species), year])) + th

pC <- ggplot(cum, aes(first_year, cumulative)) +
  geom_area(fill="#8EC3E0", alpha=.75) +
  geom_line(colour="#1F3864", linewidth=.6) +
  geom_hline(yintercept=0.9*tot, linetype=2, colour="grey35") +
  annotate("text", x=WIN[1]+8, y=0.9*tot, label="90%", vjust=-.5, size=3, colour="grey35") +
  scale_x_continuous(limits=WIN, breaks=seq(1750,2000,50)) +
  labs(x="Year of first collection", y="Species recorded to date",
       title="C   When the flora was first collected",
       subtitle=sprintf("%s species have a dated first specimen", format(tot, big.mark=","))) + th

pD <- ggplot(yr, aes(year, per_species)) +
  geom_line(colour="#E8A33A", linewidth=.5) +
  geom_smooth(method="loess", span=.2, se=FALSE, colour="#C2544D", linewidth=.8) +
  scale_x_continuous(limits=WIN, breaks=seq(1750,2000,50)) +
  labs(x="Year", y="Specimens per species",
       title="D   Effort per taxon, not just in total",
       subtitle="Rising values mean repeat collection rather than wider coverage") + th

fig <- (pA | pB) / (pC | pD) +
  plot_annotation(title="Two and a half centuries of Alpine collecting",
    theme=theme(plot.title=element_text(face="bold", size=14, colour="#1F3864")))
ggsave(file.path(OUT,"fig_collecting_curve.pdf"), fig, width=13, height=9)
ggsave(file.path(OUT,"fig_collecting_curve.png"), fig, width=13, height=9, dpi=300)
for (i in seq_along(list(pA,pB,pC,pD)))
  ggsave(file.path(OUT, paste0("fig_collecting_", LETTERS[i], ".pdf")),
         list(pA,pB,pC,pD)[[i]], width=6.5, height=4.2)
fwrite(yr,  file.path(OUT,"specimens_per_year.csv"))
fwrite(cum, file.path(OUT,"species_accumulation.csv"))

cat("\n---------------------------------------------------------------\n")
cat("WHAT THESE ARE FOR\n---------------------------------------------------------------\n")
cat("  A and B give the reader the collecting history the index is built on.\n")
cat("  They also justify the window: the counts before", WIN[1], "are listed above\n")
cat("  and are negligible, which is why those decades were excluded rather\n")
cat("  than carried as empty bins for every species.\n\n")
cat("  C is the one worth a sentence. If most species were first collected\n")
cat("  early, the record is deep for the flora as a whole and the problem is\n")
cat("  one of CONTINUITY, not of discovery - which is exactly what the index\n")
cat("  measures.\n\n")
cat("  D separates effort from coverage. A rise in A with D flat means more\n")
cat("  ground covered; a rise in both means the same taxa collected again.\n")
cat("\nwritten to", OUT, "\n"); print(Sys.time())
