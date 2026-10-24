# Herbarium Density Index (HDI) for the Alpine flora

Code to reproduce the Herbarium Density Index and all sensitivity analyses in:

> Ranjbaran Y, Canella M, Dubbini N, Singh G, Beltran-Sanz N, Petruzzellis F, Dal Grande F. A spatiotemporal prioritization framework for herbarium-based time-series research: the Herbarium Density Index applied to the Alpine flora. *New Phytologist* (in review).

Archived version: https://doi.org/10.5281/zenodo.XXXXXXX

## What is here

```
HDI_Ranjbaran_et_al.Rmd     main pipeline: grid, six inputs, HDI, PCA, results, sensitivity analyses
scripts/                    companion scripts for figures and descriptive analyses (read metrics_complete.csv)
LICENSE                     MIT licence for the code
CITATION.cff                how to cite this repository
```

## Input data

The pipeline reads four files, placed in the folder set by the environment variable `HDI_BASE` (default: the current folder):

| File | Content | Source |
|---|---|---|
| `herbaria_occurrences_filtered.csv` | herbarium specimens, clipped to the Alpine Convention perimeter | GBIF download, DOI: ⟦herbarium download DOI⟧ |
| `current_occurrences_enriched.csv` | contemporary occurrences, clipped to the same perimeter | GBIF download, DOI: ⟦modern download DOI⟧ |
| `species_index_full_info.csv` | Alpine checklist (`;`-separated) | ⟦source⟧ |
| `species_index_full_info_all_info.csv` | species attributes: taxonomy, WCVP status, establishment (`;`-separated) | ⟦source⟧ |

Raw occurrence data are not redistributed here. Download them from GBIF using the DOIs above.

## How to run

1. Install R (≥ ⟦version⟧) and the packages `data.table`, `sf`, `ggplot2`, `rmarkdown` and `knitr`.
2. Place the four input files in one folder and point to it:
   ```
   export HDI_BASE=/path/to/input/files
   ```
3. Render the main file:
   ```
   Rscript -e 'rmarkdown::render("HDI_Ranjbaran_et_al.Rmd")'
   ```
   Outputs are written to `HDI_output/`, including `metrics_complete.csv` (all six inputs, the four dimensions and HDI for every species).
4. Run the companion scripts in `scripts/` to produce the figures. They read `metrics_complete.csv` and recompute nothing.

The sensitivity analyses use 10,000 resamples each and are slow. On a cluster, set `SLURM_CPUS_PER_TASK` to use more cores. The random seed is fixed (`20260928`).

## Alpine Convention perimeter

Occurrences were clipped to the Alpine Convention Perimeter, edition v2_2018 (author: A. Piva, 2020; © Permanent Secretariat of the Alpine Convention, PSAC). The polygon is not redistributed here. Download it from the Alpine Convention Atlas: https://www.atlas.alpconv.org/layers/geonode:Alpine_Convention_Perimeter_2018_v2 (licence: CC BY 4.0).

## Credits

The first version of the report and the first implementation of the occurrence indices were written by Nevio Dubbini (Miningful Studio). The revised index construction and the sensitivity analyses were added by Yasaman Ranjbaran.

## Contact

Yasaman Ranjbaran (⟦email⟧) and Francesco Dal Grande (⟦email⟧), University of Padova.
