# Herbarium Density Index (HDI) for the Alpine flora

Code to reproduce the Herbarium Density Index and all sensitivity analyses in:

> Ranjbaran Y, Canella M, Dubbini N, Singh G, Beltran-Sanz N, Petruzzellis F, Dal Grande F. A spatiotemporal prioritization framework for herbarium-based time-series research: the Herbarium Density Index applied to the Alpine flora. *New Phytologist* (in review).

## What is here

```
HDI_Ranjbaran_et_al.Rmd     main pipeline: grid, six inputs, HDI, PCA, results, sensitivity analyses
scripts/                    companion scripts for figures and descriptive analyses
LICENSE                     MIT licence for the code
CITATION.cff                how to cite this repository
```

## Input data

The pipeline reads four files. They are not included in this repository because of their size.

| File | Content | Source |
|---|---|---|
| `herbaria_occurrences_filtered.csv` | herbarium specimens, clipped to the Alpine Convention perimeter | GBIF.org occurrence download, https://doi.org/10.15468/dl.hjwvd7 |
| `current_occurrences_enriched.csv` | contemporary occurrences, clipped to the same perimeter | GBIF.org occurrence download, https://doi.org/10.15468/dl.mvd7qu |
| `species_index_full_info.csv` | Alpine checklist (`;`-separated) | Canella et al. (2026), *Plants, People, Planet* 8: 680–692 |
| `species_index_full_info_all_info.csv` | species attributes: taxonomy, WCVP status, establishment (`;`-separated) | Canella et al. (2026), with status from the World Checklist of Vascular Plants (rWCVP) |

The GBIF download DOIs record the exact queries used. The clipped occurrence files and the two species files are available from the corresponding authors on request.

### Clipping to the Alpine perimeter

Both GBIF downloads were clipped manually in a GIS, keeping only records that fall inside the Alpine Convention Perimeter, edition v2_2018 (author: A. Piva, 2020; © Permanent Secretariat of the Alpine Convention, PSAC). The polygon is not redistributed here. Download it from the Alpine Convention Atlas: https://www.atlas.alpconv.org/layers/geonode:Alpine_Convention_Perimeter_2018_v2 (licence: CC BY 4.0).

No script exists for this step. The pipeline starts from the clipped files.

## How to run

1. Install R (tested with R 4.5.3) and the packages `data.table`, `sf`, `ggplot2`, `rmarkdown` and `knitr`. Some companion scripts also use `patchwork` and `plot3D`.
2. Place the four input files in one folder, and point to it:
   ```
   export HDI_BASE=/path/to/input/files
   ```
3. Render the main file:
   ```
   Rscript -e 'rmarkdown::render("HDI_Ranjbaran_et_al.Rmd")'
   ```
   Outputs are written to `HDI_output/`, including `metrics_complete.csv` (all six inputs, the four dimensions and HDI for every species).
4. Run the companion scripts from the same folder. They read `HDI_output/metrics_complete.csv` and write figures to `figures/`. Run `24_elevation.R` before `23_retention.R`, because the first writes `species_elevation.csv` and the second reads it. The table at the end of the Rmd lists what each script produces.

The sensitivity analyses use 10,000 resamples each and are slow. On a cluster, set `SLURM_CPUS_PER_TASK` to use more cores. The random seed is fixed (`20260928`).

## Credits

The first version of the report and the first implementation of the occurrence indices were written by Nevio Dubbini (University of Pisa). The revised index construction and the sensitivity analyses were added by Yasaman Ranjbaran.

## Contact

The corresponding authors of the article: Yasaman Ranjbaran and Francesco Dal Grande, University of Padova.
