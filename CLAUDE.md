# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is an R-based research project computing global food-related **energy and labor time footprints** by linking two multi-regional input-output (MRIO) databases:

- **FABIO** (Food and Agriculture Biomass Input-Output): 187 countries × 123 food commodities, mass-based
- **EXIOBASE3** (pxp, 2020): 49 regions × 200 sectors, monetary-based, with energy (TJ) and labor (M.hr) satellite accounts
- **Global Human Day (GHD/Fajzel)** dataset: non-economic food time (preparation, processing, growth/collection) by gender

The key research question is: how much energy and human time (by gender, domestic vs. trade-embedded) does the global food system require, normalized per capita and per calorie/protein unit?

## Script Execution Order

Scripts must be sourced in order — each builds on objects left in the R environment by the previous:

1. **`0.mrio_prep.R`** — Load raw data: FABIO matrices (L, x, Y), EXIOBASE satellites (energy, labor), region/product mappings, GHD gender-split food time data. Sets `year = 2020`, `type = "pxp"`.
2. **`1.mrio_convert.R`** — Convert FABIO mass vectors to EXIO classification; compute direct energy/labor satellites in FABIO units; compute calorie/protein nutrient flows. Saves `data/FABIO_exio_satellites_food_2020.rds`.
3. **`1.1.mrio_convert_indirect.R`** — Derive indirect (non-food sector) energy/labor satellites via Leontief inverse. Saves `data/FABIO_exio_satellites_nonfood_2020.rds`.
4. **`2.analyze_result.R`** — Load saved satellites, compute country-level footprint matrices, aggregate domestic/export/import summaries, plot results, save to `results/` and `output/`.

**`99.utils.R`** is sourced inside `1.mrio_convert.R` and provides shared utility functions (reordering, intensity conversion, plotting helpers).

**`exiobase.R`** is a standalone helper (`readExio`, `exioloop`) for constructing arbitrary EXIOBASE3 indicator matrices; not part of the main pipeline.

## External Data Paths (hardcoded)

Data paths are absolute and machine-specific. Key locations:
- `H:/MyDocuments/Data/EXIOBASE3/IOT_2020_pxp/` — EXIOBASE3 matrices (A, Y, x, unit.txt, energy/F.txt, employment/F.txt)
- `H:/MyDocuments/Data/FABIO/input/` — FABIO matrices (L, x, Y, regions.csv, items.csv, io_labels.csv, etc.)
- `H:/MyDocuments/Data/GlobalHumanDay/` — GHD output CSVs (gender split, restaurants)
- `data/EXIO_L_2020_pxp.rds` — Pre-computed EXIO Leontief inverse (large; computed once, then loaded)

## Key Dimensional Facts

| Object | Dimensions | Description |
|---|---|---|
| FABIO x | 23001 (=187×123) | Mass output vector |
| FABIO L | 23001×23001 | Leontief inverse |
| EXIO x | 9800 (=49×200) | Monetary output vector |
| EXIO L | 9800×9800 | Leontief inverse |
| `prod_map` | 123×200 | FABIO-to-EXIO product concordance |
| `p_fabio_exio` | 37400×37400 | Block-diagonal product mapping (187 countries) |
| `FABIO_x_in_EXIO` | 23001×37400 | FABIO mass mapped to EXIO sectors |
| `total_intensity_fabio` | 8575×37400 | Non-food sector intensities in FABIO space |

## Architecture: The Bridging Problem

The core challenge is that FABIO is mass-based (tonnes of food commodities) and EXIOBASE is monetary (USD). The pipeline bridges these via:

1. **`convert_mass_vecs()`** — Uses EXIO monetary shares within each FABIO-to-EXIO concordance row to split FABIO mass into EXIO sectors. Output: `exio_mass_x` (mass in EXIO classification) and `FABIO_x_in_EXIO` (the full mapping matrix).

2. **`convert_intensities()`** — Converts EXIO energy/labor intensities (per USD) to per-tonne intensities, then reorders from EXIO 49-region to FABIO 187-country space via `reorder_countries_to_FABIO()`.

3. **`reorder_countries_to_FABIO()` / `reorder_countries_to_EXIO()`** — Reshape vectors/matrices between the 49 EXIO regions and 187 FABIO countries using `FABIO_reg$EXIOBASE_code` as index. RoW regions paste single EXIO values to all FABIO countries in that aggregate.

## Footprint Types

The final footprint has two components:
- **Food-sector** (`l_int_d`): direct EXIO energy/labor mapped through FABIO supply chain
- **Non-food-sector** (`l_int_i`): indirect energy/labor from non-food EXIO sectors (packaging, transport, etc.) flowing through food production

Plus **non-economic food time** from GHD (household preparation, processing, growth/collection) added at the analysis stage.

All results are decomposed into **domestic**, **export** (origin country's footprint for goods it produces for others), and **import** (footprint embedded in imports consumed domestically).

## Key Conventions

### Naming Patterns
- `int_d` = direct intensity (energy or labor per kg at production point)
- `int_i` = indirect intensity (embodied in non-food supply chains)
- `sat_*` = satellite accounts (e.g., `sat_energy`, `sat_labor`)
- `x` = production vector/matrix, `Y` = final demand, `L` = Leontief inverse
- FABIO = food-specific IO system; EXIO = full economy IO system

### Matrix Operations
- Use `Matrix` package sparse types (`dgCMatrix`) for large IO tables
- Row/column scaling uses `Diagonal()` and `sweep()` rather than dense multiplication
- NaN and Inf produced by division (0/0 or x/0) are replaced with 0
- Kronecker products expand multi-region matrices

### Physical Units
- Energy: TJ (terajoules), converted to J/kg for intensities
- Labor: Million hours (M.hr) by gender; units are hr/kg for intensities
- Food mass: tonnes; calories: kcal; protein: grams

### Dependencies
Core packages: `mrio`, `tidyverse`, `Matrix`, `data.table`, `readxl`, `gt`, `countrypops`
