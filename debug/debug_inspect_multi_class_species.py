#!/usr/bin/env python3

import polars as pl
import json

# Read multi-class entries
multis = pl.read_csv("multi_class.csv", sep="\t")
parsed = pl.read_csv("data/bacdive_scrape_20220201.json.parsed.csv", sep="\t")

parsed_unique = parsed.unique()
multis_unique = multis.unique()

m = multis_unique.join(parsed_unique, on="Species", how="left")

# In [10]: m.sample(4)
# Out[10]: 
# shape: (4, 11)
# ┌───────────────────────────┬───────────────────────────┬────────────────────┬──────────────────────────┬─────┬──────────────────┬──────────────────────────┬──────────────────────────┬──────────────────────────┐
# │ Species                   ┆ gtdb_taxonomy             ┆ genus              ┆ sp                       ┆ ... ┆ Oxygen tolerance ┆ BacDive-ID               ┆ Oxygen tolerance_right   ┆ Oxygen tolerance         │
# │ ---                       ┆ ---                       ┆ ---                ┆ ---                      ┆     ┆ ---              ┆ ---                      ┆ ---                      ┆ raw_right                │
# │ str                       ┆ str                       ┆ str                ┆ str                      ┆     ┆ str              ┆ str                      ┆ str                      ┆ ---                      │
# │                           ┆                           ┆                    ┆                          ┆     ┆                  ┆                          ┆                          ┆ str                      │
# ╞═══════════════════════════╪═══════════════════════════╪════════════════════╪══════════════════════════╪═════╪══════════════════╪══════════════════════════╪══════════════════════════╪══════════════════════════╡
# │ Corynebacterium lowii     ┆ d__Bacteria;p__Actinobact ┆ g__Corynebacterium ┆ s__Corynebacterium lowii ┆ ... ┆ null             ┆ 140633                   ┆ anaerobe,microaerophile  ┆ obligate                 │
# │                           ┆ eriota;...                ┆                    ┆                          ┆     ┆                  ┆                          ┆                          ┆ anaerobe,microaerophile  │
# │ Globicatella              ┆ d__Bacteria;p__Firmicutes ┆ g__Globicatella    ┆ s__Globicatella          ┆ ... ┆ null             ┆ 157190,156360,151650,151 ┆ anaerobe,facultative     ┆ anaerobe,facultative     │
# │ sulfidifaciens            ┆ ;c__Bac...                ┆                    ┆ sulfidifaciens           ┆     ┆                  ┆ 649,1512...              ┆ anaerobe,mi...           ┆ anaerobe,mi...           │
# │ Caballeronia terrestris   ┆ d__Bacteria;p__Proteobact ┆ g__Caballeronia    ┆ s__Caballeronia          ┆ ... ┆ null             ┆ 133974                   ┆ aerobe,anaerobe          ┆ aerobe,anaerobe          │
# │                           ┆ eria;c_...                ┆                    ┆ terrestris               ┆     ┆                  ┆                          ┆                          ┆                          │
# │ Rhodovulum steppense      ┆ d__Bacteria;p__Proteobact ┆ g__Rhodovulum      ┆ s__Rhodovulum steppense  ┆ ... ┆ null             ┆ 13751                    ┆ aerobe,anaerobe          ┆ aerobe,anaerobe          │
# │                           ┆ eria;c_...                ┆                    ┆                          ┆     ┆                  ┆                          ┆                          ┆                          │
# └───────────────────────────┴───────────────────────────┴────────────────────┴──────────────────────────┴─────┴──────────────────┴──────────────────────────┴──────────────────────────┴──────────────────────────┘

import IPython; IPython.embed()