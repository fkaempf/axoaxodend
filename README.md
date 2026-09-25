# axoaxodend

Take a set of presynaptic neurons and a set of postsynaptic neurons in the male CNS
(`male-cns:v1.0` on neuprint). For every synapse between them, say whether it lands on the
postsynaptic neuron's **dendrite** (axo-dendritic) or its **axon** (axo-axonic). Output is a
table plus a Clio-NG (neuroglancer) scene.

Two versions:

| file | axon/dendrite split | use for |
|---|---|---|
| `axoaxo_navis.py` | `navis.split_axon_dendrite` (Schlegel), synapse flow centrality | **numbers** |
| `axoaxo.R` | single-edge flow-centrality cut (Schneider-Mizell et al. 2016) | quick R scenes |

Use the Python version for anything you report. The R version cuts at one skeleton edge. On
weakly polarised neurons that cut can land anywhere: on AVLP743m body 531383 it gives a
segregation index of 0.01 where navis gets 0.20. On AVLP749m → PAM01/PAM02 it calls 13% of
synapses axo-axonic; navis calls none.

## Setup

You need a neuprint token. Get it from https://neuprint.janelia.org → Account → Auth token.

### Python (uv)

```sh
uv sync                      # creates .venv from uv.lock, on Python 3.12 (.python-version)
cp .env.example .env         # paste your token into .env
uv run axoaxo_navis.py --pre MBON01 --post PAM01 --name mbon01_pam01
```

The script looks for the token in `NEUPRINT_APPLICATION_CREDENTIALS`, then in `.env`, then
as `NEUPRINT_TOKEN` in `~/.Renviron`.

You can also call it from Python:

```python
from axoaxo_navis import axoaxo
syn, summary = axoaxo(pre=["MBON01"], post=["PAM01"], name="mbon01_pam01")
```

### R

Packages: `neuprintr`, `nat`, `igraph`, `nabor`, `dplyr`, `jsonlite`. Put the token in
`~/.Renviron` as `NEUPRINT_TOKEN=...`.

```r
source("axoaxo.R")
axoaxo(pre = "AVLP749m", post = c("PAM02", "PAM01"), name = "test")
```

## Inputs

`pre` and `post` each take bodyids, type names, or a mix. Type names are resolved on neuprint.
A body needs at least 20 input and 20 output sites to be split (`min_sites`). Bodies below
that, or with no synapses from `pre`, are listed and skipped.

## Output

Written to `scenes_navis/<name>/` (Python) or `scenes/<name>/` (R):

| file | contents |
|---|---|
| `synapses.csv` | one row per pre→post synapse: bodies, types, compartment, domain, segregation index, xyz |
| `summary.csv` | per post body × pre type: counts, % axo-axonic, segregation index |
| `scene.json` | Clio-NG state: one segmentation layer per neuron type and one point layer per pre type × domain, each toggled on its own |
| `scene_url.txt` | the same state inline in a clio-ng URL, so it opens without hosting anything |
| `<name>.webloc` | double-click in Finder to open the scene; opened automatically after a run |

The navis version also reports `axo-linker`: synapses on the high-flow bridge between axon and
dendrite, or on the cell body fibre.

## Reading the numbers

- **Check the segregation index before trusting a split.** It runs from 0 (inputs and outputs
  fully mixed) to 1 (fully polarised). Below about 0.1 the axon/dendrite boundary is uncertain.
  Open the scene and check that the axo-axonic points sit where the axon should be.
- **Only the postsynaptic side is classified.** "Axo-axonic" means the synapse lands on the
  target's axon. The script does not check whether the presynaptic site is on the partner's
  axon.
- **Pool before you compare.** A body with 1–4 synapses swings from 0% to 100% on a single
  label, so compare pooled counts, not per-body percentages.

## Files

- `neuroglancer_scene.R`: `make_mcns_scene()`, the base Clio-NG scene (EM, segmentation, brain and VNC shells)
- `cache/`, `cache_navis/`: per-body skeletons and synapses, downloaded once (not committed)
