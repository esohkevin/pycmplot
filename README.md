# pycmplot

Multi-track **circular** and **linear** Manhattan plot generation for GWAS summary statistics.

```
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
|  PACKAGE FOR CIRCULAR AND LINEAR MANHATTAN PLOTTING  |
|                    Kevin Esoh, 2026                  |
|                    kesohku1@jh.edu                   |
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
```

---

## Installation

### From PyPI
```bash
pip install pycmplot
```


### From source
```bash
# From the repo root (use python >=3.8)
python -m venv ~/bin/pycmplot

source ~/bin/pycmplot/bin/activate

pip install --upgrade pip setuptools wheel

pip install -e .
```


### From source
```bash
# From the repo root (use python >=3.8)
python -m venv ~/bin/pycmplot

source ~/bin/pycmplot/bin/activate

pip install --upgrade pip setuptools wheel

pip install -e .
```


# Test the installation
```bash
pycmplot -h
```

### Dependencies

| Package | Purpose |
|---------|---------|
| pandas, numpy | Data loading & statistics |
| matplotlib | Plotting backend |
| pycirclize | Circular (Circos-style) tracks |
| natsort | Natural chromosome sorting |
| adjustText | Label collision avoidance |
| pyliftover | hg19 to hg38 coordinate conversion |
| Pillow | Image utilities |

---


## Command-line usage

### Linear Manhattan (default)

```bash
pycmplot \
  -s HbF.tsv.gz,MCV.tsv.gz,MCH.tsv.gz \
  -l HbF,MCV,MCH \
  --logp \
  -sig \
  -hl \
  -a GENE \
  -od ./results \
  -of png \
  --dpi 300
```

### Circular Manhattan

```bash
pycmplot \
  -s HbF.tsv.gz,MCV.tsv.gz \
  -l HbF,MCV \
  --mode cm \
  --logp \
  -sig \
  -plt "RBC Traits" \
  -od ./results
```

### Key options

| Flag | Description | Default |
|------|-------------|---------|
| `-s` | Comma-separated sumstats files | **required** |
| `-l` | Comma-separated track labels | **required** |
| `-m` | `lm` linear or `cm` circular | `lm` |
| `--logp` | Plot -log10(p) | off |
| `-sig` | Genome-wide significance line | off (auto 0.05/N) |
| `-sug` | Suggestive significance line | off |
| `-hl` | Highlight significant loci | off |
| `-a` | Annotate with `SNP` or `GENE` | `SNP` |
| `-tp` | Trim variants above this p-value for speed | off |
| `-st` | Sort tracks by `label` or `chrom_len` | input order |
| `-od` | Output directory | `.` |
| `-of` | Output format (`png`, `pdf`, `svg`, `jpg`) | `png` |

Run `pycmplot -h` for the full option list.

---

## Python API

```python
from pycmplot import multi_track_linear_manhattan
import pandas as pd

df1 = pd.read_csv("HbF.tsv.gz", sep="\t")
df2 = pd.read_csv("MCV.tsv.gz", sep="\t")

fig, axes = multi_track_linear_manhattan(
    tracks=[df1, df2],
    track_labels=["HbF", "MCV"],
    chr_col="CHR",
    pos_col="POS",
    value_col="P",
    logp=True,
    highlight=True,
    title="results/HbF_MCV.png",
    figsize=(15, 8),
)
```

---

_Under development_
