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


### From GitHub
```bash
git clone https://github.com/esohkevin/pycmplot.git

cd pycmplot

pip install -e .

# or

pip install -e . --break-system-packages
```


### Use python virtual environment if local installation is not possible
```bash
python -m venv ~/bin/pycmplot

source ~/bin/pycmplot/bin/activate

pip install --upgrade pip setuptools wheel

# then follow any of the installation steps above
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
  --sum_stats HbF.tsv.gz,MCV.txt.gz,MCH.tsv.gz \
  --labels HbF,MCV,MCH \
  --logp \
  --signif_line \
  --highlight \
  --annotate GENE \
  --output_dir ./results \
  --output_format png \
  --dpi 300
```

### Circular Manhattan

```bash
pycmplot \
  --sum_stats HbF.tsv.gz,MCV.tsv.gz \
  --labels HbF,MCV \
  --mode cm \
  --logp \
  --signif_threshold \
  --plot_title "RBC Traits" \
  --output_dir ./results
```

### Key options

| Flag | Description | Default |
|------|-------------|---------|
| `-s, --sum_stats` | Comma-separated sumstats files | **required** |
| `-l, --labels` | Comma-separated track labels | **required** |
| `-m, --mode` | `lm` linear or `cm` circular | `lm` |
| `-qq, --qq_plot` | Also generate a QQ-plot | off (coming soon...) |
| `--logp` | Plot -log10(p) | off |
| `-sig, --signif_threshold` | Genome-wide significance threshold | off (auto 0.05/N) |
| `-sigl, --signif_line` | Value for genome-wide significance line if different from `-sig` | `-sig` |
| `-sug, --suggest_threshold` | Suggestive significance line | off |
| `-hl, --highlight` | Highlight significant loci | off |
| `-a, --annotate` | Annotate with `SNP` or `GENE` | `SNP` |
| `-tp, --trim_pval` | Trim variants above this p-value for speed | off |
| `-st, --sort_track` | Sort tracks by `label` or `chrom_len` | input order |
| `-od, --output_dir` | Output directory | `.` |
| `-of, --output_format` | Output format (`png`, `pdf`, `svg`, `jpg`) | `png` |

Run `pycmplot -h` for the full option list.

---

## Python API

```python
from pycmplot import plot_linear
import pandas as pd

df1 = pd.read_csv("HbF.tsv.gz", sep="\t")
df2 = pd.read_csv("MCV.tsv.gz", sep="\t")

plot_linear(
    tracks=[df1, df2],
    track_labels=["HbF", "MCV"],
    chr_col="CHR",
    pos_col="POS",
    p_col="P",
    logp=True,
    highlight=True,
    plot_title="results/HbF_MCV.png",
    figsize=(15, 8),
)
```

---

_Under development_
