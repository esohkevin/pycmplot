# pycmplot

Multi-track **circular** and **linear** Manhattan plot generation for GWAS summary statistics.

```
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
|  PACKAGE FOR CIRCULAR AND LINEAR MANHATTAN PLOTTING  |
|                    Kevin Esoh, 2026                  |
|                    kesohku1@jh.edu                   |
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
```

This package will take any number of per SNP/variant summary statistics, be it GWAS, 
selection scans (e.g. iHS, EHH, FST), etc and generate Manhattan plots. If given a single
file, a single one-track Manhattan plot will be generated. Multiple files will result in 
the generation of a multi-track stacked Manhattan plot. 

In the process, the package will generate a **hits summary table** for variants with p-value 
(or whatever statistic for significance is used) below the user-specified significance threshold. 
This hits summary table will contain annotated gene names, in addition to other annotations, that
would then be used to annotate the plots.

Importantly, the package allows for conversion of hg19 genomic coordinates to hg38 coordinates.
This ensures that summary stats obtained using different imputation panels, for instance, can be
processed in the same run. That is, users can simply concatenate multiple summary stats files together, 
such as those for the same trait but analysed using different imputation panels. Users only need to 
add a new column specifying the genome build (hg19 or hg38) of the variants. Then the `--build_column` 
option of the package should be used to indicate the column and then the package will liftover all 
postions in hg19 to hg38 ensuring that hits table generation and plotting are done with one unified 
corrdinate system.

A key functionality of the package is its ability to auto-detect certain columns if ommited on the 
command-line or python API:
- Chromosome column: `-chr, --chrom_column` or ommited
- Basepair position column: `-pos, --pos_column` or ommited
- SNP or Marker ID column: `-snp, --snp_column` or ommited
- P-value (or whatever value) column: `-p, --pval_column` or ommited
- Build version column: `-b, --build_column` or ommited


Candidate names for each of the columns is shown below.

```python
# Resolve column names
chr_candidates = [chrom, 'CHR', 'CHROM', 'Chromosome', '#CHROM', '#CHR', 'Chrom', 'chrom', 'chr', 'chromosome', '#chr', '#chrom']
pos_candidates = [pos, 'BP', 'POS', 'bp', 'pos', 'Basepair']
snp_candidates = [snp, 'SNP', 'RSID', 'rsID', 'MarkerName', 'MarkerID', 'Predictor', 'Marker', 'SNPID', 'ID']
pvl_candidates = [pcol, 'P', 'P-value', 'Wald_P', 'pvalue', 'p_val', 'pval']
bld_candidates = [build, 'BUILD', 'Genome', 'Genome_Build', 'Genome-build']
```

> NB: Upper and lower cases of the candidates are also considered, making each candidate expanded 3 times.


Since GWAS summary stats files can be very large, to improve speed and memory efficiency, it is 
**highly recommended** to use `-tp, --trim_pval` with a value to exclude variants with p-value above a 
certain threshold, e.g. `0.01 (1e-2)` or `0.001 (1e-3)`.

A potential useful application is **comparative visualization** of results from multiple imputation panels, 
multiple populations, or multiple traits to observe shared genetic architecture.

Read more in the package documentation page: https://pycmplot.readthedocs.io/en/latest/

---

## Installation

### From PyPI
```bash
pip install pycmplot
```


### From GitHub
```bash
git clone https://github.com/esohkevin/pycmplot.git

# or with most recent updates from development branch
# git clone -b dev https://github.com/esohkevin/pycmplot.git

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
  --trim_pval 0.01 \
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
| `-b, --build` | Comma-separated genome builds of sumstats  | off |
| `-bc, --build_column` | Genome build column name (containing hg18/hg19/hg38) | off |
| `-m, --mode` | `lm` linear or `cm` circular | `lm` |
| `-qq, --qq_plot` | Also generate a QQ-plot | off |
| `-qq_thin, --qq_thin` | Thin p-values for faster QQ-plotting | off |
| `--logp` | Plot -log10(p) | off |
| `-sig, --signif_threshold` | Genome-wide significance threshold | off (auto 0.05/N) |
| `-sigl, --signif_line` | Value for genome-wide significance line if different from `-sig` | 5e-8 |
| `-sug, --suggest_threshold` | Threshold for suggestive signals | off |
| `-hl, --highlight` | Highlight significant loci | off |
| `-a, --annotate` | Annotate with `snp`, `gene`, or any column in `hits_table` | `snp` |
| `-tp, --trim_pval` | Trim variants above this p-value for speed | off |
| `-st, --sort_track` | Sort tracks by `label` or `chrom_len` | input order |
| `-od, --output_dir` | Output directory | `.` |
| `-of, --output_format` | Output format (`png`, `pdf`, `svg`, `jpg`) | `png` |

Run `pycmplot -h` for the full option list.

---

## Python API

A demonstration of how to use the python API is provided in this notebook: https://github.com/esohkevin/pycmplot/blob/main/pycmplot_python_api.ipynb


---

## Contributing

See how to contribute here https://github.com/esohkevin/pycmplot?tab=contributing-ov-file


## Contribution

- [Kevin Esoh](https://github.com/esohkevin)

