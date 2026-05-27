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

# Key features
## Column auto-detection
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

## Density-aware sub-sampling
Another key feature is density-aware sub-sampling for Manhattan-style scatter plots.
This was inspired by ``gwaslab``'s default behaviour (https://cloufield.github.io/gwaslab/). 

Every variant whose "interestingness" signal is at or above ``keep_threshold`` is preserved (so peaks, suggestive hits, genome-wide-significant hits, and extreme 
selection-scan values are kept verbatim). It uniformly sub-samples the dense bulk 
below the threshold down to at most ``max_below`` rows in total.  For a 10 M-variant 
scan with the defaults below, this typically cuts the plotted point count from 10 M 
to ~200 K + a few hundred peaks — visually indistinguishable above the suggestive 
band, but two orders of magnitude faster to render.

## Trim insignificant variants for faster plotting
An optional parameter `-tp, --trim_pval` is provided to increase speed even further. 
Set with a value to exclude variants with p-value above a certain threshold, 
e.g. `0.01 (1e-2)` or `0.001 (1e-3)`. Performed on top of the default auto-thin 
feature above, it siginificant increases speed and reduces peak memory usage. 
See benchmark figure (manuscript in preparation).

## Genome build conversion (liftover)
Conversion of a both hg18 and hg19 positions to their hg38 equivalent is included through
`pyliftover.LiftOver`.

This means you can concatenate multiple summary stats into one file and include a `BUILD` 
column to specify the genome build of each position ('hg18', 'hg19', or 'hg38') and all 
'hg18' and 'hg19' positions will be converted to 'hg38' so that all positions are plotted 
using one coordinate system. If only 'hg18' or 'hg19' positions are present, no liftover 
be necessary. Hence, liftover is only performed in cases of mixed genome builds.

## Nearest-gene annotation for GWAS lead SNPs
The package bundles GFF3 files in hg19 and hg38 coordinates processed to reduce size 
for gene annotation. Also included are UCSC chain files for coordinate conversion (liftover).
  - ``chain_hg19_hg38`` -- UCSC LiftOver chain file for hg19 to hg38
    conversion. Resolved from ``PYCMPLOT_CHAIN_HG19_HG38`` or the bundled
    ``hg19ToHg38.over.chain.gz``.
  - ``chain_hg18_hg38`` -- UCSC LiftOver chain file for hg18 to hg38
    conversion. Resolved from ``PYCMPLOT_CHAIN_HG18_HG38`` or the bundled
    ``hg18ToHg38.over.chain.gz``. Only required when any input summary
    statistics file carries a ``hg18`` build label.
  - ``geneinfo_hg38`` -- Ensembl gene-info TSV for GRCh38, used for
    nearest-gene annotation. Resolved from ``PYCMPLOT_GENEINFO_HG38`` or
    the bundled ``Homo_sapiens.GRCh38.geneinfo.tsv.gz``.
  - ``geneinfo_hg19`` -- Ensembl gene-info TSV for GRCh37, used when
    input data carry a hg19 build label. Resolved from
    ``PYCMPLOT_GENEINFO_HG19`` or the bundled
    ``Homo_sapiens.GRCh37.geneinfo.tsv.gz``.


# Application
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


## Contributors

- [Kevin Esoh](https://github.com/esohkevin)

