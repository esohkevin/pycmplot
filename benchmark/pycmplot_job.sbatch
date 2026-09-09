#!/usr/bin/env bash
#SBATCH --account awonkam1
#SBATCH --job-name=bench_pycmplot
#SBATCH --partition=shared
#SBATCH --ntasks=1
##SBATCH -w sr47
#SBATCH --cpus-per-task=1
#SBATCH --mem=100G
#SBATCH --time=12:00:00
#SBATCH --output=/home/kesohku1/projects/gwas/scripts/pycmplotbenchmark/logs/bench_pycmplot_%j.out
#SBATCH --error=/home/kesohku1/projects/gwas/scripts/pycmplotbenchmark/logs/bench_pycmplot_%j.err

# -- Environment isolation: single core, no implicit threading
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export MPLBACKEND=Agg

# -- Activate Python environment
#source "$(conda info --base)/etc/profile.d/conda.sh"
module load anaconda
conda activate gwas_bench_py

echo "=============================="
echo "Job: bench_pycmplot"
echo "Node: $(hostname)"
echo "Date: $(date)"
echo "Python: $(python --version)"
echo "=============================="

cd /home/kesohku1/projects/gwas/scripts/pycmplotbenchmark

# for SIZE in 500K 1M 2M 5M 10M; do
#   INPUT="/scratch4/awonkam1/kesoh/gwas/pycmplotbenchmark/sumstats_${SIZE}.tsv"
#   if [[ ! -f "$INPUT" ]]; then
#     echo "WARNING: $INPUT not found <E2><80><94> skipping"
#     continue
#   fi
#   for PLOT_TYPE in manhattan circular qq; do
#     echo ""
#     echo ">>> pycmplot  size=${SIZE}  plot=${PLOT_TYPE}"
#     python \
#         bench_python.py \
#         --tool pycmplot \
#         --input "$INPUT" \
#         --size "${SIZE}" \
#         --plot-type "${PLOT_TYPE}" \
#         --replicates 5 \
#         --outdir "/home/kesohku1/projects/gwas/scripts/pycmplotbenchmark/results" \
#         --figdir "/home/kesohku1/projects/gwas/scripts/pycmplotbenchmark/figures"
#   done
# done

for SIZE in 5M 10M 25M 50M; do ##500K 1M 2M 5M 10M 25M 50M 100M; do
  INPUT="/scratch4/awonkam1/kesoh/gwas/pycmplotbenchmark/sumstats_${SIZE}.tsv"
  if [[ ! -f "$INPUT" ]]; then
    echo "WARNING: $INPUT not found <E2><80><94> skipping"
    continue
  fi
  for PLOT_TYPE in manhattan; do # circular; do
    echo ""
    echo ">>> pycmplot_multitrack  size=${SIZE}  plot=${PLOT_TYPE}"
      python bench_python.py \
          --tool pycmplot_multitrack \
          --manifest /scratch4/awonkam1/kesoh/gwas/pycmplotbenchmark/sumstats_${SIZE}_3traits.manifest \
          --size ${SIZE} \
          --plot-type multitrack_${PLOT_TYPE} \
          --replicates 5 \
          --outdir "/home/kesohku1/projects/gwas/scripts/pycmplotbenchmark/results/" \
          --figdir "/home/kesohku1/projects/gwas/scripts/pycmplotbenchmark/figures/"
  done
done

for SIZE in 5M 10M 25M 50M; do ##500K 1M 2M 5M 10M 25M 50M 100M; do
  INPUT="/scratch4/awonkam1/kesoh/gwas/pycmplotbenchmark/sumstats_${SIZE}.tsv"
  if [[ ! -f "$INPUT" ]]; then
    echo "WARNING: $INPUT not found <E2><80><94> skipping"
    continue
  fi
  for PLOT_TYPE in manhattan; do # circular; do
    echo ""
    echo ">>> pycmplot_multitrack  size=${SIZE}  plot=${PLOT_TYPE}"
      python bench_python.py \
          --tool pycmplot_multitrack \
          --manifest /scratch4/awonkam1/kesoh/gwas/pycmplotbenchmark/sumstats_${SIZE}_3traits.manifest \
          --size ${SIZE} \
          --plot-type multitrack_${PLOT_TYPE} \
          --replicates 5 \
          --outdir "/home/kesohku1/projects/gwas/scripts/pycmplotbenchmark/results/" \
          --figdir "/home/kesohku1/projects/gwas/scripts/pycmplotbenchmark/figures/" \
          --trim 0.01
  done
done

echo ""
echo "Done: $(date)"
