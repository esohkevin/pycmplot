#!/usr/bin/env bash
#SBATCH --account awonkam1
#SBATCH --job-name=bench_CMplot
#SBATCH --partition=shared
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=100G
##SBATCH -w sr47
#SBATCH --time=12:00:00
#SBATCH --output=/home/kesohku1/projects/gwas/scripts/pycmplotbenchmark/logs/bench_CMplot_%j.out
#SBATCH --error=/home/kesohku1/projects/gwas/scripts/pycmplotbenchmark/logs/bench_CMplot_%j.err

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export R_DATATABLE_NUM_THREADS=1

#source "$(conda info --base)/etc/profile.d/conda.sh"
module load anaconda
conda activate gwas_bench_r

echo "=============================="
echo "Job: bench_CMplot"
echo "Node: $(hostname)"
echo "Date: $(date)"
echo "R: $(R --version | head -1)"
echo "=============================="

cd /home/kesohku1/projects/gwas/scripts/pycmplotbenchmark

OUTDIR="results"
FIGDIR="figures"
REPLICATES=5

# for SIZE in 10M; do #500K 1M 2M 5M 10M; do
#   INPUT="/scratch4/awonkam1/kesoh/gwas/pycmplotbenchmark//sumstats_${SIZE}.tsv"
#   if [[ ! -f "$INPUT" ]]; then
#     echo "WARNING: $INPUT not found <E2><80><94> skipping"
#     continue
#   fi
#   for PLOT_TYPE in manhattan circular qq; do
#     echo ""
#     echo ">>> CMplot  size=${SIZE}  plot=${PLOT_TYPE}"
# 
#     bash bench_r_timed.sh \
#       --tool CMplot --input "$INPUT" --size "$SIZE" \
#       --plot-type "$PLOT_TYPE" --replicates "$REPLICATES" \
#       --outdir "$OUTDIR" --figdir "$FIGDIR"
# 
#     #Rscript bench_r.r \
#     #  --tool CMplot \
#     #  --input "$INPUT" \
#     #  --size "${SIZE}" \
#     #  --plot-type "${PLOT_TYPE}" \
#     #  --replicates 5 \
#     #  --outdir "/home/kesohku1/projects/gwas/scripts/pycmplotbenchmark/results" \
#     #  --figdir "/home/kesohku1/projects/gwas/scripts/pycmplotbenchmark/figures"
#   done
# done

for SIZE in 1M 2M 5M 10M 25M 50M; do ##500K 1M 2M 5M 10M 25M 50M 100M; do
  INPUT="/scratch4/awonkam1/kesoh/gwas/pycmplotbenchmark/sumstats_${SIZE}.tsv"
  if [[ ! -f "$INPUT" ]]; then
    echo "WARNING: $INPUT not found <E2><80><94> skipping"
    continue
  fi
  for PLOT_TYPE in manhattan circular; do
    echo ""
    echo ">>> CMplot_multitrack  size=${SIZE}  plot=${PLOT_TYPE}"

    bash bench_r_timed.sh \
      --tool CMplot_multitrack --size "$SIZE" \
      --manifest /scratch4/awonkam1/kesoh/gwas/pycmplotbenchmark/sumstats_${SIZE}_3traits.manifest \
      --plot-type "$PLOT_TYPE" --replicates "$REPLICATES" \
      --outdir "$OUTDIR" --figdir "$FIGDIR"

     # Rscript bench_r.r \
     #     --tool CMplot_multitrack \
     #     --manifest /scratch4/awonkam1/kesoh/gwas/pycmplotbenchmark/sumstats_${SIZE}_3traits.manifest \
     #     --size ${SIZE} \
     #     --plot-type multitrack_${PLOT_TYPE} \
     #     --replicates 5 \
     #     --outdir "/home/kesohku1/projects/gwas/scripts/pycmplotbenchmark/results/" \
     #     --figdir "/home/kesohku1/projects/gwas/scripts/pycmplotbenchmark/figures/"
  done
done

echo ""
echo "Done: $(date)"
