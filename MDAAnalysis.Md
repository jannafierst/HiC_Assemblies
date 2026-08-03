First align sequences

```

#!/bin/bash

# --- SLURM HEADER ---
#SBATCH --job-name=bowtie2_parallel_haps
#SBATCH --account=acc_jfierst
#SBATCH --qos=highmem1
#SBATCH --partition=highmem1-sapphirerapids
#SBATCH --mem=128G
#SBATCH --cpus-per-task=32
#SBATCH --output=bowtie2_parallel_%j.log
#SBATCH --mail-user=jfierst@fiu.edu
#SBATCH --mail-type=ALL
#SBATCH --time=12:00:00

# --- 1. MODULES AND ENVIRONMENT ---
module load samtools miniconda3

source activate minimap2

# --- 2. PATH ARCHITECTURE DECLARATIONS ---
GENOME_DIR="/home/data/jfierst/jfierst/worm-nomics/Bunonema_JU3390/JU3390_genome_files"
BASE_DIR="/home/data/jfierst/jfierst/worm-nomics/Bunonema_JU3390/JU3390_MDA_Illumina"
ALIGN_DIR="${BASE_DIR}/alignments"
READS_DIR="${BASE_DIR}/reads"

# Ensure the local alignment tracking space is established
mkdir -p "${ALIGN_DIR}"

# Fully reactivated array targeting your pristine, centralized references
HAPLOTYPES=(
    "JU3390_hap1"
    "JU3390_hap2_renamed"
    "JU3390_hap3_renamed"
)

# --- 3. DYNAMIC COMPUTE RESOURCE ALLOCATION ENGINE ---
TOTAL_CPUS=${SLURM_CPUS_PER_TASK:-32}

# Proportional distribution across the streaming thread pipe layout
BT2_THREADS=$(( TOTAL_CPUS * 75 / 100 ))    # 75% for index-matching and seed alignments
SORT_THREADS=$(( TOTAL_CPUS * 18 / 100 ))   # 18% for background coordinate sorting
VIEW_THREADS=$(( TOTAL_CPUS * 7 / 100 ))    # Remaining overhead for compressed BAM creation

# Fallback bounds to protect low-core interactive partitions
[ "${SORT_THREADS}" -eq 0 ] && SORT_THREADS=1
[ "${VIEW_THREADS}" -eq 0 ] && VIEW_THREADS=1

echo "========================================================================"
echo "    LAUNCHING PARALLEL MULTI-HAPLOTYPE ILLUMINA MDA ALIGNMENT PIPE     "
echo "========================================================================"
echo "Total Compute Threads Allocated  : ${TOTAL_CPUS}"
echo " -> Samtools Sorting Engine      : ${SORT_THREADS} cores"
echo " -> Samtools BAM Compression     : ${VIEW_THREADS} cores"
echo "------------------------------------------------------------------------"

# --- 4. OUTER LOOP: EVALUATE HAPLOTYPE REFERENCES ---
for HAP in "${HAPLOTYPES[@]}"; do
    GENOME="${GENOME_DIR}/${HAP}.fa"

    if [ ! -f "${GENOME}" ]; then
        echo "[!] Warning: Reference structural map missing at ${GENOME}. Skipping layout."
        continue
    fi

    # --- 5. INNER LOOP: INDEPENDENT READ REPLICATE STREAMING ---
    for R1 in "${READS_DIR}"/*_R1_001.fastq.gz; do
        
        # Safe skip if directory is empty or contains non-matching patterns
        [ -e "$R1" ] || continue
        
        SAMPLE=$(basename "${R1}" _R1_001.fastq.gz)
        R2="${READS_DIR}/${SAMPLE}_R2_001.fastq.gz"
        
        # Verify read pairs are intact before firing off the engine
        if [ ! -f "${R2}" ]; then
            echo "   [!] File Error: Orphan read detected. Missing partner track R2 for ${SAMPLE}."
            continue
        fi
        
        # Output paths aligned with downstream feature-matrix collection profiles
        SORTED_BAM="${ALIGN_DIR}/${SAMPLE}_mapped_to_${HAP}_strict.bam"
        TMP_PREFIX="${ALIGN_DIR}/${SAMPLE}_${HAP}_sort_tmp"

        echo "   -> Quantifying Variant: ${SAMPLE}"

	# Use this for all your Illumina samples
	minimap2 -ax sr -N 1 --secondary=no -t "${TOTAL_CPUS}" \
    	"${GENOME_DIR}/${HAP}.fa" \
    	"${R1}" "${R2}" | \
    	samtools sort -@ "${SORT_THREADS}" -o "${SORTED_BAM}" -
        
        samtools index -@ "${TOTAL_CPUS}" "${SORTED_BAM}"
        
    done
done

```
