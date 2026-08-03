First align the hifi sequences to each of the haplotypes

#!/bin/bash

# --- SLURM HEADER ---
#SBATCH --job-name=minimap2
#SBATCH --account acc_jfierst
#SBATCH --qos highmem1
#SBATCH --partition highmem1-sapphirerapids
#SBATCH --mem=384G
#SBATCH --cpus-per-task=32
#SBATCH --output=minimap2_%j.log
#SBATCH --mail-user=jfierst@fiu.edu
#SBATCH --mail-type=ALL
#SBATCH --time=24:00:00 # Recommended time limit for this process

# --- 0. MODULES AND ENVIRONMENT ---
# Load necessary modules (adjust paths/names if needed)
module load miniconda3
module load proxy # If required for network access
module load samtools
source activate minimap2

# --- 1. USER-DEFINED VARIABLES (CHANGE THESE) ---

DIR="/home/data/jfierst/frenchworms"
ASSEMBLY_DIR="/home/data/jfierst/jfierst/worm-nomics/Bunonema_JU3390/JU3390_genome_files"

# --- 1. SETTINGS & INPUTS ---
# Define an array of your assembly files
ASSEMBLIES=(
    "${ASSEMBLY_DIR}/JU3390_hap1.fa"
    "${ASSEMBLY_DIR}/JU3390_hap2_renamed.fa"
    "${ASSEMBLY_DIR}/JU3390_hap3_renamed.fa"
)

HIFI_READS="$DIR/JulyFrenchPacBio/reads/raw_libraries/JU3390.hifi.fastq"
THREADS=${SLURM_CPUS_PER_TASK:-16}

# --- 2. LOOP THROUGH ASSEMBLIES ---
for ASM in "${ASSEMBLIES[@]}"; do
    # Extract the filename without extension to use as a prefix
    ASM_NAME=$(basename "$ASM" .fa)
    BAM_FILE="${ASM_NAME}_hifi.bam"

    echo "-------------------------------------------------------"
    echo "Processing: $ASM_NAME"
    echo "Mapping HiFi reads to $ASM..."
    echo "-------------------------------------------------------"

    # Optimized for high-precision structural variation calling
    minimap2 -ax map-hifi -t "${THREADS}" --secondary=no -N 1 "$ASM" "${HIFI_READS}" | \
    samtools view -@ 4 -bS -q 30 - | \
    samtools sort -@ "${THREADS}" -o "${BAM_FILE}"

    # Map, convert to BAM, and sort in one pipe
    # Note: Using -@ for samtools to speed up compression/sorting
    #minimap2 -ax map-hifi -t "${THREADS}" "$ASM" "${HIFI_READS}" | \
    #    samtools view -@ 4 -bS - | \
    #    samtools sort -@ "${THREADS}" -o "${BAM_FILE}"

    echo "Indexing $BAM_FILE..."
    samtools index "${BAM_FILE}"

done

Then extract coverage 

#!/bin/bash
# comprehensive_windowed_depth.sh
# Usage: ./comprehensive_windowed_depth.sh <bam_directory> <window_size> <output_csv>

ALIGN_DIR=$1
WINDOW_SIZE=${2:-50000} 
OUT_FILE=$3

# Updated header to include Median
echo "Sample,Haplotype,Scaffold,Window_Start,Window_End,Median_Depth,Mean_Depth,Variance" > "${OUT_FILE}"
echo "Running universal whole-genome extraction in ${WINDOW_SIZE}bp windows..."

for BAM_FILE in "${ALIGN_DIR}"/*.bam; do
    if [[ -f "$BAM_FILE" ]]; then
        
        SAMPLE=$(basename "${BAM_FILE}" | cut -d'_' -f1)
        
        if [[ "${BAM_FILE}" == *"unified"* ]]; then
            HAP="unified"
        else
            HAP=$(echo "${BAM_FILE}" | grep -o -E 'hap[123]' || echo "unknown")
        fi
        
        if [[ -z "$HAP" ]]; then continue; fi

        echo "Processing ${SAMPLE} mapped to ${HAP}..."
        
        samtools idxstats "${BAM_FILE}" | grep -v "^\*" | while read SCAFFOLD LENGTH MAPPED UNMAPPED; do
            
            if [ "$LENGTH" -lt "$WINDOW_SIZE" ] || [ "$MAPPED" -eq 0 ]; then
                continue
            fi
            
            # Note: Changed to 'gawk' to ensure the sorting function works
            samtools depth -a -r "${SCAFFOLD}" "${BAM_FILE}" | gawk -v sample="${SAMPLE}" -v hap="${HAP}" -v chrom="${SCAFFOL
D}" -v ws="${WINDOW_SIZE}" '
            BEGIN { 
                OFS=","; 
                current_win = 0; 
                sum = 0; 
                sumsq = 0;
                count = 0; 
            }
            {
                pos = $2;
                depth = $3;
                win_idx = int((pos - 1) / ws);
                
                if (win_idx != current_win) {
                    if (count > 0) {
                        # Calculate Mean and Variance
                        mean = sum / count;
                        var = (sumsq / count) - (mean * mean);
                        if (var < 0) var = 0; 
                        
                        # Calculate Median using asort
                        asort(vals);
                        if (count % 2 == 1) {
                            median = vals[int((count+1)/2)];
                        } else {
                            median = (vals[count/2] + vals[(count/2) + 1]) / 2.0;
                        }
                        
                        print sample, hap, chrom, current_win * ws, (current_win + 1) * ws, median, mean, var;
                    }
                    current_win = win_idx;
                    sum = depth;
                    sumsq = depth * depth;
                    count = 1;
                    delete vals;
                    vals[count] = depth;
                } else {
                    sum += depth;
                    sumsq += depth * depth;
                    count++;
                    vals[count] = depth;
                }
            }
            END {
                if (count > 0) {
                    mean = sum / count;
                    var = (sumsq / count) - (mean * mean);
                    if (var < 0) var = 0;
                    
                    asort(vals);
                    if (count % 2 == 1) {
                        median = vals[int((count+1)/2)];
                    } else {
                        median = (vals[count/2] + vals[(count/2) + 1]) / 2.0;
                    }
                    
                    print sample, hap, chrom, current_win * ws, (current_win + 1) * ws, median, mean, var;
                }
            }' >> "${OUT_FILE}"
        done
    fi
done


echo "Done! All alignments completed."
