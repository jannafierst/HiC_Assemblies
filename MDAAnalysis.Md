First align sequences to a unified reference containing Scaffold 1, X_{A} and X_{B}

```

#!/bin/bash
#SBATCH --job-name=minimap2_illumina
#SBATCH --account acc_jfierst
#SBATCH --qos highmem1
#SBATCH --partition highmem1-sapphirerapids
#SBATCH --mem=384G
#SBATCH --cpus-per-task=32
#SBATCH --output=minimap2_illumina_%j.log
#SBATCH --time=42:00:00

module load miniconda3
module load samtools
source activate minimap2

# --- 1. USER-DEFINED VARIABLES ---
READS_DIR="./reads"
ASM="/home/data/jfierst/jfierst/worm-nomics/Bunonema_JU3390/JU3390_genome_files/JU3390_unified.fa"
THREADS=${SLURM_CPUS_PER_TASK:-16}

echo "Starting Illumina alignments to unified assembly (Unfiltered)..."

for R1 in "${READS_DIR}"/*_R1_001.fastq.gz; do
    SAMPLE=$(basename "${R1}" | cut -d'_' -f1)
    R2="${R1/_R1_001.fastq.gz/_R2_001.fastq.gz}"
    BAM_FILE="${SAMPLE}_mapped.bam"

    echo "-------------------------------------------------------"
    echo "Processing: $SAMPLE"
    echo "-------------------------------------------------------"

    # Removed -q filter. Storing all alignments.
    minimap2 -ax sr -t "${THREADS}" "$ASM" "${R1}" "${R2}" | \
    samtools view -@ 4 -b - | \
    samtools sort -@ "${THREADS}" -o "${BAM_FILE}"

    echo "Indexing $BAM_FILE..."
    samtools index "${BAM_FILE}"
done

echo "Done! All alignments completed."

./diagnostic_window_coverage.sh ./ ./diagnostic_windows.bed ./diagnostic_windows.csv

# Usage: ./targeted_bed_coverage.sh <bam_directory> <diagnostic_windows.bed> <output_csv>

'''

Then extract coverage

```
#!/bin/bash

# --- SLURM HEADER ---
#SBATCH --job-name=JU3390_unified_profiler
#SBATCH --account=acc_jfierst
#SBATCH --qos=highmem1
#SBATCH --partition=highmem1-sapphirerapids
#SBATCH --mem=34G
#SBATCH --cpus-per-task=8
#SBATCH --output=unified_profiler_%j.log
#SBATCH --mail-user=jfierst@fiu.edu
#SBATCH --mail-type=ALL
#SBATCH --time=04:00:00

# Exit immediately if a command fails
set -e 

# --- 1. ENVIRONMENT & PATH ARCHITECTURE ---
module load miniconda3
module load samtools

GENOME_DIR="/home/data/jfierst/jfierst/worm-nomics/Bunonema_JU3390/JU3390_genome_files"
ALIGN_DIR="/home/data/jfierst/jfierst/worm-nomics/Bunonema_JU3390/JU3390_MDA_Illumina/alignments"
FINAL_OUTPUT="${ALIGN_DIR}/JU3390_comprehensive_population_feature_matrix.txt"

# Initialize output with header
echo -e "Sample\tHaplotype_Context\tChromosome\tRegion_Name\tMean_Depth\tMedian_Depth\tVariance\tBreadth_Coverage_Pct\tMean_M
APQ\tMedian_MAPQ\tPct_Reads_MAPQ_0\tPct_Reads_Clipped_Or_Indel" > "${FINAL_OUTPUT}"

echo "========================================================================"
echo "    STARTING UNIFIED MULTI-DIMENSIONAL POPULATION COHORT PROFILING      "
echo "========================================================================"

# --- 2. FILE-FIRST LOOP: Detects your actual .bam files ---
for BAM_PATH in "${ALIGN_DIR}"/*_strict.bam; do
    
    [ -f "$BAM_PATH" ] || continue
    
    FILENAME=$(basename "${BAM_PATH}")
    
    # Extract Haplotype Context from filename
    # Parses 'JU3390_hap1' from 'F3-10_mapped_to_JU3390_hap1_strict.bam'
    
    # 1. Chop off everything before and including "_mapped_to_", then chop off "_strict.bam"
    HAP_DETECTED=$(echo "${FILENAME}" | sed 's/.*_mapped_to_//' | sed 's/_strict.bam//')
    
    # 2. Chop off "_mapped_to_" and everything after it to get the Sample ID
    SAMPLE=$(echo "${FILENAME}" | sed 's/_mapped_to_.*//')
    
    # BED file must match HAP_DETECTED
    BED_FILE="${GENOME_DIR}/${HAP_DETECTED}_domains.bed"
    
    if [ ! -f "${BED_FILE}" ]; then
        echo "[!] Warning: No BED file found for ${HAP_DETECTED} at ${BED_FILE}. Skipping."
        continue
    fi
    
    echo ">>> PROCESSING SAMPLE: ${SAMPLE} | HAPLOTYPE: ${HAP_DETECTED} <<<"

    # --- 3. CORE COMPUTATION ENGINE ---
    while read -r CHROM START END REGION; do
        [[ -z "$CHROM" || "$CHROM" =~ ^# ]] && continue
        
        REGION_STR="${CHROM}:$((START + 1))-${END}"

        # Depth stats
        stats_depth=$(samtools depth -a -Q 0 -G SUPPLEMENTARY -r "${REGION_STR}" "${BAM_PATH}" | awk '
            { sum += $3; sq_sum += $3 * $3; count++; freq[$3]++; if ($3 > 0) covered++ }
            END {
                if (count == 0) { printf "0.00\t0.00\t0.00\t0.00%%"; exit }
                avg = sum / count; var = (count > 1) ? (sq_sum - (sum^2 / count)) / (count - 1) : 0.0
                target = count / 2; cum = 0; med = 0
                for (d = 0; d <= 5000; d++) { cum += freq[d]; if (cum >= target) { med = d; break } }
                printf "%.2f\t%.2f\t%.2f\t%.2f%%", avg, med, (var < 0 ? 0 : var), (covered / count) * 100
            }')

        # Alignment stats (Filtering for MAPQ >= 1 to get valid scores)
        stats_align=$(samtools view -q 1 "${BAM_PATH}" "${REGION_STR}" | awk '
            { total_reads++; sum_mapq += $5; freq_mq[$5]++; if ($5 == 0) mapq_0++; if ($6 ~ /[SID]/) clipped++ }
            END {
                if (total_reads == 0) { printf "0.00\t0.00\t0.00%%\t0.00%%"; exit }
                target = total_reads / 2; cum = 0; med_mq = 0
                for (q = 0; q <= 60; q++) { cum += freq_mq[q]; if (cum >= target) { med_mq = q; break } }
                printf "%.2f\t%.2f\t%.2f%%\t%.2f%%", sum_mapq / total_reads, med_mq, (mapq_0 / total_reads) * 100, (clipped /
 total_reads) * 100
            }')

        echo -e "${SAMPLE}\t${HAP_DETECTED}\t${CHROM}\t${REGION}\t${stats_depth}\t${stats_align}" >> "${FINAL_OUTPUT}"
    done < "${BED_FILE}"
done

```

Next align sequences to each phased haplotype-resolved sequence individually

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



Extract coverage across domain seams

```

#!/bin/bash
# Usage: ./chimera_seam_hap1.sh <bam_directory> <output_csv>

ALIGN_DIR=$1
OUT_FILE=$2

# Define the exact physical seams on Hap1
SEAM_1="scaffold_1:1708900-1709100" # Domain A -> Domain B
SEAM_2="scaffold_1:3650900-3651100" # Domain B -> Domain C

echo "Sample,Reference,Seam,Total_Reads,Continuous,Continuous_Pct,Clipped,Clipped_Pct" > "${OUT_FILE}"
echo "Running Hap1 (Ancestral X) Seam Test..."

analyze_seam() {
    local region=$1
    local seam_name=$2
    local sample=$3
    local bam=$4

    samtools view -F 2308 "${bam}" "${region}" | awk -v sample="${sample}" -v seam="${seam_name}" '
    BEGIN { continuous=0; clipped=0; OFS=","; }
    {
        if ($6 ~ /[SH]/) { clipped++; } else { continuous++; }
    }
    END {
        total = continuous + clipped;
        if (total > 0) {
            pct_cont = (continuous / total) * 100;
            pct_clip = (clipped / total) * 100;
            printf "%s,Hap1,%s,%d,%d,%.2f,%d,%.2f\n", sample, seam, total, continuous, pct_cont, clipped, pct_clip;
        } else {
            printf "%s,Hap1,%s,0,0,0.00,0,0.00\n", sample, seam;
        }
    }' >> "${OUT_FILE}"
}

for BAM_FILE in "${ALIGN_DIR}"/*hap1*.bam; do
    if [[ -f "$BAM_FILE" ]]; then
        # Extracts "F3-10" from "F3-10_mapped_to_JU3390_hap1_strict.bam"
        SAMPLE_NAME=$(basename "${BAM_FILE}" | cut -d'_' -f1)
        echo "Processing ${SAMPLE_NAME}..."
        
        analyze_seam "${SEAM_1}" "Domain_A_to_B" "${SAMPLE_NAME}" "${BAM_FILE}"
        analyze_seam "${SEAM_2}" "Domain_B_to_C" "${SAMPLE_NAME}" "${BAM_FILE}"
    fi
done

echo "Done. Results written to ${OUT_FILE}"

```

Calculate concordance score in R and plot results

```
(library(tidyverse))
# We will use facet_grid to stack plots cleanly without needing extra packages

setwd("~/IDrive-Sync/Projects/Bunonema/JU3390_analyses/MDA_analyses/seam_analyses")

# 1. Load Data
matrix_df <- read_table("JU3390_comprehensive_unified_matrix.txt", show_col_types = FALSE) %>%
  filter(Sample != "F3-12") %>%
  mutate(Pct_Reads_MAPQ_0 = as.numeric(gsub("%", "", Pct_Reads_MAPQ_0))) %>%
  filter(Region_Name %in% c("Domain_A_s1_colinear", "Domain_C_XB", "Domain_C_XA")) %>%
  select(Sample, Region_Name, Pct_Reads_MAPQ_0, Mean_Depth, Variance) %>%
  pivot_wider(names_from = Region_Name, values_from = c(Pct_Reads_MAPQ_0, Mean_Depth, Variance))

seam_hap1 <- read_csv("hap1_seam_results.csv", show_col_types = FALSE) %>%
  filter(Sample != "F3-12", Seam == "Domain_A_to_B") %>%
  mutate(Seam_Backbone = ifelse(Total_Reads >= 2, Continuous_Pct, 0)) %>%
  select(Sample, Seam_Backbone)

seam_hap2 <- read_csv("hap2_seam_results.csv", show_col_types = FALSE) %>%
  filter(Sample != "F3-12", Seam %in% c("XB_Fusion", "XA_Inversion")) %>%
  mutate(Valid_Pct = ifelse(Total_Reads >= 2, Continuous_Pct, 0)) %>%
  select(Sample, Seam, Valid_Pct) %>%
  pivot_wider(names_from = Seam, values_from = Valid_Pct, names_prefix = "Seam_")

# 2. Compute Scores & Assign Specific Groups
scored_df <- matrix_df %>%
  inner_join(seam_hap1, by = "Sample") %>%
  inner_join(seam_hap2, by = "Sample") %>%
  mutate(
    `1_Backbone_SCS` = (100 - Pct_Reads_MAPQ_0_Domain_A_s1_colinear) * (Seam_Backbone / 100),
    `2_XA_SCS` = (100 - Pct_Reads_MAPQ_0_Domain_C_XA) * (Seam_XA_Inversion / 100),
    `3_XB_SCS` = (100 - Pct_Reads_MAPQ_0_Domain_C_XB) * (Seam_XB_Fusion / 100),
    `4_Stability_Variance` = Variance_Domain_C_XB / (Mean_Depth_Domain_C_XB + 0.01),
    
    # Apply detailed phenotypic/genotypic labels
    Group = case_when(
      Sample %in% c("F3-6", "F3-16") ~ "Putative Female",
      Sample == "M3-3" ~ "Putative Male",
      Sample == "M2-8" ~ "Male lacking XA and XB",
      Sample %in% c("M2-1", "M2-6") ~ "XB-carrying Male",
      str_starts(Sample, "F") ~ "Female",
      str_starts(Sample, "M") ~ "Male",
      TRUE ~ "Unknown"
    )
  ) %>%
  # Fix the legend order
  mutate(Group = factor(Group, levels = c(
    "Female", "Putative Female", 
    "XB-carrying Male", "Putative Male", 
    "Male", "Male lacking XA and XB"
  ))) %>%
  # Sort so Females are strictly on the left, Males on the right. 
  # Within those major blocks, order by XB presence, then Backbone strength.
  arrange(
    desc(str_starts(Sample, "F")), 
    desc(`3_XB_SCS`), 
    desc(`1_Backbone_SCS`)
  ) %>%
  mutate(Sample = factor(Sample, levels = unique(Sample)))

# 3. Reshape for faceted plotting with clean labels
plot_data <- scored_df %>%
  select(Sample, Group, `1_Backbone_SCS`, `2_XA_SCS`, `3_XB_SCS`, `4_Stability_Variance`) %>%
  pivot_longer(cols = starts_with("1_") | starts_with("2_") | starts_with("3_") | starts_with("4_"), 
               names_to = "Metric", values_to = "Value") %>%
  mutate(
    # Using plotmath syntax: atop() stacks text vertically to create multiline facet labels
    Metric_Label = case_when(
      str_detect(Metric, "1_Backbone") ~ "atop('Scaffold 1', 'Concordance')",
      str_detect(Metric, "2_XA") ~ "atop(X[A], 'Concordance')",
      str_detect(Metric, "3_XB") ~ "atop(X[B], 'Concordance')",
      str_detect(Metric, "4_Stability") ~ "atop(X[B]~'Dispersion', '(Variance/Mean)')"
    )
  ) %>%
  # Lock in the factor levels to guarantee exact top-to-bottom panel order
  mutate(Metric_Label = factor(Metric_Label, levels = c(
    "atop('Scaffold 1', 'Concordance')", 
    "atop(X[A], 'Concordance')", 
    "atop(X[B], 'Concordance')", 
    "atop(X[B]~'Dispersion', '(Variance/Mean)')"
  )))

# Define exact custom colors requested
custom_colors <- c(
  "Female" = "#E41A1C",                   # Standard Red
  "Putative Female" = "#FB9A99",          # Light Red
  "Male" = "#4DAF4A",                     # Standard Green
  "Putative Male" = "#B2DF8A",            # Light Green
  "Male lacking XA and XB" = "#005A32",   # Dark Green
  "XB-carrying Male" = "#377EB8"          # Standard Blue
)

# Use expression mapping to render subscripts properly in the legend
custom_breaks <- c(
  "Female", 
  "Putative Female", 
  "XB-carrying Male", 
  "Putative Male", 
  "Male", 
  "Male lacking XA and XB"
)

custom_labels <- expression(
  "Female", 
  "Putative Female", 
  X[B]*"-carrying Male", 
  "Putative Male", 
  "Male", 
  "Male lacking "*X[A]*" and "*X[B]
)

# 4. Generate Plot
p <- ggplot(plot_data, aes(x = Sample, y = Value, fill = Group, color = Group)) +
  # Bars for Concordance, Points/Lines for Variance
  geom_col(data = filter(plot_data, !str_detect(Metric, "Stability")), color = "black", linewidth = 0.2) +
  geom_point(data = filter(plot_data, str_detect(Metric, "Stability")), size = 5, alpha = 0.8) +
  geom_segment(data = filter(plot_data, str_detect(Metric, "Stability")), 
               aes(xend = Sample, yend = 0), linetype = "dotted", alpha = 0.5, linewidth = 0.8) +
  # labeller = label_parsed parses the math expressions (like atop and subscripts) for facets
  facet_grid(Metric_Label ~ ., scales = "free_y", switch = "y", labeller = label_parsed) +
  scale_fill_manual(values = custom_colors, breaks = custom_breaks, labels = custom_labels) +
  scale_color_manual(values = custom_colors, breaks = custom_breaks, labels = custom_labels) +
  theme_bw(base_size = 18) + # Increased base size for publication readiness
  theme(
    plot.title = element_text(size = 22, hjust = 0.5, margin = margin(b = 15)),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 14),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_text(face = "bold", size = 18, margin = margin(t = 15)),
    strip.placement = "outside",
    strip.text.y.left = element_text(face = "bold", size = 16, angle = 0),
    strip.background = element_rect(fill = "white", color = "white"),
    legend.position = "bottom",
    legend.title = element_text(face = "bold", size = 16),
    legend.text = element_text(size = 15),
    legend.key.size = unit(1.5, "lines"),
    # Remove all internal gridlines
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  labs(
    title = expression(bolditalic("Bunonema")~bold("JU3390 X Variants")),
    x = "",
    y = NULL,
    fill = "Cohort", color = "Cohort"
  )

p

output_pdf <- "Integrated_Alignment.pdf"
ggsave(output_pdf, plot = p, width = 12, height = 10, dpi = 300)
message(paste("Successfully generated:", output_pdf))

```
