## Procedure and FIU HPC scripts for producing assembled genome sequences from Pacific Biosystems Hifi and Arima HiC libraries ###

***Libraries and best practices***

We assembled two sets of sequence libraries, one high depth Hifi (~40-400x coverage) and one high depth Arima HiC (~100-1000x coverage). We didn't know how large the genomes were going in so ended up with very high sequencing coverage in some cases.

Each of these sequencing libraries is stored on 2 servers that are air-gapped (physically separated). Until the sequences are submitted to the NCBI sequence read archive (SRA) we need to be sure we have access to them in case of any catastrophic events. This can range from an accidental rm * to weather, flooding, etc. Our lab server on the FIU Roary cluster has a cloud back-up as well but that is a little harder to access and has a delay.

Species we are studying:

JU760 	*Pelodera teres*

JU763 	*Poikilolaimus oxycercus*

JU3390 	*Bunonema sp.*

JU3391	*Bunonema sp.*

JU4118	*Caenorhabditis sp. 65*

JU4421	*Caenorhabditis sp. 70*

JU3778	*Caenorhabditis sp. 59*

JU4112	*Caenorhabditis sp. 62*

JU3779	*Caenorhabditis sp. 58*

JU4110	*Caenorhabditis sp. 61*

JU4113	*Caenorhabditis sp. 63*

JU4121	*Caenorhabditis sp. 66*

***

<details>
  <summary><b>Hifiasm assembly</b></summary>

Hifi and HiC libraries were assembled with Hifiasm using default parameters. I experimented intensely with Hifiasm parameters and different workflows for JU760 (*Pelodera teres*) as there was a lot of residual heterozygosity in the assembly. The highest quality assembled sequence (as measured by BUSCO completion and duplication statistics) resulted from Hifiasm default parameters and purge_dups post-assembly allelic removal so I chose to use that approach for each of these species.

https://hifiasm.readthedocs.io/en/latest/

https://github.com/dfguan/purge_dups

https://busco.ezlab.org/

#!/bin/bash
#SBATCH --account=acc_jfierst
#SBATCH --qos=highmem1
#SBATCH --partition=highmem1-sapphirerapids
#SBATCH --mem=384G
#SBATCH --cpus-per-task=32
#SBATCH --time=48:00:00
#SBATCH --output=logs/asm_%A_%a.out
#SBATCH --mail-user=jfierst@fiu.edu
#SBATCH --mail-type=ALL
#SBATCH --array=1-11%3

# 1. Environment & Paths
module load hifiasm
module load miniconda3
source /home/data/jfierst/miniconda3/etc/profile.d/conda.sh

BASE_DIR="/home/data/jfierst/frenchworms"
HIFI_DIR="$BASE_DIR/JulyFrenchPacBio/reads/raw_libraries"
ARIMA_DIR="$BASE_DIR/arima"
LINEAGE="/home/data/jfierst/nematoda_odb12"

# 2. Setup (Comment out after first run)
if [ "$SLURM_ARRAY_TASK_ID" -eq 1 ]; then
    mkdir -p "$BASE_DIR/logs"
    mkdir -p "$BASE_DIR/assemblies"
    ls "$ARIMA_DIR"/*_R1_001.fastq.gz | xargs -n 1 basename | sed 's/_R1_001\.fastq\.gz//' > "$BASE_DIR/species_list.txt"
fi
sleep 10 

# 3. Identify Task Species
ID=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$BASE_DIR/species_list.txt")
HIFI_FILE="$HIFI_DIR/${ID}.hifi.fastq"
R1_FILE="$ARIMA_DIR/${ID}_R1_001.fastq.gz"
R2_FILE="$ARIMA_DIR/${ID}_R2_001.fastq.gz"

mkdir -p "$BASE_DIR/assemblies/$ID"
cd "$BASE_DIR/assemblies/$ID" || exit

# STEP 1: hifiasm Phased (Restart-aware)
# If *.hic.lk.bin or *.hic.tlb.bin are deleted, this re-runs just the Hi-C part.
hifiasm -o "$ID" -t "$SLURM_CPUS_PER_TASK" \
    --h1 "$R1_FILE" \
    --h2 "$R2_FILE" \
    "$HIFI_FILE"

# STEP 2 & 3: Process both Haplotypes
for HAP in hap1 hap2; do
    GFA="${ID}.hic.${HAP}.p_ctg.gfa"
    FA="${ID}.${HAP}.p_ctg.fa"

    if [ -f "$GFA" ]; then
        echo "Processing $GFA..."
        # Extract sequences from GFA to FASTA
        awk '/^S/{print ">"$2;print $3}' "$GFA" > "$FA"
        
        # BUSCO Analysis
        echo "Running BUSCO for $ID $HAP..."
        source activate busco
        busco -c "$SLURM_CPUS_PER_TASK" -m genome -i "$FA" -o "${ID}_${HAP}_busco" --offline --lineage_dataset "$LINEAGE"
        conda deactivate

        # QUAST Analysis
        echo "Running QUAST for $ID $HAP..."
        source activate quast
        quast "$FA" -o "${ID}_${HAP}_quast"
        conda deactivate
    else
        echo "ERROR: Haplotype file $GFA not found for $ID."
    fi
done
