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
  <summary><b>Hifiasm assembly, file conversion and BUSCO and QUAST evaluation</b></summary>

Hifi and HiC libraries were assembled with Hifiasm using default parameters. I experimented intensely with Hifiasm parameters and different workflows for JU760 (*Pelodera teres*) as there was a lot of residual heterozygosity in the assembly. The highest quality assembled sequence (as measured by BUSCO completion and duplication statistics) resulted from Hifiasm default parameters and purge_dups post-assembly allelic removal so I chose to use that approach for each of these species.

https://hifiasm.readthedocs.io/en/latest/

https://github.com/dfguan/purge_dups

https://busco.ezlab.org/

```
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

# STEP 2 & 3: Process both Haplotypes and Diploid
for TYPE in p_ctg hap1 hap2; do
    if [ "$TYPE" == "p_ctg" ]; then
        GFA="${ID}.hic.p_ctg.gfa"
        FA="${ID}.p_ctg.fa"
    else
        GFA="${ID}.hic.${TYPE}.p_ctg.gfa"
        FA="${ID}.${TYPE}.p_ctg.fa"
    fi

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
```

</details>

<details>
  <summary><b>Decontamination</b></summary>

</details>

<details>
  <summary><b>Removing retained heterozygosity</b></summary>

I tried both purge_dups and purge_haplotigs. I found that purge_dups removed more sequence when compared with purge_haplotigs and retained slightly better BUSCO statistics and decided to use purge_dups for these genomes.

```
#!/bin/bash
#SBATCH --job-name=PurgeAll
#SBATCH --account=acc_jfierst
#SBATCH --qos=highmem1
#SBATCH --partition=highmem1-sapphirerapids
#SBATCH --mem=384G
#SBATCH --cpus-per-task=32
#SBATCH --output=logs/purge_triple_%A_%a.out
#SBATCH --mail-user=jfierst@fiu.edu
#SBATCH --mail-type=ALL
#SBATCH --array=1-11%3
#SBATCH --time=12:00:00

# --- 0. ENVIRONMENT ---
module load miniconda3
module load samtools
module load minimap2
source /home/data/jfierst/miniconda3/etc/profile.d/conda.sh
conda activate purge_dups

# --- 1. PATHS & SPECIES IDENTIFICATION ---
BASE_DIR="/home/data/jfierst/frenchworms"
HIFI_DIR="$BASE_DIR/JulyFrenchPacBio/reads/raw_libraries"
ID=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$BASE_DIR/species_list.txt")
HIFI_READS="$HIFI_DIR/${ID}.hifi.fastq"
ASM_DIR="$BASE_DIR/assemblies/$ID"

cd "$ASM_DIR" || exit

# --- 2. LOOP THROUGH ALL THREE TARGETS ---
# This will process p_ctg, hap1, and hap2 sequentially
for TYPE in p_ctg hap1 hap2; do
    TARGET_FA="${ID}.${TYPE}.fa"
    
    if [ ! -f "$TARGET_FA" ]; then
        echo "Skipping $TYPE for $ID: File not found."
        continue
    fi

    echo "################################################"
    echo "Starting Purge for $ID - Assembly Type: $TYPE"
    echo "################################################"

    # Create a unique subdir for this specific version
    WORK_DIR="purge_${TYPE}"
    mkdir -p "$WORK_DIR"
    cd "$WORK_DIR" || exit

    # STEP A: Mapping Reads (Necessary for depth)
    minimap2 -x map-hifi -t "${SLURM_CPUS_PER_TASK}" "../$TARGET_FA" "$HIFI_READS" > "${TYPE}.hifi.paf"

    # STEP B: Calculate Depth Stats
    pbcstat "${TYPE}.hifi.paf"
    calcuts PB.stat > cutoffs 2> calcults.log

    # STEP C: Self-Alignment (To find duplications)
    split_fa "../$TARGET_FA" > "${TYPE}.split.fa"
    minimap2 -xasm5 -DP -t "${SLURM_CPUS_PER_TASK}" "${TYPE}.split.fa" "${TYPE}.split.fa" | gzip -c - > "${TYPE}.split.self.paf.gz"

    # STEP D: Run Purge
    purge_dups -2 -T cutoffs -c PB.base.cov "${TYPE}.split.self.paf.gz" > dups.bed 2> purge_dups.log

    # STEP E: Extract
    get_seqs -e dups.bed "../$TARGET_FA"
    
    # Rename for easy identification later
    mv purged.fa "${ID}_${TYPE}_purged.fa"
    mv hap.fa "${ID}_${TYPE}_haplotigs.fa"

    echo "Completed $TYPE for $ID"
    
    # Step back up to species dir for the next type in the loop
    cd ..
done
```
https://github.com/dfguan/purge_dups

https://bitbucket.org/mroachawri/purge_haplotigs/src/master/

</details>

<details>
  <summary><b>Scaffolding</b></summary>  

The workflow here is a little wonky. I used juicer for alignment and processing then YaHS for scaffolding then juicer for assembly visualization and assessment. It is complex but seems to work well.

Set up juicer for the FIU Roary cluster

https://github.com/aidenlab/juicer

```
module load miniconda3

source activate hic_analysis
# Configure conda channels in order of priority
# This only needs to be done once to configure and install software
conda config --add channels defaults       # Standard packages
conda config --add channels bioconda       # Bioinformatics packages
conda config --add channels conda-forge    # Community-maintained packages
conda config --set channel_priority strict # Preve
module load proxy

# Install essential tools for Hi-C analysis
conda install wget
conda install git
conda install samtools
conda install sra-tools
conda install fastqc
conda install trim-galore
conda install bwa
conda install java-jdk
conda install gcc
conda install make
conda install gawk
conda install parallel

mkdir -p hic_tools
cd hic_tools/

# Clone the Juicer repository
git clone https://github.com/aidenlab/juicer.git
cd juicer

# Create necessary directory structure for Juicer 2.0
mkdir -p scripts/common
cp CPU/*.* scripts/common
cp CPU/common/* scripts/common

# Download Juicer Tools JAR file
# Note - use the current stable release here
wget https://github.com/aidenlab/Juicebox/releases/download/v2.17.00/juicer_tools_2.17.00.jar   

# Create symbolic link with the expected name
mv juicer_tools_2.17.00.jar scripts/common/juicer_tools.jar
```

Use the Arima site positions python script to generate cut sites across the genome

```
#!/bin/bash

srun -p highmem1-sapphirerapids --account=acc_jfierst --qos=highmem1 -c 4 --mem=64G --pty bash
module load miniconda3
source activate juicer
module load proxy

GENOME=YOUR_GENOME_HERE  

python generate_site_positions_Arima.py -i ${GENOME}.fa -e GATC GANTC TTAA CTNAG -o ${GENOME}_GATC_GANTC_TTAA_CTNAG.txt
```

This is the generate_site_positions_Arima.py script

```
# Generate cut site positions in genome from given restriction enzymes (Multiple enzymes and 'N's are supported!)
# Author: Xiang Zhou
# Affiliation: Arima Genomics Inc.

# Usage: generate_site_positions_Arima.py [-h] -i INPUT -e ENZYMES [ENZYMES ...] -o OUTPUT
# Arguments:
#   -h, --help                show this help message and exit
#   -i INPUT                  Input FASTA filename
#   -e ENZYMES [ENZYMES ...]  Enzyme sequences
#   -o OUTPUT                 Output filename

# python generate_site_positions_Arima.py -i XXXXX.fa -e GATC GANTC TTAA CTNAG -o XXXXX_GATC_GANTC_TTAA_CTNAG.txt

import sys
import re
import argparse
import os.path
from Bio import SeqIO

def main():
	parser = argparse.ArgumentParser(description = "Generate cut site positions in genome from given restriction enzymes (Multiple enzymes and 'N's are supported!)")

	parser.add_argument('-i', metavar='INPUT', type=str, dest='input', help='Input FASTA filename', required=True)
	parser.add_argument('-e', metavar='ENZYMES', type=str.upper, nargs='+', dest='enzymes', help='Enzyme sequences', required=True)
	parser.add_argument('-o', metavar='OUTPUT', type=str, dest='output', help='Output filename', required=True)

	args = parser.parse_args()

	with open(args.input, 'r') as FASTA, open(args.output, 'w') as OUTFILE:
		for record in SeqIO.parse(FASTA, "fasta"):
			length = len(record.seq)
			positions = find_re_sites(str(record.seq), args.enzymes)

			if positions:
				OUTFILE.write("{} {} {}\n".format(record.id, " ".join(map(str, positions)), length))
			else:
				OUTFILE.write("{} {}\n".format(record.id, length))

def find_re_sites(sequence, enzymes):
	#positions = []
	regex = "|".join(enzymes)
	regex = regex.replace("N", "[ACGT]")
	print("Finding cut sites matching pattern: {} ...".format(regex))
	iter = re.finditer(regex, sequence, re.IGNORECASE)
	#print([1+m.start() for m in iter])

	return [1+m.start() for m in iter]

if __name__ == '__main__':
	main()
```

Once you have generated the site positions you can generate the preliminary file structure and accessory files. Juicer requires each directory and file follows the specific naming conventions.

```
#!/bin/bash

GENOME=YOUR_GENOME_HERE
LIBRARY=YOUR_HIC_LIBRARY_HERE

bwa index ${GENOME}.fa
samtools faidx ${GENOME}.fa
cut -f1,2 ${GENOME}.fa.fai > ${GENOME}.chrom.sizes
mkdir qc_reports
cd qc_reports/
fastqc ${LIBRARY} -o . -t 8

trim_galore --paired --quality 20 --stringency 3 --length 20 --fastqc --cores 4 --output_dir ./trimmed_data/ ${LIBRARY}_R1.fastq.gz ${LIBRARY}_R2.fastq.gz

mkdir fastq

cp ./trimmed_data/* ./fastq/

gunzip ./fastq/*
```

Run the juicer pipeline

```
#!/bin/bash

#SBATCH --account acc_jfierst
#SBATCH --qos highmem1
#SBATCH --partition highmem1-sapphirerapids
#SBATCH --mem=64G
#SBATCH --cpus-per-task=16
#SBATCH --output=out_hifi.log
#SBATCH --mail-user=jfierst@fiu.edu  #insert your own email
#SBATCH --mail-type=ALL

module load miniconda3

module load proxy

source activate hic_analysis

# STEP 7: Run Juicer pipeline
#-----------------------------------------------
 
# Set key variables for the analysis
DIR="/home/data/jfierst/frenchworms"
GENOME="JU760_purged"
REFERENCE_DIR="$DIR/hic_analyses/JU760/juicer/"
RESTRICTION_SITE="Arima"
RESTRICTION_FILE="$DIR/hic_analyses/JU760/juicer/JU760_purged_GATC_GANTC_TTAA_CTNAG.txt"
CHROM_SIZES="$DIR/hic_analyses/JU760_purged/JU760_purged.chrom.sizes"
THREADS=16
 
# Run the complete Juicer pipeline
bash $DIR/hic_tools/juicer/scripts/common/juicer.sh -D $DIR/hic_tools/juicer -d $DIR/hic_analyses/JU760/juicer/ -g ${GENOME} -s ${RESTRICTION_SITE} -p ${CHROM_SIZES} -z ${REFERENCE_DIR}/JU760_pur
ged.fa -t ${THREADS} -y ${RESTRICTION_FILE} --assembly            # Restriction enzyme cut sites
 
```

<details>
  <summary><b>Assembly analysis and visualization</b></summary>
