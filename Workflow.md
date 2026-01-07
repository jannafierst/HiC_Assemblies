## Procedure and FIU HPC scripts for producing assembled genome sequences from Pacific Biosystems Hifi and Arima HiC libraries ###

***Libraries and best practices***

We assembled two sets of sequence libraries, one high depth Hifi (~40-400x coverage) and one high depth Arima HiC (~100-1000x coverage). We didn't know how large the genomes were going in for the non-Caenorhabditis species so ended up with very high sequencing coverage in some cases.

Each of these sequencing libraries is stored on 2 servers that are air-gapped (physically separated). Until the sequences are submitted to the NCBI sequence read archive (SRA) we need to be sure we have access to them in case of any catastrophic events. This can range from an accidental rm * to weather, flooding, etc. Our lab server on the FIU Roary cluster has a cloud back-up as well.

***

<details>
  <summary><b>Hifiasm assembly</b></summary>
