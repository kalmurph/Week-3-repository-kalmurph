# Instructions

## Edit the Script:

1. Open the script in a text editor (`vi`):

    ```bash
    vi sra-download.slurm
    ```

2. Replace `your-email@example.com` with your actual email address.

3. Select an SRA accession number from the **1000 Genomes Project**:
   - Visit the [1000 Genomes Project Bioproject Page](https://www.ncbi.nlm.nih.gov/bioproject/59771).
   - Scroll down to explore different **human populations** that have been sequenced.
   - Identify an **SRR accession number** that interests you.

4. Replace `SRX123456` in the script with your chosen **SRR accession number**.

5. Replace `your-username` in the `STORAGE_DIR` variable with your actual username.

6. Save and exit `vi`:
   - Press `Esc`, then type `:wq` and press **Enter**.

## Submit the Job:

1. Navigate to the directory where the script is located (if not already there):

    ```bash
    cd /path/to/your/script
    ```

2. Submit the job to the SLURM scheduler:

    ```bash
    sbatch sra-download.slurm
    ```

## Check Job Status:

1. Monitor the status of your job:

    ```bash
    squeue -u your-username
    ```

## Review the Output:

1. Once the job completes, check the output log file:

    ```bash
    cat /ocean/projects/agr250001p/your-username/download_sra.log
##IF YOU GET AN SRA CONFIGURATION ISSUE

1. Run
``` bash
vdb-config --interactive
```

