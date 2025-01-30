https://www.ncbi.nlm.nih.gov/bioproject/59771

Download human genome DNA sequence when FTP link is available
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.fna.gz

wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.gff.gz #Downloads genome annotations

wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_protein.faa.gz #Download human proteome
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_translated_cds.faa.gz #download coding sequences only

When you do not have and FTP link, but an SRA accession.
module load sra-toolkit  #the sratoolkit is installed in the HPC, but before we can use it, we should call it and load it

slurm #!/bin/bash
#SBATCH --job-name=download_sra
#SBATCH --output=download_sra.log
#SBATCH --time=00:10:00          # Adjust time limit as needed
#SBATCH --ntasks=1
#SBATCH --mem=2G                 # Minimal memory required
#SBATCH --partition=standard     # Use a partition with internet access

# Load the SRA Toolkit module
module load sratoolkit

# Download the dataset
fastq-dump --split-files SRX123456

wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR741/SRR741411/SRR741411_2.fastq.gz #download lowpass genome data from human female British

1000 genomes project: sra accessions, explore each to find most interesting

https://www.ncbi.nlm.nih.gov/sra?linkname=bioproject_sra_all&from_uid=41223


wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.fna.gz

gunzip GCF_000001405.40_GRCh38.p14_cds_from_genomic.fna.gz

module load BLAST

makeblastdb -in GCF_000001405.40_GRCh38.p14_cds_from_genomic.fna -dbtype nucl -out human_cds_db

#Locate unknown sequence with your number
blastn -query unk.fasta -db human_genome -out results.txt -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle"


https://useast.ensembl.org/index.html

Tyoe


samtools faidx GCF_000001405.40_GRCh38.p14_genomic.fna
samtools faidx GCF_000001405.40_GRCh38.p14_genomic.fna NC_000001.11:223795811-223796603 > mismatchregion.txt
