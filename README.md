# Detecting bifunctional antibiotic resistance genes in metagenomics data
## Master thesis at Dr.Zhang's lab - Computer Science - Virginia Tech - May 2022

### Install and activate the environment:
    conda env create -f environment.yml
    conda activate bifunc
### Get database
    wget https://bifunc.s3.amazonaws.com/database.tar.gz
    tar -xvf database.tar.gz

### Using
    -i or --input: Path to input directory (fastq or fasta files)
    -d or --database: Path to database file (fastq or fasta files)
    -a or --annotation: Path to annotation file (fastq or fasta files)
    -sc or --subjectcover: Subject cover parameter for alginer (default: 70)
    -s or --distance: Minimal distance of aligned sequence in a cluster (default: 20)
    -l or --log: Log directory (default: log/)
    -o or --output: Output directory (default: result/)