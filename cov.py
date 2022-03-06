from Bio.Seq import Seq
import subprocess
from pathlib import Path
import pandas as pd
from pyfaidx import Fasta

database_path = "database/CARD/broadstreet/nucleotide_fasta_protein_homolog_model.fasta"
query_path = "results/Hospital_Wastewater/03-01--16-08-45/orfs.fa"
out = Path("cov").joinpath("hospital_water")
Path.mkdir(out, parents=True, exist_ok=True)

subprocess.run(["bowtie2-build", database_path,
               str(out.joinpath("bowtie_db_build"))])
subprocess.run(["bowtie2", "-x", str(out.joinpath("bowtie_db_build")), "-f", "-U", query_path,
               "-S", str(out.joinpath("bowtie2_aligned.sam")), "--end-to-end"])
subprocess.run(["samtools", "view", "-S", "-b", str(out.joinpath("bowtie2_aligned.sam")),
               "-o", str(out.joinpath("bowtie2_aligned.bam"))])

subprocess.run(["samtools", "sort", "-o", str(out.joinpath(
    "bowtie2_aligned_sorted.bam")), str(out.joinpath("bowtie2_aligned.bam"))])

subprocess.run(["samtools", "coverage", str(out.joinpath(
    "bowtie2_aligned_sorted.bam")), "-o", str(out.joinpath("bowtie2_aligned.cov"))])
subprocess.run(["samtools", "coverage", "-m", str(out.joinpath(
    "bowtie2_aligned_sorted.bam")), "-o", str(out.joinpath("bowtie2_aligned.his"))])
