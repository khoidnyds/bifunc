from Bio.Seq import Seq
import subprocess
from pathlib import Path
import pandas as pd
from pyfaidx import Fasta

i1_path = "data/Hospital_Wastewater/SRR17068071_1.fastq"
i2_path = "data/Hospital_Wastewater/SRR17068071_2.fastq"
contig_path = "results/Hospital_Wastewater/03-01--16-08-45/megahit/final.contigs.fa"
out = Path("cov").joinpath("hospital_water")
Path.mkdir(out, parents=True, exist_ok=True)

subprocess.run(["bowtie2-build", contig_path,
               str(out.joinpath("bowtie_db_build"))])

subprocess.run(["bowtie2", "-p", "32", "-x", str(out.joinpath("bowtie_db_build")),
               "-1", i1_path, "-2", i2_path, "-S", str(out.joinpath("bowtie2_aligned.sam"))])

subprocess.run(["samtools", "view", "-S", "-b", str(out.joinpath("bowtie2_aligned.sam")),
               "-o", str(out.joinpath("bowtie2_aligned.bam"))])

subprocess.run(["samtools", "sort", str(out.joinpath(
    "bowtie2_aligned.bam")), "-o", str(out.joinpath("bowtie2_aligned_sorted.bam"))])

subprocess.run(["samtools", "index", str(
    out.joinpath("bowtie2_aligned_sorted.bam"))])

subprocess.run(["samtools", "depth", str(out.joinpath("bowtie2_aligned_sorted.bam")),
               "-o", str(out.joinpath("bowtie2_aligned_sorted.coverage"))])

subprocess.run(["samtools", "coverage", str(out.joinpath(
    "bowtie2_aligned_sorted.bam")), "-o", str(out.joinpath("bowtie2_aligned.cov"))])

subprocess.run(["samtools", "coverage", "-m", str(out.joinpath(
    "bowtie2_aligned_sorted.bam")), "-o", str(out.joinpath("bowtie2_aligned.his"))])
