import subprocess
import logging
from pyfaidx import Fasta
from pathlib import Path


class Preprocess():
    def __init__(self, input_dir, output):
        self.input_dir = input_dir
        self.output = output

    def execute(self):
        logging.info(f"Run the QC steps: remove adapter and bad read")
        subprocess.run(
            ["orfipy", self.query, "--dna", "orfs.fa", "--outdir", self.out])
        count = len(Fasta(f"{self.out}/orfs.fa"))
        logging.info(f"Found {count} ORFs")
        if count == 0:
            raise Exception("No ORFs")
