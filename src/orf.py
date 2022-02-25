import subprocess
import logging
from pyfaidx import Fasta

class OrfFinder():
    def __init__(self, query, out):
        self.query = query
        self.out = out

    def orf_finder(self):
        logging.info(f"Run the ORF finder")
        subprocess.run(
            ["orfipy", self.query, "--dna", "orfs.fa", "--outdir", self.out])
        count = len(Fasta(f"{self.out}/orfs.fa"))
        logging.info(f"Found {count} ORFs")
        if count == 0:
            raise Exception("No ORFs")