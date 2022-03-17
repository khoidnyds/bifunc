import subprocess
import logging
from pyfaidx import Fasta
from pathlib import Path


class OrfFinder():
    def __init__(self, query):
        self.query = query

    def orf_finder(self):
        out_path = self.query.parent.parent.joinpath(Path("orfs.fa"))

        logging.info(f"Run the ORF finder")
        subprocess.run(
            ["orfipy", self.query, "--dna", "orfs.fa", "--outdir", self.query.parent.parent])

        count = len(Fasta(str(out_path)).keys())
        logging.info(f"Found {count} ORFs")
        if count == 0:
            raise Exception("No ORFs")

        return out_path
