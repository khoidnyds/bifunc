import subprocess
import logging
from pathlib import Path
from pyfaidx import Fasta


class OrfFinder():
    def __init__(self, query):
        self.query = query

    def orf_finder(self):
        out_path = self.query.parent.parent.joinpath(Path("orfs.fa"))
        return out_path

        logging.info(f"Run the ORF finder")
        subprocess.run(
            ["orfipy", self.query, "--dna", "orfs.fa", "--outdir", self.query.parent.parent, "--min", "200"])

        count = len(Fasta(str(out_path)).keys())
        logging.info(f"Found {count} ORFs")
        if count == 0:
            raise Exception("No ORFs")

        return out_path
