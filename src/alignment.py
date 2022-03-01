import subprocess
from pathlib import Path
import logging
import pandas as pd


class Alignment():
    def __init__(self, query, database, subject_cover):
        self.query = query
        self.subject_cover = str(subject_cover)
        self.database = database

    def align(self):
        out_path = self.query.parent.joinpath("aligned.tsv")

        database_built_path = Path("database").joinpath(
            "built").joinpath(self.database.stem)
        if not Path(str(database_built_path)+".dmnd").is_file():
            subprocess.run(
                ["diamond", "makedb", "--in", self.database, "-d", database_built_path])
        logging.info(f"Run the DIAMOND alignment")
        subprocess.run(
            ["diamond", "blastx", "--db", database_built_path, "--query", self.query, "--verbose", "--sensitive", "--subject-cover", self.subject_cover, "-o", out_path])
        # "--outfmt", "101"

        try:
            count = len(pd.read_csv(out_path, sep='\t',  comment='@'))
            logging.info(f"Found {count} aligned sequences")
        except:
            raise Exception("No aligned sequences")

        return out_path
