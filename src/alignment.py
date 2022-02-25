import subprocess
from pathlib import Path
import logging
import pandas as pd
class Alignment():
    def __init__(self, query, database, subject_cover, out):
        self.query = query
        self.subject_cover = subject_cover
        self.database = database
        self.out = out

    def align(self):
        database_built_path = Path("database").joinpath(
            "built").joinpath(Path(self.database).stem)
        if not Path(str(database_built_path)+".dmnd").is_file():
            subprocess.run(
                ["diamond", "makedb", "--in", self.database, "-d", database_built_path])
        logging.info(f"Run the DIAMOND alignment")
        subprocess.run(
            ["diamond", "blastx", "--db", database_built_path, "--query", self.query, "--verbose", "--ultra-sensitive", "--subject-cover", self.subject_cover, "-o", self.out, "--outfmt", "101"])

        try:
            count = len(pd.read_csv(self.out, sep='\t'))
            logging.info(f"Found {count} aligned sequences")
        except:
            raise Exception("No aligned sequences")
        
       