from sklearn.cluster import DBSCAN
from pathlib import Path
import pandas as pd


class ARGSCAN():
    def __init__(self, alignment, database, out):
        self.alignment_path = Path(alignment)
        self.database_path = Path(database)
        self.out = Path(out)

    def fit(self):

        labels = DBSCAN(eps=0.5, min_samples=5).fit()

    def get_bifunctional_arg(self):
        pass


ARGSCAN("results/ERR3307076/11-12--12-05-55/aligned.tsv",
        "database/CARD/broadstreet/nucleotide_fasta_protein_homolog_model.fasta")
