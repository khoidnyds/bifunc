import matplotlib.pyplot as plt
import pandas as pd
import logging
from pathlib import Path
import numpy as np


class Visualization():
    def __init__(self, query, clusters, annontation):
        self.annontation = annontation
        self.query = query
        self.clusters = clusters

    def generate_graph(self):
        name = None
        if not self.annontation is None:
            name = pd.read_csv(self.annontation, sep="\t").fillna(np.NaN)

        df = pd.read_csv(self.clusters)
        if df.empty:
            return
        logging.info(
            f"Drawing {df['Cluster'].iloc[-1]} clusters")

        out_path = self.query.parent.joinpath("viz")
        Path.mkdir(out_path, parents=True, exist_ok=True)

        for cluster_idx, cluster in df.groupby("Cluster"):
            length = cluster['Length'].iloc[0]
            plt.figure(figsize=(5, len(cluster)))
            plt.xlim(0, length)
            plt.axline((0, 0), (length, 0), c='r', linewidth=4)
            plt.title(
                f"{cluster['Query accession'].iloc[0]} - Length:{length}")
            locs, labels = [], []
            for i in range(len(cluster)):
                y = i + 1
                locs.append(y)
                target_start = cluster["Target start"].iloc[i]
                target_end = cluster["Target end"].iloc[i]
                x_start = cluster["Query start"].iloc[i]
                x_end = cluster["Query end"].iloc[i]

                raw_label = cluster['Target accession'].iloc[i]
                drug_class = name[name['Protein Accession'] == raw_label.split(
                    '|')[1]]['Drug Class'].values[0]
                if drug_class is np.NaN:
                    arg_name = "N/A"
                else:
                    arg_name = drug_class.replace(" antibiotic", "")
                label = raw_label.split("|")[-2:]
                label.append(arg_name)
                label = " ".join(label)

                labels.append(
                    f"Target location: {target_start} {target_end}\nQuery location: {x_start} {x_end}\n{label}")
                plt.plot((x_start, x_end), (y, y))
            plt.yticks(locs, labels)
            plt.savefig(out_path.joinpath(
                f"cluster_{cluster_idx}.png"), bbox_inches='tight')
            plt.close()


