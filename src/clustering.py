import pandas as pd
import numpy as np
from pyfaidx import Fasta
import logging
from tqdm import tqdm
from pathlib import Path


class Clustering():
    def __init__(self, query, alignment, distance):
        self.alignment_path = alignment
        self.query_path = query
        self.distance = distance

    def cluster(self):
        labels = ["Query accession", "Target accession", "Sequence identity", "Length", "Mismatches",
                  "Gap openings", "Query start", "Query end", "Target start", "Target end", "E-value", "Bit score"]
        query = Fasta(str(self.query_path))
        aligned = pd.read_csv(self.alignment_path, sep='\t', names=labels)

        groups = aligned.groupby("Query accession")
        results = pd.DataFrame(columns=["Cluster", "Query accession", "Length",
                                        "Query start", "Query end", "Target accession", "Target start", "Target end"])
        # iterate through read
        logging.info(f'Set distance of {self.distance}, find clusters')
        num_cluster = 1
        for group in tqdm(groups, total=len(groups)):
            # fix start and end position in order
            group = group[1]

            group['Start'] = np.where(
                group['Query start'] < group['Query end'], group['Query start'], group['Query end'])
            group['End'] = np.where(group['Query start'] < group['Query end'],
                                    group['Query end'], group['Query start'])
            group = group[["Query accession",
                           "Target accession",  "Start", "End", "Target start", "Target end"]]

            cluster = pd.DataFrame(columns=["Cluster", "Query accession", "Length",
                                            "Query start", "Query end", "Target accession", "Target start", "Target end", "Range"])

            # query seq must be hitted more than 1 times
            if len(group) > 1:
                group = group.sort_values(by=["Start"])
                group = group.reset_index(drop=True)

                start, end, start_cluster, end_cluster = 0, 0, 0, 0
                for seq in group.iterrows():
                    if cluster.shape[0] == 0:
                        start_cluster, end_cluster = seq[1]["Start"], seq[1]["End"]
                        cluster.loc[cluster.shape[0]] = [num_cluster,
                                                         seq[1]["Query accession"],
                                                         len(str(
                                                             query[seq[1]["Query accession"]])),
                                                         seq[1]["Start"],
                                                         seq[1]["End"],
                                                         seq[1]["Target accession"],
                                                         seq[1]["Target start"],
                                                         seq[1]["Target end"],
                                                         seq[1]["Start"]-seq[1]["End"]]

                    else:
                        start, end = seq[1]["Start"], seq[1]["End"]
                        # whenever new seq fall outsize the range of cluster, store cluster and create new empty one
                        if(start > end_cluster + self.distance):
                            cluster = cluster.sort_values("Range").drop_duplicates(
                                "Target accession", keep="first")
                            cluster = cluster.drop_duplicates(
                                ["Query start", "Query end"], keep="first")

                            if cluster.shape[0] > 1:
                                results = pd.concat(
                                    [results, cluster])
                                num_cluster += 1
                            cluster = cluster.iloc[0:0]
                            start_cluster, end_cluster = seq[1]["Start"], seq[1]["End"]
                            cluster.loc[cluster.shape[0]] = [num_cluster,
                                                             seq[1]["Query accession"],
                                                             len(str(
                                                                 query[seq[1]["Query accession"]])),
                                                             seq[1]["Start"],
                                                             seq[1]["End"],
                                                             seq[1]["Target accession"],
                                                             seq[1]["Target start"],
                                                             seq[1]["Target end"],
                                                             seq[1]["Start"]-seq[1]["End"]]
                            continue

                        if start > (start_cluster+(end_cluster-start_cluster)*0.9):
                            # add seq to cluster and expand the range of cluster
                            cluster.loc[cluster.shape[0]] = [num_cluster,
                                                             seq[1]["Query accession"],
                                                             len(str(
                                                                 query[seq[1]["Query accession"]])),
                                                             seq[1]["Start"],
                                                             seq[1]["End"],
                                                             seq[1]["Target accession"],
                                                             seq[1]["Target start"],
                                                             seq[1]["Target end"],
                                                             seq[1]["Start"]-seq[1]["End"]]
                            start_cluster = start if start < start_cluster else start_cluster
                            end_cluster = end if end > end_cluster else end_cluster

                        # end of group
                        if(seq[0] == group.shape[0]-1):
                            cluster = cluster.sort_values("Range").drop_duplicates(
                                "Target accession", keep="first")
                            cluster = cluster.drop_duplicates(
                                ["Query start", "Query end"], keep="first")
                            if cluster.shape[0] > 1:
                                results = pd.concat(
                                    [results, cluster])
                                num_cluster += 1
                            cluster = cluster.iloc[0:0]

        if results.empty:
            raise Exception("No clusters found")
        logging.info(
            f"Found {results['Cluster'].iloc[-1]} clusters")

        out_path = self.query_path.parent.joinpath("clusters.csv")
        results.to_csv(out_path, index=False)
        return out_path


Clustering(Path("results/ILLUMINA_LANE_1/02-28--09-38-32/orfs.fa"), Path("results/ILLUMINA_LANE_1/02-28--09-38-32/aligned.tsv"), 20).cluster()