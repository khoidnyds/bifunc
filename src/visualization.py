import matplotlib.pyplot as plt
import pandas as pd
import logging
from pathlib import Path
import numpy as np
from pyfaidx import Fasta
import re
import pandas as pd
from scipy import stats
from textwrap import wrap


class Visualization():
    def __init__(self, query, clusters, annontation, coverage):
        self.annontation = annontation
        self.query = query
        self.clusters = clusters
        self.coverage = coverage

    def get_window_size(self, sample_size):
        # This size was chosen because it reduces the size, but still
        # looks nice over a large range of values.
        window_size = np.floor(sample_size / 1000 * 3)

        # Window size must be odd.
        window_size = window_size + 1 if window_size // 2 == 0 else window_size
        return window_size

    def generate_graph(self):
        name = None
        if not self.annontation is None:
            name = pd.read_csv(self.annontation, sep="\t").fillna(np.NaN)

        df = pd.read_csv(self.clusters)
        if df.empty:
            return
        logging.info(
            f"Drawing {df['Cluster'].iloc[-1]} clusters")

        out_path = self.query.parent.joinpath("viz_clusters")
        Path.mkdir(out_path, parents=True, exist_ok=True)
        orfs = Fasta(str(self.query))

        for _, cluster in df.groupby("Cluster"):
            plt.figure(figsize=(5, len(cluster)))

            # draw cluster
            contig_name = cluster["Query accession"].iloc[0]
            pos_start, pos_end = map(int, re.findall(
                r'\[.*?\]', orfs[contig_name].long_name)[0][1:-1].replace("-", ",").split(","))

            contig_key = contig_name.split("_")[:-1]
            contig_key.append("split")
            contig_key.append("00001")
            contig_key = "_".join(contig_key)

            length = cluster['Length'].iloc[0]

            # draw coverage
            cov_path = self.query.parent.joinpath(
                "viz_coverage").joinpath(contig_key+".txt")
            cov = pd.read_csv(cov_path, sep='\t')
            window_size = self.get_window_size(len(cov))

            cov = cov.groupby('nt_position', as_index=False)['coverage'].sum()
            smooth_cov_ls = []
            for row_idx, _ in cov.iterrows():
                start = 0 if row_idx < window_size//2 else row_idx
                end = start + window_size//2 + 1 if start + \
                    window_size//2 + 1 < len(cov) else len(cov)
                smooth_cov_ls.append(
                    np.mean(cov[int(start):int(end)].coverage.values))
            cov['smooth_cov'] = smooth_cov_ls

            ks_data = []
            # draw clusters
            y_max = np.max(cov['smooth_cov'][pos_start:pos_end].values)
            for i in range(len(cluster)):
                y = (i + 1)*y_max/(len(cluster)*2)
                target_start = cluster["Target start"].iloc[i]
                target_end = cluster["Target end"].iloc[i]
                x_start = cluster["Query start"].iloc[i]+pos_start
                x_end = cluster["Query end"].iloc[i]+pos_start

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

                label_wrap = '\n'.join(wrap(label, 50))
                label_wrap = f"Query location: {x_start} {x_end}\n"+label_wrap
                # labels.append(
                #     f"Target location: {target_start} {target_end}\nQuery location: {x_start} {x_end}\n{label_wrap}")
                plt.plot((x_start, x_end), (y, y), label=label_wrap)

                ks_data.append(cov[int(cluster["Query start"].iloc[i]) -
                               1:int(cluster["Query end"].iloc[i])-1]['smooth_cov'].values)

            if len(range(pos_start, pos_end)) != len(cov['smooth_cov'][pos_start:pos_end].values):
                logging.info(f"Something's wrong {contig_key}")
                continue
            plt.plot(range(pos_start, pos_end),
                     cov['smooth_cov'][pos_start:pos_end].values,  'k--',
                     alpha=0.5,
                     label="Coverage")

            ks_test = stats.ks_2samp(ks_data[0], ks_data[1])
            plt.xlim(pos_start-10, pos_end+10)
            plt.axline((0, 0), (length, 0), c='r', linewidth=4)
            mean_cov = np.mean(np.concatenate((ks_data[0], ks_data[1])))
            plt.title(
                f"{cluster['Query accession'].iloc[0]} | Length:{length} | Mean coverage:{mean_cov:.2f} | Kolmogorov-Smirnov test: value:{ks_test.statistic:.2f}, p-value:{ks_test.pvalue:.1E}")
            plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
            plt.xlabel("Query sequence")
            plt.ylabel("# of mapped read")
            plt.savefig(out_path.joinpath(
                contig_name+".png"), bbox_inches='tight')
            plt.close()
