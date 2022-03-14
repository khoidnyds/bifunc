import subprocess
import logging
from pyfaidx import Fasta
from pathlib import Path


class Preprocess():
    def __init__(self, input_dir, output):
        self.input_dir = input_dir
        self.output = output

    def execute(self):
        reference = '/beegfs/projects/ciwars/databases/contaminants_database/merged_ref_5081444234059102403.fa.gz'

        logging.info(
            f"Run the QC steps: remove adapters and low-quality reads")

        files = sorted([x for x in self.input_dir.glob(
            '*') if x.is_file() and ".fastq" in x.suffixes])

        if len(files) % 2 != 0:
            raise Exception("Miss pair-end files")

        pairs = zip(files[::2], files[1::2])

        pairs_1 = []
        pairs_2 = []
        count = 0
        for pair in pairs:
            # if count == 2:
            #     break
            out_f_1 = self.output.joinpath(
                Path(f"fastp_{pair[0].name}"))
            out_f_2 = self.output.joinpath(
                Path(f"fastp_{pair[1].name}"))
            out_b_1 = self.output.joinpath(
                Path(f"bbduk_{pair[0].name}"))
            out_b_2 = self.output.joinpath(
                Path(f"bbduk_{pair[1].name}"))
            out_stat = self.output.joinpath(
                Path(f"bbduk_stat"))
            # fastp: remove adapter and read with quality score < 10
            logging.info(f"Fastp of {pair[0]} and {pair[1]}")
            # subprocess.run(
            #     ["fastp", "-i", pair[0], "-I", pair[1],  "-o", out_f_1, "-O", out_f_2, "--average_qual", "10"])
            # bbduk
            # print(
            #     f"bbduk.sh ref={reference} in={out_f_1} in2={out_f_2} out={out_b_1} out2={out_b_2} k=31 hdist=1 ftm=5 stats={out_stat}")
            # subprocess.run(
            #     f"bbduk.sh ref={reference} in={out_f_1} in2={out_f_2} out={out_b_1} out2={out_b_2} k=31 hdist=1 ftm=5 stats={out_stat}")

            pairs_1.append(str(out_f_1))
            pairs_2.append(str(out_f_2))
        subprocess.run(["megahit", "-1", ",".join(pairs_1), "-2",
                        ",".join(pairs_2), "--presets", "meta-large", "-o", str(self.output.joinpath("megahit"))])

        out_path = self.output.joinpath("megahit").joinpath("final.contigs.fa")
        count = len(Fasta(str(out_path)))
        logging.info(f"Found {count} contigs")
        return out_path
