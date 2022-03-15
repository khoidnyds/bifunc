import subprocess
import logging
from pyfaidx import Fasta
from pathlib import Path
from utils import run_subprocess


class Preprocess():
    def __init__(self, input_dir, output):
        self.input_dir = input_dir
        self.output = output

    def execute(self):
        reference = 'database/merged_ref_5081444234059102403.fa.gz'

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
            if count == 2:
                break
            count += 1

            out_f_1 = self.output.joinpath(
                Path(f"fastp_{pair[0].name}"))
            out_f_2 = self.output.joinpath(
                Path(f"fastp_{pair[1].name}"))
            out_f_html = self.output.joinpath(
                Path(f"fastp_{pair[0].stem}.html"))
            out_j_html = self.output.joinpath(
                Path(f"fastp_{pair[0].stem}.json"))

            out_b_at_1 = self.output.joinpath(
                Path(f"bbduk_at_{pair[0].name}"))
            out_b_at_2 = self.output.joinpath(
                Path(f"bbduk_at_{pair[1].name}"))
            out_b_at = self.output.joinpath(
                Path(f"bbduk_adapter_trimming_stats_{pair[0].stem}.txt"))

            # fastp: remove adapter and read with quality score < 10
            logging.info(f"Preprocessing {pair[0]} and {pair[1]}")
            fastp_command = f"fastp -i {pair[0]} -I {pair[1]} -o {out_f_1} -O {out_f_2} -h {out_f_html} -j {out_j_html} --average_qual 10"
            # bbduk
            bbduk_trim_adapters_command = f"bbduk.sh ref={reference} in={out_f_1} in2={out_f_2} out={out_b_at_1} out2={out_b_at_2} -Xmx500g ktrim=r stats={out_b_at}"

            run_subprocess(fastp_command)
            run_subprocess(bbduk_trim_adapters_command)

            pairs_1.append(str(out_b_at_1))
            pairs_2.append(str(out_b_at_2))

        subprocess.run(["megahit", "-1", ",".join(pairs_1), "-2",
                        ",".join(pairs_2), "--presets", "meta-large", "-o", str(self.output.joinpath("megahit"))])

        out_path = self.output.joinpath("megahit").joinpath("final.contigs.fa")
        count = len(Fasta(str(out_path)))
        logging.info(f"Found {count} contigs")
        return out_path
