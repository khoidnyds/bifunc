from pathlib import Path
from utils import run_subprocess
from pathlib import Path


class CoveragePlot():
    def __init__(self, input_dir, contigs_ls, output):
        self.input_dir = input_dir
        self.output = output
        self.contigs_ls = contigs_ls

    def plot(self):
        contigs_raw = self.output
        self.output = self.output.parent.parent
        contigs_fa = self.output.joinpath("contigs.fa")
        contigs_db = self.output.joinpath("contigs.db")

        profile_merged_dir = self.output.joinpath("samples_merged")
        profile_merged = profile_merged_dir.joinpath('PROFILE.db')
        plot_dir = self.output.joinpath('viz_coverage')
        Path.mkdir(plot_dir, parents=True, exist_ok=True)

        run_subprocess(
            f'reformat.sh in={contigs_raw} out={contigs_fa} trimreaddescription=t')
        run_subprocess(
            f'anvi-gen-contigs-database -f {contigs_fa} -T 32 -o {contigs_db}')

        run_subprocess(f'anvi-run-hmms -T 32 -c {contigs_db} --just-do-it')
        # run_subprocess(f'anvi-setup-ncbi-cogs -T 32')
        # run_subprocess(f'anvi-run-ncbi-cogs -T 32 -c {contigs_db}')

        files = sorted([x for x in self.input_dir.glob(
            '*') if x.is_file() and ".fastq" in x.suffixes])
        pairs = zip(files[::2], files[1::2])
        for pair in pairs:
            p1 = self.output.joinpath(pair[0].name)
            p2 = self.output.joinpath(pair[1].name)

            accession_id = "_".join(pair[0].stem.split("_")[:-1])
            fa_1 = p1.with_suffix('').with_suffix('.fa')
            fa_2 = p2.with_suffix('').with_suffix('.fa')
            sam_path = self.output.joinpath(f"{accession_id}.sam")
            bam_path = self.output.joinpath(f"{accession_id}_raw.bam")
            bam_sort_path = self.output.joinpath(f"{accession_id}.bam")

            run_subprocess(
                f'reformat.sh in={pair[0]} in2={pair[1]} out={fa_1} out2={fa_2} trimreaddescription=t')
            run_subprocess(
                f'minimap2 -ax sr {contigs_fa} {fa_1} {fa_2} -t 32 -o {sam_path}')
            run_subprocess(
                f'samtools view -h {sam_path} -@ 32 -o {bam_path}')
            run_subprocess(
                f'anvi-init-bam {bam_path} -T 32 -o {bam_sort_path}')
            run_subprocess(
                f'anvi-profile -i {bam_sort_path} -T 32 -c {contigs_db}')

        profile_ls = [str(x) for x in self.output.glob(
            '**/*') if x.is_file() and x.stem == 'PROFILE']
        run_subprocess(
            f'anvi-merge {" ".join(profile_ls)} -o {profile_merged_dir} -c {contigs_db}')

        run_subprocess(
            f'anvi-get-split-coverages -p {profile_merged} -c {contigs_db} --list-splits')

        for contig in self.contigs_ls:
            cov_txt = plot_dir.joinpath(contig).with_suffix('.txt')
            cov_pdf = plot_dir.joinpath(
                contig).with_suffix('').with_suffix('.pdf')
            run_subprocess(
                f'anvi-get-split-coverages -p {profile_merged} -c {contigs_db} --split-name {contig} -o {cov_txt}')
            run_subprocess(
                f'anvi-script-visualize-split-coverages -i {cov_txt} -o {cov_pdf}')
