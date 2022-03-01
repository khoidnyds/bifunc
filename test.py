import subprocess
from pathlib import Path

path = Path('/beegfs/projects/ciwars/ILLUMINA_LANE_1/')
files = sorted([str(x) for x in path.glob('*') if x.is_file() and ".fastq" in x.suffixes])
if len(files) % 2 != 0:
    raise Exception("Miss pair-end files")

limit = 2
pairs_1 = files[::2][:limit]
pairs_2 = files[1::2][:limit]

subprocess.run(["megahit", "-1", ",".join(pairs_1), "-2", ",".join(pairs_2), "--presets", "meta-large"])