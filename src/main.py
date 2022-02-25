import argparse
import logging

from alignment import Alignment
from clustering import Clustering
from pathlib import Path
from myLog import Log
from datetime import datetime
from visualization import Visualization
from orf import OrfFinder
import time
from pyfaidx import Fasta
from preprocessing import Preprocess


class Bifunc():
    """
    Main class: take the input directory and run the pipeline
    """

    def __init__(self, input_dir, database, annotation, subject_cover, distance):
        self.input = input_dir

        query = Preprocess.execute(input_dir)
        self.query = Path(query)
        self.database = database
        self.distance = distance
        self.annotation = annotation
        self.subject_cover = subject_cover

        out_dir = Path("results")\
            .joinpath(self.query.stem)\
            .joinpath(datetime.today().strftime("%m-%d--%H-%M-%S"))
        Path.mkdir(out_dir, parents=True, exist_ok=True)

        start = time.time()
        self.out_alignment = out_dir.joinpath("aligned.tsv")
        self.out_clustering = out_dir.joinpath("clusters.csv")
        self.out_visualization = out_dir.joinpath("viz")
        self.out_orf = out_dir
        self.pipeline()

        running_time = time.time() - start
        logging.info(
            f"Total running time: {int(running_time//60)}:{int(running_time%60):0>2d}\n")

    def pipeline(self):
        query = Fasta(str(self.query))
        logging.info(f"Querry file has {len(query)} reads")
        try:
            OrfFinder(str(self.query), str(self.out_orf)).orf_finder()
            orf = self.out_orf.joinpath("orfs.fa")
            Alignment(str(orf), self.database, self.subject_cover,
                      self.out_alignment).align()

            Clustering(str(orf), self.out_alignment,
                       self.distance, self.out_clustering).cluster()

            Visualization(str(orf), self.out_clustering,
                          self.out_visualization, self.annotation).generate_graph()
        except Exception as e:
            logging.info(e)


def arg_parse():
    """
    Parsing arguments function
    """
    parser = argparse.ArgumentParser(
        description='Finding bifunctional ARGs')
    parser.add_argument('-i', '--input',
                        type=str,
                        required=True,
                        help='path to sequence file')
    parser.add_argument('-d', '--database',
                        type=str,
                        default="database/CARD/broadstreet/protein_fasta_protein_homolog_model.fasta",
                        help='path to database file (default: %(default)s)')
    parser.add_argument('-a', '--annotation',
                        type=str,
                        default="database/CARD/broadstreet/aro_categories_index.tsv",
                        help='path to annotation file (default: %(default)s)')
    parser.add_argument('-sc', '--subjectcover',
                        type=str,
                        default=70,
                        help='subject cover parameter for alginer (default: %(default)s)')
    parser.add_argument('-s', '--distance',
                        type=int,
                        default=20,
                        help='minimal distance (default: %(default)s)')
    parser.add_argument('-l', '--log',
                        type=str,
                        default="log",
                        help='log directory (default: %(default)s)')
    parser.add_argument('-o', '--output',
                        type=str,
                        default="results",
                        help='output directory (default: %(default)s)')
    args = parser.parse_args()
    return args


def main(args):
    """
    Main function
    """
    # logging set up
    logging = Log(path=Path(args.log).joinpath(Path(args.input).stem))
    logging.info(f"Running {args.input} against {args.database}")
    logging.info("----- Parameters -----")
    logging.info(f"Query file: {args.input}")
    logging.info(f"Database file: {args.database}")
    logging.info(f"Annotation file: {args.annotation}")
    logging.info(f"Subject cover: {args.subjectcover}")
    logging.info(f"Distance: {args.distance}")
    logging.info(
        f"Logging directory: {Path(args.log).joinpath(Path(args.input).stem)}")
    logging.info(
        f"Results directory: {Path(args.output).joinpath(Path(args.input).stem)}")
    logging.info("----------------------")

    # report file set up
    outfile = Path(args.output)
    Path.mkdir(outfile, parents=True, exist_ok=True)

    Bifunc(args.input, args.database, args.annotation,
           args.subjectcover, args.distance)
    ###################################################################


class Range(object):
    def __init__(self, start, end):
        self.start = start
        self.end = end

    def __eq__(self, other):
        return self.start <= other <= self.end


if __name__ == "__main__":
    args = arg_parse()
    main(args)
