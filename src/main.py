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
from preprocessing import Preprocess


class Bifunc():
    """
    Main class: take the input directory and run the pipeline
    """

    def __init__(self, input_dir, database, annotation, subject_cover, distance):
        self.input = Path(input_dir)
        self.database = Path(database)
        self.distance = distance
        self.annotation = annotation
        self.subject_cover = subject_cover

        today = datetime.today().strftime("%m-%d--%H-%M-%S")
        today = "03-14--19-18-25"
        self.out_dir = Path("results")\
            .joinpath(self.input.stem)\
            .joinpath(today)
        Path.mkdir(self.out_dir, parents=True, exist_ok=True)

        start = time.time()
        self.pipeline()

        running_time = time.time() - start
        logging.info(
            f"Total running time: {int(running_time//60)}:{int(running_time%60):0>2d}\n")

    def pipeline(self):
        try:
            query = Preprocess(self.input, self.out_dir).execute()
            orf = OrfFinder(query).orf_finder()
            aligned = Alignment(orf, self.database, self.subject_cover).align()
            clusters = Clustering(orf, aligned, self.distance).cluster()
            Visualization(orf, clusters, self.annotation).generate_graph()
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
