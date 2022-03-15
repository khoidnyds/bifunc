from subprocess import Popen, PIPE
from Bio import SeqIO
import csv
import logging


def read_fastq_length(filename: str):
    # read fasta files
    records = {}
    for record in SeqIO.parse(filename, "fastq"):
        length = int(record.description[record.description.find("=")+1:])
        records[record.id] = length
    return records


def write_list_to_file(filename, results: list):
    # save list to file
    my_file = open(filename, 'w')
    for element in results:
        my_file.write(str(element))
        my_file.write('\n')
    my_file.close()


def write_dict_to_file(filename, results: dict):
    # save dict to file
    my_file = open(filename, "w")
    writer = csv.writer(my_file)
    for key, value in results.items():
        writer.writerow([key, value])
    my_file.close()


def run_subprocess(cmd: str, get_stdout=False):
    if get_stdout:
        p = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)
        out, err = p.communicate()
        out = out.decode().strip()
        err = err.decode().strip()
        if out != "":
            return out
        elif err != "":
            return err
        else:
            return ""
    else:
        p = Popen(cmd, shell=True)
        p.wait()
