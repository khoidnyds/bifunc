from Bio import SeqIO
import csv


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
