"""
author = Ronald de Jongh
WUR #: 930323409080

Notes: currently only writes files to the current directory!

input:
    [1] out_name: name for the outfiles (will add .fasta and .fpkm)
    [2] in_file: name for the infile
    Optional:
        [3] data_path: path to the genes.attr_table file (from a cuffnorm run)
        [4] genome_path_transcript: path to the genome folder that contains a gff/gtf transcript file
        [5] genome_path_annotation: path to the genome folder that contains a fasta annotation file
        [6] fpkm_path: path to the genome folder that contains the fpkm_table file
output:
    fasta file: >cro_id xloc_id scaffold:position \n sequence \n etc...

Usage: python id_finder.py out_name in_file [optional: data_path genome_path_transcript genome_path_annotation]
"""
from sys import argv
import os


def parse_in_file(infile_name):
    """
    Returns a list of XLOC id's from a file

    :param infile_name: the name of the infile, should be newline delimited XLOC id's
    :return: a python list of files
    """
    out_list = []
    infile = open(infile_name)
    for line in infile:
        line = line.strip()
        out_list.append(line)
    return out_list


def get_cro_id(data_path, genome_path_annotation):
    """
    This script will take in a list of XLOC id's and output a file with sequences that belong to those ID's

    :param out_name: name for the outfile
    :param in_list: name for the infile
    :param data_path: path to the genes.attr_table file (from a cuffnorm run)
    :param genome_path_transcript: path to the transcripts from the genome file.
    :param genome_path_annotation: path to the annotations from the genome file.
    """

    print "Opening attributes file"
    attr_file = open(data_path)
    loc_dict = {}
    for line in attr_file:
        line = line.strip()
        line_list = line.split('\t')
        if line_list[0] in in_list:
            loc_dict[line_list[6]] = line_list[0]  # {pos => xloc}
    attr_file.close()
    print "Opening annotation file"
    new_dict = {}
    annotation_file = open(genome_path_annotation)
    for line in annotation_file:
        line = line.strip()
        line_list = line.split('\t')
        if line.startswith("#"):
            continue
        start = int(line_list[3]) - 1  # index error weirdness between files
        end = int(line_list[4])
        pos = line_list[0] + ":" + '%s-%s' % (start, end)
        if pos in loc_dict.keys() and line_list[2] == 'mRNA':
            cro_id = line_list[8].split(';')[0].split('=')[1]
            new_dict[cro_id] = [pos, loc_dict[pos]] # {cro_id => [pos, xloc]}
    annotation_file.close()
    return new_dict

def write_fasta_file(genome_path_transcript, out_name, pos_dict):
    print "Opening transcripts file"
    transcripts_file = open(genome_path_transcript)
    outfile = open(out_name+".fasta", "w")
    record = False
    for line in transcripts_file:
        line = line.strip()
        cro_id = line[1:]
        if record and line.startswith(">"):
            record = False
        if cro_id in pos_dict.keys() and line.startswith(">"):
            record = True
            outfile.write(">%s %s %s\n" % (cro_id, pos_dict[cro_id][0], pos_dict[cro_id][1]))
        elif record:
            outfile.write(line + "\n")
    print "Wrote fasta file"
    outfile.close()
    transcripts_file.close()


def write_fpkm_file(pos_dict, fpkm_path, out_name):
    print "Opening fpkm file"
    fpkm_file = open(fpkm_path)
    outfile = open(out_name+".fpkm", "w")
    xloc_dict = {}  # {xloc => cro_id}
    for key in pos_dict:
        xloc_dict[pos_dict[key][1]] = key
    for line in fpkm_file:
        line = line.strip()
        line_list = line.split('\t')
        xloc = line_list[0]
        if line.startswith('tracking_id'):  # Hardcoded for dealing with cuffnorm output
            outfile.write(line+ "\n")
        if line_list[0] in xloc_dict.keys():
            outfile.write(xloc_dict[line_list[0]]+"\t"+"\t".join(line_list[1:])+"\n")
    print "Wrote fpkm file"
    outfile.close()
    fpkm_file.close()


def main(out_name, in_list, data_path, fpkm_path, genome_path_transcript, genome_path_annotation):
    """
    This function was designed to be importable, it will write a fasta

    **DEPRECATED**
    And a transformed fpkm table to the current directory.
    **************

    :param out_name: name for the outfile
    :param in_list: name for the infile
    :param data_path: path to the genes.attr_table file (from a cuffnorm run)
    :param fpkm_path: path the genes.count_table file (from a cuffnorm run)
    :param genome_path_transcript: path to the transcripts from the genome file.
    :param genome_path_annotation: path to the annotations from the genome file.
    """

    pos_dict = get_cro_id(data_path, genome_path_annotation)
    write_fasta_file(genome_path_transcript, out_name, pos_dict)
    #write_fpkm_file(pos_dict, fpkm_path, out_name)


if __name__ == '__main__':
    assert len(argv) > 1
    out_name = argv[1]
    infile_name = argv[2]

    try:
        data_path = argv[3]
    except IndexError:
        data_path = '/local/data/BIF30806_2015_2/project/groups/go/Data/Tests/Cuff_test_Data/Cuffnorm_Data/Test/' \
                    'genes.attr_table'
        print "Using default data path!"
    try:
        genome_path_transcript = argv[4]
    except IndexError:
        genome_path_transcript = '/local/data/BIF30806_2015_2/project/genomes/Catharanthus_roseus/' \
                                 'cro_std_maker_anno.final.transcripts.fasta'
        print "Using default transcript path!"
    try:
        genome_path_annotation = argv[5]
    except IndexError:
        genome_path_annotation = '/local/data/BIF30806_2015_2/project/genomes/Catharanthus_roseus/' \
                                 'cro_std_maker_anno.final.gff3'
        print "Using default annotation path!"
    try:
        fpkm_path = argv[6]
    except IndexError:
        fpkm_path = '/local/data/BIF30806_2015_2/project/groups/go/Data/Tests/Cuff_test_Data/Cuffnorm_Data/Test/' \
                    'genes.fpkm_table'
        print "Using default fpkm path!"
    assert os.path.exists(data_path)
    assert os.path.exists(fpkm_path)
    assert os.path.exists(genome_path_transcript)
    assert os.path.exists(genome_path_annotation)
    # assert os.path.exists(infile_name)

    in_list = parse_in_file(infile_name)
    main(out_name, in_list, data_path, fpkm_path, genome_path_transcript, genome_path_annotation)
