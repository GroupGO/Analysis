"""
author = Ronald de Jongh
WUR #: 930323409080

    [1] out_name: name of the file to be output (will add .fpkm_table)
    optional:
        [2] data_path: path to the genes.attr_table (cuffnorm output)
        [3] fpkm_path: path to the fpkm table (cuffnorm output)
        [4] genome_path_annotation: path to the genome annotation file (gtf/gff format)
        [5] samples_table_path: path to the samples table (cuffnorm output)

"""
from sys import argv
from sys import stdout
import os
import pickle

def get_names(samples_table_path):
    """
    This script will take in sample_table file and output a dictionary of codes with names

    :param samples_table_path: opened file to the samples_table
    :return: name_dict: q_id => sample_name
    """
    name_dict = {}  # q_id => path
    samples_table_file = open(samples_table_path)
    for line in samples_table_file:
        line_list = line.split()
        name = line_list[1].split('/')[-1]
        parsed_name = name.split('.')[0]
        name_dict[line_list[0]] = parsed_name

    samples_table_file.close()
    return name_dict


def get_cro_id(data_path, genome_path_annotation, time=False):
    """
    This script will take in a list of XLOC id's and output a file with sequences that belong to those ID's

    :param data_path: path to the genes.attr_table file (from a cuffnorm run)
    :param genome_path_annotation: path to the annotations from the genome file.
    :param time: default False, if True shows the progress to the output
    """

    print "Opening attributes file"
    attr_file = open(data_path)
    loc_dict = {}
    for line in attr_file:
        line = line.strip()
        line_list = line.split('\t')
        loc_dict[line_list[6]] = line_list[0]  # {pos => xloc}
    attr_file.close()
    print "Opening annotation file"
    new_dict = {}
    annotation_file = open(genome_path_annotation)
    count = 0
    for line in annotation_file:
        line = line.strip()
        line_list = line.split('\t')
        if line.startswith("#"):
            continue
        start = int(line_list[3]) - 1  # index error weirdness between files
        end = int(line_list[4])
        pos = line_list[0] + ":" + '%s-%s' % (start, end)
        if pos in loc_dict.keys() and line_list[2] == 'mRNA':
            count += 1
            if count % 50 == 0 and time:
                stdout.write("XLOC codes processed: {0}\r".format(":".join(map(str, [(count), len(loc_dict.keys())]))))
                stdout.flush()
            cro_id = line_list[8].split(';')[0].split('=')[1]
            new_dict[loc_dict[pos]] = cro_id  # {xloc => cro_id}
    annotation_file.close()
    return new_dict


def write_fpkm_file(xloc_dict, fpkm_path, out_name, name_dict, time=False):
    """
    writes the fpkm file with a header and

    :param xloc_dict: dict of xloc => cro_id
    :param fpkm_path: path to the fpkm table file
    :param out_name: name for the file to be output
    :param name_dict: dict of q_id => sample_name
    :param time: default False, if True shows the progress to the output
    """
    print "Opening fpkm file"
    fpkm_file = open(fpkm_path)
    outfile = open(out_name + ".fpkm_table", "w")
    count = 0
    for line in fpkm_file:
        count += 1
        line = line.strip()
        line_list = line.split('\t')
        xloc = line_list[0]
        if line.startswith('tracking_id'):  # Hardcoded for dealing with cuffnorm output
            line_correct = line_list
            for index in range(1, len(line_list[1:])):
                line_correct[index] = name_dict[line_list[index]]
            outfile.write('tracking_id\t{0}\n'.format("\t".join(line_correct)))
        if xloc in xloc_dict.keys():
            count += 1
            if count % 50 == 0 and time:
                stdout.write("fpkm lines processed: {0}\r".format(":".join(map(str, [(count), len(xloc_dict.keys())]))))
                stdout.flush()
            outfile.write(xloc_dict[xloc] + "\t" + "\t".join(line_list[1:]) + "\n")
        else:
            print "This XLOC code was not in the dict: ", xloc
            check_raw = raw_input("Does that make sense? Y/N")
            if check_raw == "Y" or check_raw == "y":
                outfile.write(xloc_dict[xloc] + "\t" + "\t".join(line_list[1:]) + "\n")
            else:
                quit()
    print "Wrote fpkm file"
    outfile.close()
    fpkm_file.close()


def main(out_name, data_path, fpkm_path, genome_path_annotation, samples_path):
    if not os.path.exists('./name_dict_picklejar'):
        name_dict = get_names(samples_path)
        name_dict_pickle = open('name_dict_picklejar', 'wb')
        pickle.dump(name_dict, name_dict_pickle)
        name_dict_pickle.close()
    else:
        name_dict_pickle = open('name_dict_picklejar', 'rb')
        name_dict = pickle.load(name_dict_pickle)

    if not os.path.exists('./cro_id_picklejar'):
        cro_dict = get_cro_id(data_path, genome_path_annotation)
        cro_dict_pickle = open('cro_id_picklejar', 'wb')
        pickle.dump(cro_dict, cro_dict_pickle)
        cro_dict_pickle.close()
    else:
        cro_dict_pickle = open('cro_id_picklejar', 'rb')
        cro_dict = pickle.load(cro_dict_pickle)

    write_fpkm_file(cro_dict, fpkm_path, out_name, name_dict)


if __name__ == '__main__':
    assert len(argv) > 1
    out_name = argv[1]
    try:
        data_path = argv[2]
    except IndexError:
        data_path = '/local/data/BIF30806_2015_2/project/groups/go/Data/Tests/Cuff_test_Data/Cuffnorm_Data/Test/' \
                    'genes.attr_table'
        print "Using default data path!"
    try:
        fpkm_path = argv[3]
    except IndexError:
        fpkm_path = '/local/data/BIF30806_2015_2/project/groups/go/Data/Tests/Cuff_test_Data/Cuffnorm_Data/Test/' \
                    'genes.fpkm_table'
        print "Using default fpkm path!"
    try:
        genome_path_annotation = argv[4]
    except IndexError:
        genome_path_annotation = '/local/data/BIF30806_2015_2/project/genomes/Catharanthus_roseus/' \
                                 'cro_std_maker_anno.final.gff3'
        print "Using default annotation path!"
    try:
        samples_table_path = argv[5]
    except IndexError:
        samples_table_path = '/local/data/BIF30806_2015_2/project/groups/go/Data/Tests/Cuff_test_Data/Cuffnorm_Data/' \
                             'Test/samples.table'
        print "Using default annotation path!"

    assert os.path.exists(data_path)
    assert os.path.exists(fpkm_path)
    assert os.path.exists(genome_path_annotation)
    assert os.path.exists(samples_table_path)

    main(out_name, data_path, fpkm_path, genome_path_annotation, samples_table_path)
