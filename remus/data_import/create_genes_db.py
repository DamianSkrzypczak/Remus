import argparse
import os
import sqlite3

import pandas as pd

BED_FORMAT_COLUMNS_ORDER = pd.Index([
    'chrom', 'txStart', 'txEnd', 'strand', 'cdsStart', 'cdsEnd', 'exonCount',
    'exonStarts', 'exonEnds', 'name', 'geneSymbol'
])
SEPARATOR = "\t"

EXPECTED_FILE_EXTENSION = ".tsv"

CANONICAL_CHROMOSOMES = ['chr'+str(i+1) for i in range(22)] + ['chrM', 'chrX', 'chrY']


def main(data_path, db_path):
    """
    Main function of script.

    :param data_path: directory containing *.tsv source files (script will try to use every *.tsv file)
    :param db_path: full (including filename) destination path of result database
    """
    source_files = [filename for filename in os.listdir(data_path) if filename.endswith(EXPECTED_FILE_EXTENSION)]
    genomes_frames = [extract_data_frame_from_file(data_path, filename) for filename in source_files]
    frames_files_pairs = zip(genomes_frames, source_files)
    reformatted_data_frames = [reformat_data_frame(data_frame, filename) for data_frame, filename in frames_files_pairs]
    filtered_frames = [filter_transcripts(df) for df in reformatted_data_frames]
    summary_df = pd.concat(filtered_frames, ignore_index=True)
    create_db_from_data_frame(db_path, summary_df)


def extract_data_frame_from_file(data_path, filename):
    """
    :param data_path: path of source file directory
    :param filename: source file name
    :return: pandas dataframe with file content
    """
    file_path = os.path.join(data_path, filename)
    return pd.read_csv(file_path, skiprows=0, sep=SEPARATOR)


def filter_transcripts(df):
    """
    :param df: data frame with gene coordinates
    :return: data frame with filtered gene coordinates
    """
    #return prune_noncanonical_transcripts(prune_noncanonical_chromosomes(df))
    return prune_noncanonical_chromosomes(df)


def prune_noncanonical_transcripts(df):
    """
    :param df: data frame with gene coordinates
    :return: data frame with canonical transcripts only - identified by last column being not NA
    """
    # select gene symbols with and without canonical transcripts
    with_canonical = df.loc[df["name"]==df["transcript"]]

    withno_canonical = df.loc[~df["geneSymbol"].isin(with_canonical["geneSymbol"])
                              & df["transcript"].isna()]

    return with_canonical.append(withno_canonical)


def prune_noncanonical_chromosomes(df):
    """
    :param df: data frame with gene coordinates
    :return: data frame with gene coordinates on canonical chromosomes (1-22, chrX, chrY, chrM)
    """
    chr_column = 'chrom'
    return df.loc[df[chr_column].isin(CANONICAL_CHROMOSOMES)]


def reformat_data_frame(data_frame, filename):
    """
    :param data_frame: pandas dataframe
    :param filename: name of source file (WARNING!!!: FILE SHOULD BE NAMED AFTER GENOME, FOR EXAMPLE "hg37.tsv")
    :return: reformatted dataframe
    """
    reformat_column_names(data_frame)
    data_frame = order_columns_for_bed(data_frame)
    data_frame = set_genome_value(filename, data_frame)
    return data_frame


def reformat_column_names(single_genome_data_frame):
    """
    Function removing prefixes from column names,
    in this case, prefix is everything from start to last dot of name.

    :param single_genome_data_frame: source dataframe
    """
    single_genome_data_frame.columns = [col_name.split(".")[-1] for col_name in single_genome_data_frame.columns]


def order_columns_for_bed(single_genome_data_frame):
    """
    Orders dataframe to provide compatibility with BED format

    :param single_genome_data_frame: unordered dataframe
    :return: ordered dataframe
    """
    return single_genome_data_frame[BED_FORMAT_COLUMNS_ORDER]


def set_genome_value(filename, single_genome_data_frame):
    """
    Adding genome information into table.

    :param filename: name of source file (WARNING!!!: FILE SHOULD BE NAMED AFTER GENOME, FOR EXAMPLE "hg38.tsv")
    :param single_genome_data_frame: source dataframe
    """
    return single_genome_data_frame.assign(genome = os.path.splitext(filename)[0])


def create_db_from_data_frame(db_path, summary_df):
    """
    Creates database from  whole content which was collected and reformatted as dataframe.

    :param db_path: full path of future database
    :param summary_df: content for future database
    """
    conn = sqlite3.connect(db_path)
    summary_df.to_sql("genes", conn, if_exists="replace")


def get_input():
    """
    Very simple CLI provider.
    :return: user input
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("-i")
    parser.add_argument("-o")
    return parser.parse_args()


if __name__ == '__main__':
    args = get_input()
    print("input dir {} output file {}".format(args.i, args.o))
    main(args.i, args.o)
