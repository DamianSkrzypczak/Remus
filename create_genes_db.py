import argparse
import os
import sqlite3

import pandas as pd

HARDCODED_COLUMNS_ORDER = pd.Index(['chrom',
                                    'txStart',
                                    'txEnd',
                                    'strand',
                                    'cdsStart',
                                    'cdsEnd',
                                    'exonCount',
                                    'exonStarts',
                                    'exonEnds',
                                    'name',
                                    'geneSymbol'])


def main(data_path, db_path):
    genomes_frames = []
    for filename in os.listdir(data_path):
        print("collecting {} to database".format(filename))
        if filename.endswith(".tsv"):
            file_path = os.path.join(data_path, filename)
            single_genome_df = pd.read_csv(file_path, skiprows=0, sep="\t")
            single_genome_df.columns = [col_name.split(".")[-1] for col_name in single_genome_df.columns]
            single_genome_df = single_genome_df[HARDCODED_COLUMNS_ORDER]
            single_genome_df["genome"] = os.path.splitext(filename)[0]
            genomes_frames.append(single_genome_df)
    conn = sqlite3.connect(db_path)
    summary_df = pd.concat(genomes_frames, ignore_index=True)
    summary_df.to_sql("genes", conn, if_exists="replace")


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-i")
    parser.add_argument("-o")
    args = parser.parse_args()
    print("input dir {} output file {}".format(args.i, args.o))
    main(args.i, args.o)
