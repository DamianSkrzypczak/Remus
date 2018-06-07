# import os
import os
import sqlite3

import pandas as pd

DB_PATH = 'genes.db'
DATA_PATH = "genes_db_sources"
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
GENOMES_FRAMES = []
for filename in os.listdir(DATA_PATH):
    if filename.endswith(".tsv"):
        file_path = os.path.join(DATA_PATH, filename)
        single_genome_df = pd.read_csv(file_path, skiprows=0, sep="\t")
        single_genome_df.columns = [col_name.split(".")[-1] for col_name in single_genome_df.columns]
        single_genome_df = single_genome_df[HARDCODED_COLUMNS_ORDER]
        single_genome_df["genome"] = os.path.splitext(filename)[0]
        GENOMES_FRAMES.append(single_genome_df)
conn = sqlite3.connect(DB_PATH)
summary_df = pd.concat(GENOMES_FRAMES, ignore_index=True)
summary_df.to_sql("genes", conn, if_exists="replace")
