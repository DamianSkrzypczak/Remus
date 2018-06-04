# import os
import os
import sqlite3

import pandas as pd

db_path = 'genes.db'

data_path = "genomes"
genomes_frames = []
for filename in os.listdir(data_path):
    if filename.endswith(".tsv"):
        file_path = os.path.join(data_path, filename)
        single_genome_df = pd.read_csv(file_path, skiprows=0, sep="\t")
        single_genome_df.columns = [col_name.split(".")[-1] for col_name in single_genome_df.columns]
        single_genome_df["genome"] = os.path.splitext(filename)[0]
        genomes_frames.append(single_genome_df)
conn = sqlite3.connect(db_path)
summary_df = pd.concat(genomes_frames, ignore_index=True)
summary_df.to_sql("genes", conn, if_exists="replace")

