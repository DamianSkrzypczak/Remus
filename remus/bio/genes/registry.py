import os
import sqlite3

import pandas as pd

from remus.bio.bed.beds_loading import BedLoader


class GenesDBRegistry:
    query_genomes = """SELECT DISTINCT genome FROM genes"""
    query_genes = """SELECT DISTINCT geneSymbol FROM genes 
                     WHERE genome =='{genome}' AND geneSymbol LIKE '{pattern}%' 
                     ORDER BY geneSymbol
                     LIMIT '{limit_to}'"""
    query_gene_sources = """SELECT DISTINCT * FROM genes 
                            WHERE genome=='{genome}' AND geneSymbol=='{gene}'
                            ORDER BY chrom,txStart,txEnd"""

    def __init__(self, genes_db=os.path.join("data", "genes", "genes.db")):
        self.conn = sqlite3.connect(genes_db)

    @property
    def available_genomes(self):
        query = self.query_genomes
        genomes_df = pd.read_sql_query(query, self.conn)
        return genomes_df["genome"].tolist()

    def get_matching_genes(self, genome, pattern, limit):
        query = self.query_genes.format(genome=genome, pattern=pattern, limit_to=limit)
        matching_genes_df = pd.read_sql_query(query, self.conn)
        return matching_genes_df["geneSymbol"].tolist()

    def get_bed(self, genome, gene_name):
        query = self.query_gene_sources.format(genome=genome, gene=gene_name)
        genes_df = pd.read_sql_query(query, self.conn)
        sources = self._extract_sources(genes_df)
        loader = BedLoader(src="\n".join(sources), from_string=True)
        return loader.bed

    @staticmethod
    def _extract_sources(genes_df):
        sources_df = genes_df.iloc[:, [1, 2, 3, 10, 7, 4]]
        strings_df = sources_df.apply(lambda x: '\t'.join([str(i) for i in x]), axis=1)
        return strings_df.tolist()

    @staticmethod
    def _extract_names(genes_df):
        return genes_df.iloc[:, -1]

    def teardown_registry(self):
        self.conn.close()
