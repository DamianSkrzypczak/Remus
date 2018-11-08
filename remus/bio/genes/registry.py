import os
import re
import sqlite3
import logging

import pandas as pd

from remus.bio.bed.beds_loading import BedLoader


def convert_genome_build(genome, hg19_expected="hg37", hg38_expected="hg38"):
    if re.match("(hg37|hg19|b37)", genome, re.IGNORECASE):
        return hg19_expected
    elif re.match("(hg38|grch38|b38)", genome, re.IGNORECASE):
        return hg38_expected
    raise InvalidGenomeBuildException(genome)
    

class GenesDBRegistry:
    query_genomes = """SELECT DISTINCT genome FROM genes"""
    query_genes = """SELECT DISTINCT geneSymbol FROM genes 
                     WHERE genome =='{genome}' AND geneSymbol LIKE '{pattern}%' 
                     ORDER BY geneSymbol
                     LIMIT '{limit_to}'"""
    query_gene_sources = """SELECT DISTINCT * FROM genes 
                            WHERE genome=='{genome}' AND geneSymbol=='{gene}'
                            ORDER BY chrom,txStart,txEnd"""

    instance = None

    @staticmethod
    def get_instance():
        if not GenesDBRegistry.instance:
            GenesDBRegistry.instance = GenesDBRegistry()
        return GenesDBRegistry.instance

    def __init__(self, genes_db=os.path.join("data", "genes", "genes.db")):
        self.logger = logging.getLogger(self.__class__.__name__)
        self.logger.info("Created GenesDBRegistry")
        self.db_url = genes_db

    @property
    def available_genomes(self):
        genomes_df = self._query_db(self.query_genomes)
        return genomes_df["genome"].tolist()

    def get_matching_genes(self, genome, pattern, limit):
        genome = convert_genome_build(genome)
        query = self.query_genes.format(genome=genome, pattern=pattern, limit_to=limit)
        matching_genes_df = self._query_db(query)
        return matching_genes_df["geneSymbol"].tolist()

    def get_bed(self, genome, gene_name):
        genome = convert_genome_build(genome)
        self.logger.info("Querring genome %s for gene %s" % (genome, gene_name))
        query = self.query_gene_sources.format(genome=genome, gene=gene_name)
        genes_df = self._query_db(query)
        sources = self._extract_sources(genes_df)
        self.logger.info("Returned %s records" % len(sources))
        loader = BedLoader(src="\n".join(sources), from_string=True)
        return loader.bed

    def _query_db(self, query):
        with sqlite3.connect(self.db_url) as conn:
            df = pd.read_sql_query(query, conn)
            return df

    @staticmethod
    def _extract_sources(genes_df):
        sources_df = genes_df.iloc[:, [1, 2, 3, 10, 7, 4]]
        strings_df = sources_df.apply(lambda x: '\t'.join([str(i) for i in x]), axis=1)
        return strings_df.tolist()

    @staticmethod
    def _extract_names(genes_df):
        return genes_df.iloc[:, -1]


    
class InvalidGenomeBuildException(Exception):
    pass

