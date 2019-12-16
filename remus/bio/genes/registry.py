import os
import re
import sqlite3
import logging

import pandas as pd

from remus.bio.bed.beds_loading import BedLoader


def convert_genome_build(genome, hg19_expected="hg19", hg38_expected="hg38"):
    if re.match("(hg37|hg19|b37)", genome, re.IGNORECASE):
        return hg19_expected
    elif re.match("(hg38|grch38|b38)", genome, re.IGNORECASE):
        return hg38_expected
    raise InvalidGenomeBuildException(genome)
    

class GenesDBRegistry:
    query_genomes = """SELECT DISTINCT genome FROM genes"""
    query_gene_symbols = """SELECT DISTINCT geneSymbol FROM genes 
                     WHERE genome =='{genome}' AND geneSymbol LIKE '{pattern}%' 
                     ORDER BY geneSymbol
                     LIMIT '{limit_to}'"""
    query_genes = """SELECT DISTINCT * FROM genes 
                            WHERE genome=='{genome}' AND geneSymbol IN ({genes})
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
        query = self.query_gene_symbols.format(genome=genome, pattern=pattern, limit_to=limit)
        matching_genes_df = self._query_db(query)
        return matching_genes_df["geneSymbol"].tolist()

    def get_bed(self, genome, gene_names):
        genome = convert_genome_build(genome)
        self.logger.info("Querying genome %s for genes: %s" % (genome, gene_names))
        genes_list = ','.join(["'"+g+"'" for g in gene_names])
        query = self.query_genes.format(genome=genome, genes=genes_list)
        genes_df = self._query_db(query)
        coordinates = self._extract_gene_coordinates(genes_df)
        self.logger.info("Returned %s records" % (len(coordinates)))
        loader = BedLoader(src="\n".join(coordinates), from_string=True)
        return loader.bed

    def _query_db(self, query):
        with sqlite3.connect(self.db_url) as conn:
            df = pd.read_sql_query(query, conn)
            return df

    def _extract_gene_coordinates(self, genes_df):
        #'chrom', 'txStart', 'txEnd', 'strand', 'cdsStart', 'cdsEnd', 'exonCount', 'exonStarts', 'exonEnds', 'name', 'geneSymbol'
        sources_df = genes_df.iloc[:, [1, 2, 3, 11, 10, 7, 4]]
        strings_df = sources_df.apply(self.convert_to_bed_record, axis=1)
        return strings_df.tolist()

    @staticmethod
    def convert_to_bed_record(row):
        return '\t'.join([str(e) for e in row[:3]] +
                         [str(row[3]) + "(" + str(row[4]) + ")"] +
                         [str(e) for e in row[5:]])


    
class InvalidGenomeBuildException(Exception):
    pass

