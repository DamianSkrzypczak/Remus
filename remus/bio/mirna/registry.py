import os
import sqlite3

import pandas as pd

class MiRNATargetRegistry:
    
    #OPERATION_MAP = {"mirtarbase": MiRNATargetRegistry._get_mirnas_targetting_gene_from_mirtarbase}
    
    QUERY_MAP = {"mirtarbase": "SELECT mirna FROM MirTarBaseInteractions WHERE target_gene=='{gene}' AND (support_type=='Functional MTI' OR support_type=='{other_support_type}')"}


    def __init__(self, db_filedb=os.path.join("data", "mirna", "targets.db")):
        self.conn = sqlite3.connect(db_file)
    
    @staticmethod
    def get_mirna_gene_symbols(mirnas):
        return [ MiRNATargetRegistry._get_mirna_gene_symbol(mir) for mir in mirnas ]
        
    @staticmethod
    def _get_mirna_gene_symbol(mir):
        """ hsa-([miR|let])-(\d)+(-[35]p)?   --->  MIR(\d)+ """
        
        is_let = mir.find('hsa-let')==0
        
        suffix = mir[len('hsa-mir-'):]
        if suffix[-1] == '*': 
            suffix = suffix[:-1]
            
        
        if suffix[-3:]=='-3p' or suffix[-3:]=='-5p':
            suffix = suffix[:-3]
        
        suffix = ''.join(suffix.split('-'))
            
        return 'MIR' + ('LET' if is_let else '') + suffix.upper()
    
    def get_mirnas_targetting_gene(self, gene_symbol, source, **args):
        """ shared endpoint for querying miRNAs targetting a gene """
        
        if source not in self.OPERATION_MAP: 
            raise InvalidMiRNASourceException()
        
        operation = self.OPERATION_MAP[source]
        return self.operation(gene_symbol, **args)
        
            
    def _get_mirnas_targetting_gene_from_mirtarbase(self, gene_symbol, include_weak_support=False):
        other_support_type = 'Functional MTI (Weak)' if include_weak_support else 'Functional MTI'
        query = self.QUERY_MAP['mirtarbase'].format(gene=gene_symbol, other_support_type=other_support_type)
        
        mirnas_df = pd.read_sql_query(query, self.conn)
        
        return mirnas_df.tolist()
        
        
    def teardown_registry(self):
        self.conn.close()

        
        
class InvalidMiRNASourceException(Exception):
    pass
        
        
