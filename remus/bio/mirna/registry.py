import os
import sqlite3

import pandas as pd


class GenericMiRNATargetRegistry:
    """ Interface for MiRNA target sources """
    
    def get_mirnas_targetting_gene(self, gene_symbol, **args):
        """ abstract method to fill in in subclasses """
        raise NotImplementedError("Please implement method [get_mirnas_targetting_gene] in class [%s]" % type(self).__name__)

    @staticmethod
    def get_mirna_gene_symbols(mirnas):
        return [ GenericMiRNATargetRegistry._get_mirna_gene_symbol(mir) for mir in mirnas ]
        
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

        

class MirTarBaseRegistry(GenericMiRNATargetRegistry):
    

    def __init__(self, db_file = os.path.join("data", "mirna", "targets.db")):
        self.conn = sqlite3.connect(db_file)
        self.query = "SELECT mirna FROM mirtarbase " + \
                     "WHERE target_gene=='{gene}' AND " + \
                     "(support_type=='Functional MTI' OR support_type=='{other_support_type}')"
        
            
            
    def get_mirnas_targetting_gene(self, gene_symbol, include_weak_support=False):
        other_support_type = 'Functional MTI (Weak)' if include_weak_support else 'Functional MTI'
        query = self.query.format(gene=gene_symbol, other_support_type=other_support_type)
        
        mirnas_df = pd.read_sql_query(query, self.conn)
        
        return mirnas_df['mirna'].values.tolist()
        
        
    def teardown_registry(self):
        self.conn.close()

        
        
        
