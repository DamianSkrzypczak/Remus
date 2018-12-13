import os
import sqlite3
import logging 
import pandas as pd


class GenericMiRNATargetRegistry:
    """ Interface for MiRNA target sources """
    
    def __init__(self, db_file):
        self.logger = logging.getLogger(self.__class__.__name__)
        self.logger.info("Created miRNA target registry")
        self.db_url = db_file


    
    def get_mirnas_targetting_gene(self, gene_symbol, **kwargs):
        """ Convenience method for queryying one gene """
        return self.get_mirnas_targetting_genes([gene_symbol], **kwargs)

    def get_mirnas_targetting_genes(self, gene_symbol_list, **kwargs):
        """ abstract method to fill in in subclasses """
        raise NotImplementedError("Please implement method [get_mirnas_targetting_genes] in class [%s]" % type(self).__name__)

   
    def _query_db(self, query):
        with sqlite3.connect(self.db_url) as conn:
            return pd.read_sql_query(query, conn)
 

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
            
        mir_symbol = 'MIR' + ('LET' if is_let else '') + suffix.upper()
        
        logging.getLogger("GenericMiRNATargetRegistry").info("Coverted %s to %s" % (mir, mir_symbol))
        
        return mir_symbol



class MirTarBaseRegistry(GenericMiRNATargetRegistry):
    

    def __init__(self, db_file = os.path.join("data", "mirna", "targets.db")):
        super(MirTarBaseRegistry, self).__init__(db_file)
        self.query = "SELECT DISTINCT mirna FROM mirtarbase " + \
                     "WHERE target_gene IN ({genes}) AND " + \
                     "(support_type=='Functional MTI' OR support_type=='{other_support_type}')"

    
    def get_mirnas_targetting_genes(self, gene_symbol_list, include_weak_support=False):

        self.logger.info("Querrying registry for [%s] (weak_support: %s)" % (gene_symbol_list, include_weak_support))

        other_support_type = 'Functional MTI (Weak)' if include_weak_support else 'Functional MTI'
        genes = ','.join(["'"+g+"'" for g in gene_symbol_list])
        query = self.query.format(genes=genes, other_support_type=other_support_type)
        
        mirnas_df = self._query_db(query)
        mirna_list = mirnas_df['mirna'].values.tolist()
                
        self.logger.info("Returning %s mirna ids: [%s]..." % (len(mirna_list), mirna_list[:3]))
        
        return mirna_list



        
class MirWalkRegistry(GenericMiRNATargetRegistry):
    
    def __init__(self, db_file = os.path.join("data", "mirna", "targets.db")):
        super(MirWalkRegistry, self).__init__(db_file)
        self.query = "SELECT DISTINCT mirna FROM mirwalk_3UTR " + \
                     "WHERE target_gene IN ({genes}) AND " + \
                     "(confidence>='{minimal_confidence}')"

    def get_mirnas_targetting_genes(self, gene_symbol_list, min_confidence=0.9):
        
        self.logger.info("Querrying registry for [%s] (min_confidence: %s)" % (gene_symbol_list, min_confidence))

        genes = ','.join(["'"+g+"'" for g in gene_symbol_list])
        query = self.query.format(genes=genes, minimal_confidence=min_confidence)        
        mirnas_df = self._query_db(query)
        mirna_list = mirnas_df['mirna'].values.tolist()
        
        self.logger.info("Returning %s mirna ids: [%s]..." % (len(mirna_list), mirna_list[:3]))

        return mirna_list

        
class MiRNATargetRegistryFactory:
    
    MIRTARBASE_KEY = 'mirtarbase'
    MIRWALK_KEY = 'mirwalk'
    
    AVAILABLE_REGISTRIES = {MIRTARBASE_KEY: MirTarBaseRegistry, 
                            MIRWALK_KEY: MirWalkRegistry}
                            
    instances = None
    
    @staticmethod
    def get_instance(name):
        if name not in MiRNATargetRegistryFactory.AVAILABLE_REGISTRIES:
            raise InvalidMiRNATargetRegistryException(name)
        if not MiRNATargetRegistryFactory.instances: 
            MiRNATargetRegistryFactory.instances = {}
        if name not in MiRNATargetRegistryFactory.instances:
            MiRNATargetRegistryFactory.instances[name] = MiRNATargetRegistryFactory.AVAILABLE_REGISTRIES[name]()
        return MiRNATargetRegistryFactory.instances[name]
    

class InvalidMiRNATargetRegistryException(Exception):
    pass
