import logging
from collections import OrderedDict

from remus.bio.regulatory_regions.registry import RegulatoryRegionsFilesRegistry
from remus.bio.genes.registry import GenesDBRegistry
from remus.bio.mirna.registry import MiRNATargetRegistryFactory

from remus.bio.bed.beds_operations import BedOperations, BedOperationResult
   

class BedsProcessor:

    def logger():
        return logging.getLogger("BedsProcessor")
    
    def log_count(msg, bed):
        BedsProcessor.logger().info(msg +" had [%s] records." % len(bed))
        #BedsProcessor.logger().debug(msg +" had [%s] records." % len(bed))
        
    def log_bed(bed):
        BedsProcessor.logger().debug(":\n%s" % str(bed))

    
    @staticmethod
    def get_genes_bed(genes, genome, *args):
        
        BedsProcessor.logger().info("Querying gene database for %s" % genes)
        
        registry = GenesDBRegistry.get_instance()
        gene_beds = [registry.get_bed(genome, gene) for gene in genes]
        
        BedsProcessor.logger().info("Got [%s] non-empty BED files" % len([b for b in gene_beds if b]))
        BedsProcessor.log_bed(gene_beds)
        
        result = BedOperations.union(gene_beds).result
        
        BedsProcessor.log_count("Union of the BED files",result)
        BedsProcessor.log_bed(result)
        
        return [ result ]


    @staticmethod
    def get_tss_fantom5_bed(genes, tissues, genome, combine_mode, upstream, downstream, *args):

        BedsProcessor.logger().info("Extracting F5 TSS for genes (%s): %s; tissues: %s; and combine_mode: %s" % (genome, genes, tissues, combine_mode))
                
        flanked_genes = BedsProcessor._get_gene_promoter_sites(genes, genome, upstream, downstream)
        beds = BedsProcessor._get_regulatory_regions_bed(genome, tissues, RegulatoryRegionsFilesRegistry.FANTOM5_TSS_KEY)
        
        BedsProcessor.log_count("Flanked genes' promoters BED", flanked_genes)
        BedsProcessor.log_bed(flanked_genes)
        BedsProcessor.logger().info("F5 TSS BEDs list:\n%s" % str(beds))
        
        if beds and flanked_genes:
            joined_f5_tss = BedsProcessor._combine_beds(beds, combine_mode)
            BedsProcessor.log_count("Combined F5 TSS BED", joined_f5_tss)
            BedsProcessor.log_bed(joined_f5_tss)
            
            result = BedOperations.intersect([joined_f5_tss, flanked_genes], **{"u": True}).result           
            BedsProcessor.log_count("F5 TSS BED intersected with genes' promoters",result)
            BedsProcessor.log_bed(result)
            
            return [result]
        else:
            BedsProcessor.logger().info("Returning empty TSS list") 
            return []


    @staticmethod
    def get_enhancers_fantom5_bed(genes, tissues, genome, combine_mode, upstream, downstream, *args):
        
        BedsProcessor.logger().info("Extracting F5 enhancers for genes (%s): %s; tissues: %s; and combine_mode: %s" % (genome, genes, tissues, combine_mode))

        flanked_genes = BedsProcessor._get_gene_promoter_sites(genes, genome,
                                                            int(float(upstream) * 1000),
                                                            int(float(downstream) * 1000))
        beds = BedsProcessor._get_regulatory_regions_bed(genome, tissues, RegulatoryRegionsFilesRegistry.FANTOM5_ENHANCERS_KEY)
        
        BedsProcessor.log_count("Flanked genes' promoters BED", flanked_genes)
        BedsProcessor.log_bed(flanked_genes)
        BedsProcessor.logger().info("F5 enhancer BEDs list:\n%s" % str(beds))        
        
        if beds and flanked_genes:
            joined_f5_enh = BedsProcessor._combine_beds(beds, combine_mode)
            BedsProcessor.log_count("Combined F5 enhancers BED", joined_f5_enh)
            BedsProcessor.log_bed(joined_f5_enh)
            
            result = BedOperations.intersect([joined_f5_enh, flanked_genes], **{"u": True}).result
            BedsProcessor.log_count("F5 enhancers BED intersected with genes' promoters", result)
            BedsProcessor.log_bed(result)

            return [ result ]
        else:
            BedsProcessor.logger().info("Returning empty enhancers list") 
            return []

    @staticmethod
    def get_enhancers_encode_bed(genes, tissues, genome, combine_mode, upstream, downstream, *args):
        
        BedsProcessor.logger().info("Extracting ENCODE enhancers for genes (%s): %s; tissues: %s; and combine_mode: %s" % (genome, genes, tissues, combine_mode))

        flanked_genes = BedsProcessor._get_gene_promoter_sites(genes, genome,
                                                            int(float(upstream) * 1000),
                                                            int(float(downstream) * 1000))
        beds = BedsProcessor._get_regulatory_regions_bed(genome, tissues, RegulatoryRegionsFilesRegistry.ENCODE_ENHANCERS_KEY)
        
        BedsProcessor.log_count("Flanked genes' promoters BED", flanked_genes)
        BedsProcessor.log_bed(flanked_genes)
        BedsProcessor.logger().info("ENCODE enhancer BEDs list:\n%s" % str(beds))        

        if beds and flanked_genes:
            joined_enc_enh = BedsProcessor._combine_beds(beds, combine_mode)
            BedsProcessor.log_count("Combined ENCODE enhancers BED", joined_enc_enh)
            BedsProcessor.log_bed(joined_enc_enh)

            result = BedOperations.intersect([joined_enc_enh, flanked_genes], **{"u": True}).result
            BedsProcessor.log_count("ENCODE enhancers BED intersected with genes' promoters", result)
            BedsProcessor.log_bed(result)

            return [ result ]
        else:
            BedsProcessor.logger().info("Returning empty enhancers list") 
            return []

    @staticmethod
    def get_accessible_chromatin_bed(genes, tissues, genome, combine_mode, upstream, downstream, *args):
        
        BedsProcessor.logger().info("Extracting accessible chromatin for genes (%s): %s; tissues: %s; and combine_mode: %s" % (genome, genes, tissues, combine_mode))

        flanked_genes = BedsProcessor._get_gene_promoter_sites(genes, genome,
                                                            int(float(upstream) * 1000),
                                                            int(float(downstream) * 1000))
        beds = BedsProcessor._get_regulatory_regions_bed(genome, tissues, RegulatoryRegionsFilesRegistry.ENCODE_CHROMATIN_KEY)
        
        BedsProcessor.log_count("Flanked genes' promoters BED", flanked_genes)
        BedsProcessor.log_bed(flanked_genes)
        BedsProcessor.logger().info("Accessible chromatin BEDs list:\n%s" % str(beds))        
        
        if beds and flanked_genes:
            combined_chromatin = BedsProcessor._combine_beds(beds, combine_mode)
            BedsProcessor.log_count("Combined accessible chromatin BED", combined_chromatin)
            BedsProcessor.log_bed(combined_chromatin)

            result = BedOperations.intersect([combined_chromatin, flanked_genes], **{"u": True}).result
            BedsProcessor.log_count("Accessible chromatin BED intersected with genes' promoters", result)
            BedsProcessor.log_bed(result)

            return [ result ]
        else:
            BedsProcessor.logger().info("Returning empty open chromatin list") 
            return []

    
    @staticmethod
    def get_mirnas_targetting_genes_from_mirtarbase(genes, tissues, genome, combine_mode, include_weak_support, *args):
        """ """
        BedsProcessor.logger().info("Extracting miRTarBase miRNAs targetting genes (%s): %s ; in tissues: %s ; and combine_mode: %s" % (genome, genes, tissues, combine_mode))
        
        mirtarbase_registry = MiRNATargetRegistryFactory.get_instance(MiRNATargetRegistryFactory.MIRTARBASE_KEY)
        mirna_symbols = BedsProcessor._get_mirnas_targetting_genes(genes, mirtarbase_registry, include_weak_support=True if include_weak_support else False)
        BedsProcessor.log_count("miRNA symbol list", mirna_symbols)

        accessible_mirna = BedsProcessor._get_accessible_mirnas(mirna_symbols, tissues, genome, combine_mode)
        BedsProcessor.log_count("Accessible miRNA", accessible_mirna)
        BedsProcessor.log_bed(accessible_mirna)

        return [ accessible_mirna ]

    @staticmethod
    def get_mirnas_targetting_genes_from_mirwalk(genes, tissues, genome, combine_mode, min_confidence, *args):
        """ """
        BedsProcessor.logger().info("Extracting miRWalk miRNAs targetting genes (%s): %s ; in tissues: %s ; and combine_mode: %s" % (genome, genes, tissues, combine_mode))

        mirwalk_registry = MiRNATargetRegistryFactory.get_instance(MiRNATargetRegistryFactory.MIRWALK_KEY)
        mirna_symbols = BedsProcessor._get_mirnas_targetting_genes(genes, mirwalk_registry, min_confidence=min_confidence)
        BedsProcessor.log_count("miRNA symbol list", mirna_symbols)
        
        accessible_mirna = BedsProcessor._get_accessible_mirnas(mirna_symbols, tissues, genome, combine_mode)
        BedsProcessor.log_count("Accessible miRNA", accessible_mirna)
        BedsProcessor.log_bed(accessible_mirna)

        return [ accessible_mirna ]

  

    #
    # private methods below
    #

    def _get_mirnas_targetting_genes(genes, registry, **args):
        mirs = set()
        for gene in genes:
            mirs.update(registry.get_mirnas_targetting_gene(gene, **args))
        return registry.get_mirna_gene_symbols(list(mirs))

    def _get_accessible_mirnas(mirna_symbols, tissues, genome, combine_mode):
        
        mirna_bed = BedsProcessor.get_genes_bed(mirna_symbols, genome)
                
        # intersect beds with accessible chromatin in tissues
        accessible_chromatin = BedsProcessor._get_regulatory_regions_bed(genome, tissues, RegulatoryRegionsFilesRegistry.ENCODE_CHROMATIN_KEY)
        accessible_chromatin_aggregate = BedsProcessor._combine_beds(accessible_chromatin, combine_mode)
        print(mirna_bed)
        accessible_mirna = BedOperations.intersect(mirna_bed + [accessible_chromatin_aggregate], merge=False).result
        #print(accessible_mirna)
        #accessible_mirna=mirna_bed
        
        return accessible_mirna
        
     
    def _combine_beds(beds, combine_mode, merge=False):
        if combine_mode == "all":
            return BedOperations.intersect(beds, merge=merge).result
        elif combine_mode == "any":
            return BedOperations.union(beds, merge=merge).result
        else: 
            raise InvalidCombineOperationException

    def _get_gene_promoter_sites(genes, genome, upstream, downstream):
        genes_bed = BedsProcessor.get_genes_bed(genes, genome)[0]
        promoters = BedOperations.get_promoter_region(genes_bed, upstream, downstream)
        return promoters.result

    def _get_regulatory_regions_bed(genome, tissues, reg_feature_type):
        registry = RegulatoryRegionsFilesRegistry.get_registry(genome)
        results = [registry.get_bed(tissue, reg_feature_type) for tissue in tissues]
        return [i for i in results if i]



class InvalidCombineOperationException(Exception):
    pass


class BedsCollector:
    genes_params = ["genes", "genome"]

    tss_fantom5_params = [
        "genes", "tissues", "genome",
        "transcription-fantom5-combine-mode",
        "transcription-fantom5-kbs-upstream",
        "transcription-fantom5-kbs-downstream",
        "transcription-fantom5-used"
    ]

    enhancers_fantom5_params = [
        "genes", "tissues", "genome",
        "enhancers-fantom5-combine-mode",
        "enhancers-fantom5-kbs-upstream",
        "enhancers-fantom5-kbs-downstream",
        "enhancers-fantom5-used"
    ]

    enhancers_encode_params = [
        "genes", "tissues", "genome",
        "enhancers-encode-combine-mode",
        "enhancers-encode-kbs-upstream",
        "enhancers-encode-kbs-downstream",
        "enhancers-encode-used"
    ]

    accessible_chromatin_encode_params = [
        "genes", "tissues", "genome",
        "accessible-chromatin-encode-combine-mode",
        "accessible-chromatin-encode-kbs-upstream",
        "accessible-chromatin-encode-kbs-downstream",
        "accessible-chromatin-encode-used"
    ]
    
    mirna_target_mirtarbase_params = [
        "genes", "tissues", "genome",
        "mirna-targets-combine-mode", 
        "mirna-mirtarbase-include-weak",
        "mirna-mirtarbase-used"
    ]

    mirna_target_mirwalk_params = [
        "genes", "tissues", "genome",
        "mirna-targets-combine-mode", 
        "mirna-mirwalk-minimal-confidence",
        "mirna-mirwalk-used"
    ]

    def __init__(self, data):
        self._data = data
        self._logger = logging.getLogger(self.__class__.__name__)

    def collect_bed_files(self):
        bed_files = OrderedDict([
            ("genes",
             self._get_bed_files(
                 self.genes_params,
                 BedsProcessor.get_genes_bed)
             ),
            ("transcription-fantom5",
             self._get_bed_files(
                 self.tss_fantom5_params,
                 BedsProcessor.get_tss_fantom5_bed)
             ),
            ("enhancers-fantom5",
             self._get_bed_files(
                 self.enhancers_fantom5_params,
                 BedsProcessor.get_enhancers_fantom5_bed)
             ),
            ("enhancers-encode",
             self._get_bed_files(
                 self.enhancers_encode_params,
                 BedsProcessor.get_enhancers_encode_bed)
             ),
            ("accessible-chromatin-fantom5",
             self._get_bed_files(
                 self.accessible_chromatin_encode_params,
                 BedsProcessor.get_accessible_chromatin_bed)
             ),
            ("mirna-targets-mirtarbase",
             self._get_bed_files(
                 self.mirna_target_mirtarbase_params,
                 BedsProcessor.get_mirnas_targetting_genes_from_mirtarbase,
                 [e for e in self.mirna_target_mirtarbase_params if e != "mirna-mirtarbase-include-weak"])
             ),
            ("mirna-targets-mirwalk",
             self._get_bed_files(
                 self.mirna_target_mirwalk_params,
                 BedsProcessor.get_mirnas_targetting_genes_from_mirwalk)
             )
        ])
        return bed_files

    def _get_bed_files(self, params, getter_method, required_params=None):
        
        params_values = [self._data.get(p) for p in params]
        required_params_values = [self._data.get(p) for p in required_params] if required_params else params_values
        if all(required_params_values):
            self._logger.info("All required values provided: {}. Running {}".format(params_values, getter_method.__name__))
            return getter_method(*params_values)
        else:
            self._logger.info("NOT all required ({}) values provided for {} => values:{}".format(required_params, getter_method.__name__, params_values))
            return []
