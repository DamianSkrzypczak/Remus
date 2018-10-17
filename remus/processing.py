import logging
import re
from collections import OrderedDict

from flask import g

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

        genome = convert_genome_name(genome, desirable_older_format="hg37")
        
        BedsProcessor.logger().info("Querying gene database for %s" % genes)
        
        gene_beds = [g.genes_registry.get_bed(genome, gene) for gene in genes]
        
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
        beds = BedsProcessor._get_tss_fantom5_beds(tissues)  
        
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
        beds = BedsProcessor._get_enhancers_fantom5_beds(tissues)
        
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
        beds = BedsProcessor._get_enhancers_encode_beds(tissues)
        
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
        beds = BedsProcessor._get_accessible_chromatin_encode_beds(tissues)
        
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
        
        mirna_symbols = BedsProcessor._get_mirnas_targetting_genes(genes, g.mirna_target_registries['mirtarbase'], include_weak_support=include_weak_support)
        BedsProcessor.log_count("miRNA symbol list", mirna_symbols)

        accessible_mirna = BedsProcessor._get_accessible_mirnas(mirna_symbols, tissues, genome, combine_mode)
        BedsProcessor.log_count("Accessible miRNA", mirna_symbols)
        BedsProcessor.log_bed(accessible_mirna)

        return accessible_mirna

    @staticmethod
    def get_mirnas_targetting_genes_from_mirwalk(genes, tissues, genome, combine_mode, min_confidence, *args):
        """ """
        BedsProcessor.logger().info("Extracting miRWalk miRNAs targetting genes (%s): %s ; in tissues: %s ; and combine_mode: %s" % (genome, genes, tissues, combine_mode))

        mirna_symbols = BedsProcessor._get_mirnas_targetting_genes(genes, g.mirna_target_registries['mirwalk'], min_confidence=min_confidence)
        BedsProcessor.log_count("miRNA symbol list", mirna_symbols)
        
        accessible_mirna = BedsProcessor._get_accessible_mirnas(mirna_symbols, tissues, genome, combine_mode)
        BedsProcessor.log_count("Accessible miRNA", mirna_symbols)
        BedsProcessor.log_bed(accessible_mirna)

        return accessible_mirna

  

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
        #accessible_chromatin = BedsProcessor._get_accessible_chromatin_encode_beds(tissues)
        #accessible_chromatin_aggregate = BedsProcessor._combine_beds(accessible_chromatin, combine_mode)
        #accessible_mirna = BedsProcessor._combine_beds([mirna_bed] + accessible_chromatin_aggregate, combine_mode)
        #print(accessible_mirna)
        accessible_mirna=mirna_bed
        
        return accessible_mirna
        
     
    def _combine_beds(beds, combine_mode, merge=False):
        if combine_mode == "all":
            return BedOperations.intersect(beds, merge=merge).result
        elif combine_mode == "any":
            return BedOperations.union(beds, merge=merge).result
        else:
            return []

    def _get_gene_promoter_sites(genes, genome, upstream, downstream):
        genome = convert_genome_name(genome, desirable_older_format="hg19")
        genes_bed = BedsProcessor.get_genes_bed(genes, genome)[0]
        promoters = BedOperations.get_promoter_region(genes_bed, upstream, downstream, genome)
        return promoters.result

    def _get_enhancers_fantom5_beds(tissues):
        results = [g.tissues_registry.get_bed(tissue, "ENH_F5") for tissue in tissues]
        return [i for i in results if i]

    def _get_tss_fantom5_beds(tissues):
        results = [g.tissues_registry.get_bed(tissue, "TSS_F5") for tissue in tissues]
        return [i for i in results if i]
        
        
    def _get_enhancers_encode_beds(tissues):
        results = [g.tissues_registry.get_bed(tissue, "ENH_EN") for tissue in tissues]
        return [i for i in results if i]

    def _get_accessible_chromatin_encode_beds(tissues):
        results = [g.tissues_registry.get_bed(tissue, "CHRM") for tissue in tissues]
        return [i for i in results if i]


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
                 BedsProcessor.get_mirnas_targetting_genes_from_mirtarbase)
             ),
            ("mirna-targets-mirwalk",
             self._get_bed_files(
                 self.mirna_target_mirwalk_params,
                 BedsProcessor.get_mirnas_targetting_genes_from_mirwalk)
             )
        ])
        return bed_files

    def _get_bed_files(self, params, getter_method):
        params_values = [self._data.get(p) for p in params]
        if all(params_values):
            self._logger.info("All values provided, running {}".format(getter_method.__name__))
            return getter_method(*params_values)
        else:
            self._logger.info("NOT all values provided for {} => values:{}".format(getter_method.__name__, params_values))
            return []


def get_matching_genes(pattern, genome_name, limit):
    genome_name = convert_genome_name(genome_name)
    if pattern and genome_name and (genome_name in g.genes_registry.available_genomes):
        return g.genes_registry.get_matching_genes(genome_name, pattern, limit)
    else:
        return []


def get_matching_tissues(pattern, limit):
    return g.tissues_registry.get_matching_tissues(pattern, limit)


def convert_genome_name(genome, desirable_older_format="hg37"):
    if re.match("(hg37|hg19)", genome, re.IGNORECASE):
        return desirable_older_format
    else:
        return genome.lower()

