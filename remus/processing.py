import logging
import re
from collections import OrderedDict

from flask import g

from remus.bio.bed.beds_operations import BedOperations, BedOperationResult


class BedsProcessor:
    
    @staticmethod
    def get_genes_bed(genes, genome, *args):
        genome = convert_genome_name(genome, desirable_older_format="hg37")
        gene_beds = [g.genes_registry.get_bed(genome, gene) for gene in genes]
        return [ BedOperations.union(gene_beds).result ]

    @staticmethod
    def get_tss_fantom5_bed(genes, tissues, genome, combine_mode, upstream, downstream, *args):
        flanked_genes = BedsProcessor._get_gene_promoter_sites(genes, genome, upstream, downstream)
        
        beds = BedsProcessor._get_tss_fantom5_beds(tissues)  
        
        if beds and flanked_genes:
            
            joined_f5_tss = BedsProcessor._combine_beds(beds, combine_mode)
            return [ BedOperations.intersect([joined_f5_tss, flanked_genes], **{"u": True}).result ]
        else:
            return []

    @staticmethod
    def get_enhancers_fantom5_bed(genes, tissues, genome, combine_mode, upstream, downstream, *args):
        flanked_genes = BedsProcessor._get_gene_promoter_sites(genes, genome,
                                                            int(float(upstream) * 1000),
                                                            int(float(downstream) * 1000))
        beds = BedsProcessor._get_enhancers_fantom5_beds(tissues)
        if beds and flanked_genes:
            joined_f5_enh_tissues = BedsProcessor._combine_beds(beds, combine_mode)
            return [ BedOperations.intersect([joined_f5_enh_tissues, flanked_genes], **{"u": True}).result ]
        else:
            return []

    @staticmethod
    def get_enhancers_encode_bed(genes, tissues, genome, combine_mode, upstream, downstream, *args):
        flanked_genes = BedsProcessor._get_gene_promoter_sites(genes, genome,
                                                            int(float(upstream) * 1000),
                                                            int(float(downstream) * 1000))
        beds = BedsProcessor._get_enhancers_encode_beds(tissues)
        if beds and flanked_genes:
            joined_f5_enh_tissues = BedsProcessor._combine_beds(beds, combine_mode)
            return [ BedOperations.intersect([joined_f5_enh_tissues, flanked_genes], **{"u": True}).result ]
        else:
            return []

    @staticmethod
    def get_accessible_chromatin_bed(genes, tissues, genome, combine_mode, upstream, downstream, *args):
        flanked_genes = BedsProcessor._get_gene_promoter_sites(genes, genome,
                                                            int(float(upstream) * 1000),
                                                            int(float(downstream) * 1000))
        beds = BedsProcessor._get_accessible_chromatin_encode_beds(tissues)
        if beds and flanked_genes:
            joined_f5_enh_tissues = BedsProcessor._combine_beds(beds, combine_mode)
            return [ BedOperations.intersect([joined_f5_enh_tissues, flanked_genes], **{"u": True}).result ]
        else:
            return []

    @staticmethod
    def get_mirnas_targetting_genes(genes, tissues, genome, combine_mode, include_weak_interactions, *args):
        mirna_symbols = BedsProcessor._get_mirnas_targetting_genes(genes)
        mirna_bed = BedsProcessor.get_genes_bed(mirna_symbols, genome)
        
        print(args)
        print(mirna_bed)
        
        # intersect beds with accessible chromatin in tissues
        accessible_chromatin = BedsProcessor._get_accessible_chromatin_encode_beds(tissues)
        accessible_chromatin_aggregate = BedsProcessor._combine_beds(accessible_chromatin, combine_mode)
        accessible_mirna = BedsProcessor._combine_beds([mirna_bed] + accessible_chromatin_aggregate, combine_mode)
        
        print(accessible_mirna)
        
        return accessible_mirna
        
        #return mirna_bed

    #
    # private methods below
    #
    
    @staticmethod
    def _get_mirnas_targetting_genes(genes):
        mirs = []
        for gene in genes:
            mirs += g.mirna_target_registry.get_mirnas_targetting_gene(gene)
        
        return g.mirna_target_registry.get_mirna_gene_symbols(list(set(mirs)))

    
    @staticmethod
    def _combine_beds(beds, combine_mode, merge=False):
        if combine_mode == "all":
            return BedOperations.intersect(beds, merge=merge).result
        elif combine_mode == "any":
            return BedOperations.union(beds, merge=merge).result
        else:
            return []

    @staticmethod
    def _get_gene_promoter_sites(genes, genome, upstream, downstream):
        genome = convert_genome_name(genome, desirable_older_format="hg19")
        genes_bed = BedsProcessor.get_genes_bed(genes, genome)[0]
        promoters = BedOperations.get_promoter_region(genes_bed, upstream, downstream, genome)
        return promoters.result

    @staticmethod
    def _get_enhancers_fantom5_beds(tissues):
        results = [g.tissues_registry.get_bed(tissue, "ENH_F5") for tissue in tissues]
        return [i for i in results if i]

    @staticmethod
    def _get_tss_fantom5_beds(tissues):
        results = [g.tissues_registry.get_bed(tissue, "TSS_F5") for tissue in tissues]
        return [i for i in results if i]
        
        
    @staticmethod
    def _get_enhancers_encode_beds(tissues):
        results = [g.tissues_registry.get_bed(tissue, "ENH_EN") for tissue in tissues]
        return [i for i in results if i]

    @staticmethod
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
        "mirna-mirwalk-minimal-confidence"
        "mirna-mirwalk-used"
    ]

    def __init__(self, data):
        self._data = data

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
                 BedsProcessor.get_mirnas_targetting_genes)
             )
        ])
        return bed_files

    def _get_bed_files(self, params, getter_method):
        params_values = [self._data.get(p) for p in params]
        if all(params_values):
            logging.info("All values provided, running {}".format(getter_method.__name__))
            return getter_method(*params_values)
        else:
            logging.error("NOT all values provided for {} => values:{}".format(getter_method.__name__, params_values))
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

