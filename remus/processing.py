import re
from collections import OrderedDict

from flask import g

from remus.bio.bed.beds_operations import BedsMutualOperation, BedsFlanker


class BedGetters:
    @staticmethod
    def get_genes_beds(genes, genome):
        genome = convert_genome_name(genome, desirable_older_format="hg37")
        genes_beds = [g.genes_registry.get_bed(genome, gene) for gene in genes]
        return [BedsMutualOperation(genes_beds, operation="union").result]

    @staticmethod
    def get_transcription_starting_sites_fantom5_beds(genes, genome, flank_range, upstream, downstream):
        tss_genes = BedGetters._get_joined_flanked_genes(genes, genome, upstream, downstream)
        promoters = g.tss_registry.get_bed()
        return [
            BedsMutualOperation([tss_genes, promoters], operation="intersection", **{"wa": True}).result]

    @staticmethod
    def get_enhancers_fantom5_beds(genes, tissues, genome, flank_range, upstream, downstream):
        tss_genes = BedGetters._get_joined_flanked_genes(genes, genome, upstream, downstream)
        f5_enh_tissues = BedGetters._get_tissues_beds(tissues)
        joined_f5_enh_tissues = BedGetters._process_with_overlapping(flank_range, f5_enh_tissues).result
        return [
            BedsMutualOperation([tss_genes, joined_f5_enh_tissues], operation="intersection", **{"wa": True}).result]

    @staticmethod
    def _process_with_overlapping(flank_range, beds):
        if flank_range == "all":
            return BedsMutualOperation(beds, operation="intersection")
        elif flank_range == "any":
            return BedsMutualOperation(beds, operation="union")
        else:
            return []

    @staticmethod
    def _get_joined_flanked_genes(genes, genome, upstream, downstream):
        genome = convert_genome_name(genome, desirable_older_format="hg19")
        genes_beds = BedGetters.get_genes_beds(genes, genome)
        flanked_genes_beds = BedsFlanker(genes_beds, downstream, upstream, genome).results
        return BedsMutualOperation(flanked_genes_beds, operation="union").result

    @staticmethod
    def _get_tissues_beds(tissues):
        return [g.tissues_registry.get_bed(tissue) for tissue in tissues]

    @staticmethod
    def get_enhancers_encode_beds(tissues, genome, flank_range, upstram, downstream):
        return []

    @staticmethod
    def get_accessible_chromatin_encode_beds(tissues, genome, flank_range, upstream, downstream):
        return []


class BedsCollector:
    genes_params = ["genes", "genome"]

    transcription_fantom5_params = [
        "genes", "genome",
        "transcription-fantom5-range",
        "transcription-fantom5-kbs-upstream",
        "transcription-fantom5-kbs-downstream"
    ]

    enhancers_fantom5_params = [
        "genes", "tissues", "genome",
        "enhancers-fantom5-range",
        "enhancers-fantom5-kbs-upstream",
        "enhancers-fantom5-kbs-downstream"
    ]

    enhancers_encode_params = [
        "tissues", "genome",
        "enhancers-encode-range",
        "enhancers-encode-kbs-upstream",
        "enhancers-encode-kbs-downstream"
    ]

    accessible_chromatin_encode_params = [
        "tissues", "genome",
        "accessible-chromatin-encode-range",
        "accessible-chromatin-encode-kbs-upstream",
        "accessible-chromatin-encode-kbs-downstream"
    ]

    def __init__(self, data):
        self._data = data

    def collect_bed_files(self):
        bed_files = OrderedDict([
            ("genes",
             self._get_bed_files(
                 self.genes_params,
                 BedGetters.get_genes_beds)
             ),
            ("transcription-fantom5",
             self._get_bed_files(
                 self.transcription_fantom5_params,
                 BedGetters.get_transcription_starting_sites_fantom5_beds)
             ),
            ("enhancers-fantom5",
             self._get_bed_files(
                 self.enhancers_fantom5_params,
                 BedGetters.get_enhancers_fantom5_beds)
             ),
            ("enhancers-encode",
             self._get_bed_files(
                 self.enhancers_encode_params,
                 BedGetters.get_enhancers_encode_beds)
             ),
            ("accessible-chromatin-fantom5",
             self._get_bed_files(
                 self.accessible_chromatin_encode_params,
                 BedGetters.get_accessible_chromatin_encode_beds)
             )
        ])
        return bed_files

    def _get_bed_files(self, params, getter_method):
        params_values = [self._data.get(p) for p in params]
        if all(params_values):
            return getter_method(*params_values)
        else:
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

# TODO tu flankujesz promotory, a nam chodzi o znalezienie promotorów w rejonach flankujacych geny,
# TODO trzebaby dodac jeszcze jeden arg dla tej funkcji - genes/transcripts
# TODO tu flankujesz enhancery, a nam chodzi o znalezienie enhancerów w rejonach flankujacych geny,
# TODO wiec trzebaby dodac jeszcze jeden arg dla tej funkcji - genes/transcripts
