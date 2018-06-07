import re
from collections import OrderedDict

from flask import g

from remus.bio.bed.beds_operations import BedsMutualOperation, BedsFlanker


class BedGetters:
    @staticmethod
    def get_genes_beds(genes, genome):
        genome = convert_genome_name(genome, desirable_older_format="hg37")
        return [g.genes_registry.get_bed(genome, gene) for gene in genes]
        # return BedsMutualOperation(result, operation="union").result

    @staticmethod
    def get_tissues_beds(tissues):
        return [g.tissues_registry.get_bed(tissue) for tissue in tissues]

    @staticmethod
    def get_transcription_starting_sites_fantom5_beds(tissues, genome, flank_range, upstream, downstream):
        genome = convert_genome_name(genome, desirable_older_format="hg19")
        return BedsFlanker([g.tss_registry.get_bed()], downstream, upstream, genome).results

    @staticmethod
    def get_enhancers_fantom5_beds(tissues, genome, flank_range, upstream, downstream):
        flanked_tissues_beds = []
        genome = convert_genome_name(genome, desirable_older_format="hg19")
        tissues_beds = BedGetters.get_tissues_beds(tissues)
        flanked_tissues_beds.extend(BedsFlanker(tissues_beds, downstream, upstream, genome).results)
        if len(flanked_tissues_beds) <= 1:
            return flanked_tissues_beds
        elif flank_range == "all":
            return BedsMutualOperation(flanked_tissues_beds, operation="intersection").result
        elif flank_range == "any":
            return BedsMutualOperation(flanked_tissues_beds, operation="union").result
        else:
            return []

    @staticmethod
    def get_enhancers_encode_beds(tissues, genome, flank_range, upstram, downstream):
        return []

    @staticmethod
    def get_accessible_chromatin_encode_beds(tissues, genome, flank_range, upstream, downstream):
        return []


class BedsCollector:
    genes_params = [
        "genes", "genome"
    ]
    tissues_params = [
        "tissues"]
    transcription_fantom5_params = [
        "tissues", "genome",
        "transcription-fantom5-range",
        "transcription-fantom5-kbs-upstream",
        "transcription-fantom5-kbs-downstream"
    ]
    enhancers_fantom5_params = [
        "tissues", "genome",
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
            ("tissues",
             self._get_bed_files(
                 self.tissues_params,
                 BedGetters.get_tissues_beds)
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
    if pattern:
        return g.tissues_registry.get_matching_tissues(pattern, limit)
    else:
        return []


def convert_genome_name(genome, desirable_older_format="hg37"):
    if re.match("(hg37|hg19)", genome, re.IGNORECASE):
        return desirable_older_format
    else:
        return genome.lower()
