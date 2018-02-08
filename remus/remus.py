import logging
import os

import pandas as pd
from flask import Flask, g, render_template, jsonify, request

from remus.bio.bed.beds_loading import BedLoader
from remus.bio.bed.beds_operations import BedsOperation
from remus.bio.genes.genes_registry import GenesDBRegistry
from remus.bio.tissues.tissues_registry import TissuesFilesRegistry

app = Flask(__name__)


@app.before_request
def setup_registries():
    g.tissues_registry = TissuesFilesRegistry()
    g.genes_registry = GenesDBRegistry()


@app.after_request
def teardown_registries(response):
    g.genes_registry.teardown_registry()
    return response


@app.route("/")
def index():
    return render_template('index.html')


@app.route("/genomes")
def list_genomes():
    available_genomes = g.genes_registry.available_genomes
    return jsonify(available_genomes)


@app.route("/genes")
def list_matching_genes():
    genome_name = request.args.get("genome")
    pattern = request.args.get("pattern")
    genes_names = []
    if pattern and (genome_name in g.genes_registry.available_genomes):
        limit = request.args.get("limit", default=100, type=int)
        genes_names = g.genes_registry.get_matching_genes(genome_name, pattern, limit)
    return jsonify(genes_names)


@app.route("/tissues")
def list_tissues():
    return jsonify(g.tissues_registry.available_tissues)


@app.route("/operations")
def list_available_operations():
    return jsonify(list(BedsOperation.operations.keys()))


@app.route("/perform", methods=["GET", "POST"])
def perform_operation():
    params = get_operation_params()
    if all(params.values()):
        genes_beds = [get_genes_bed(params["genome"], gene) for gene in params["genes"] if gene]
        tissues_beds = [get_tissue_bed(tissue) for tissue in params["tissues"] if tissue]
        all_beds = genes_beds + tissues_beds
        try:
            processor = BedsOperation(all_beds, operation=params["operation"])
        except Exception as e:
            logging.exception("")
            return "Problem occured: " + str(e)
        data = {
            "Time elapsed (s)": processor.time_elapsed,
            "Result length (lines)": len(processor.result),
            "Coverage (bps)": processor.result.total_coverage()
        }
        summary = pd.DataFrame(data, index=[0])
        return summary.to_html()


def get_tissue_bed(tissue_name):
    tissue_path = os.path.join("db/tissues", "{}.bed".format(tissue_name))
    loader = BedLoader(tissue_path)
    return loader.bed


def get_genes_bed(genome, gene):
    return g.genes_registry.get_bed(genome, gene)


def get_operation_params():
    return {
        'genome': request.form.get("genome", None),
        'genes': request.form.getlist("genes", None),
        'tissues': request.form.getlist("tissues", None),
        'operation': request.form.get("operation", None)
    }
