#!/usr/bin/env python3

import os
import requests

"""
available_args = ["genome", "genes", "tissues"]
available_args += ["genes-select-include-gene-transcripts"]
available_args += ["transcription-fantom5-used",
                        "transcription-fantom5-combine-mode",
                        "transcription-fantom5-kbs-upstream",
                        "transcription-fantom5-kbs-downstream"]

available_args += ["mirna-mirtarbase-used",
                        "mirna-mirwalk-used",
                        "mirna-targets-combine-mode",
                        "mirna-mirtarbase-include-weak",
                        "mirna-mirwalk-minimal-confidence"]
"""

def get_url(endpoint, operation):
    suffix_dict = {'genes':'api/genes',
                   'tissues':'api/tissues',
                   'perform':'api/perform',
                   'download':'api/download_last'}
    return os.path.join(endpoint, suffix_dict[operation])


def query(endpoint, request_args):
    url = get_url(endpoint,"perform")
    resp = requests.post(url, request_args)
    if resp is None or resp.text == "Error occurred" or resp.headers.get('Set-cookie') is None:
        raise Exception("Erronous response for query to %s, with args: %s" % (url, request_args))

    cookies = {e[0]: e[1] for e in [e.split("=") for e in resp.headers.get("Set-Cookie").split(";")] if len(e) > 1}
    resp = requests.get(get_url(endpoint, 'download'), cookies=cookies)

    return resp.content


def query_default(endpoint, genome, genes, tissues,
                  enh_enc=True, enh_f5=True, chrom=True):

    request_args = {'genes': genes, 'tissues': tissues, 'genome': genome,
                    "genes-select-include-gene-transcripts": None}
    if enh_enc:
        request_args.update({"enhancers-encode-used":          "yes",
                             "enhancers-encode-combine-mode":  "any",
                             "enhancers-encode-kbs-upstream":   500,
                             "enhancers-encode-kbs-downstream": 500})
    if enh_f5:
        request_args.update({"enhancers-fantom5-used":           "yes",
                             "enhancers-fantom5-combine-mode":   "any",
                             "enhancers-fantom5-kbs-upstream":   500,
                             "enhancers-fantom5-kbs-downstream": 500})
    if chrom:
        request_args.update({"accessible-chromatin-encode-used":          "yes",
                             "accessible-chromatin-encode-combine-mode":  "any",
                             "accessible-chromatin-encode-kbs-upstream":   500,
                             "accessible-chromatin-encode-kbs-downstream": 500})

    return query(endpoint, request_args)


def get_all_tissues(endpoint, genome):
    return requests.get(get_url(endpoint, 'tissues'), {"genome": genome}).json()


def get_matching_tissues(endpoint, genome, patterns):
    tissues = set()
    for p in patterns:
        ts = requests.get(get_url(endpoint, 'tissues'),
                          {"pattern": p, "genome": genome}).json()
        tissues.update(ts)

    return list(tissues)


def get_matching_tissues_strict(endpoint, genome, patterns):
    """
    Requires perfect match of the tissue name (dataset information in parenthesis is stripped)
    """
    tissues = get_matching_tissues(endpoint, genome, patterns)
    return [t for t in tissues if (t.split("(")[0]).strip() in patterns]


def get_tissues_from_file(f):
    return [l.strip() for l in f.readlines() if l.strip() != ""]


if __name__ == '__main__':

    import sys, argparse
    parser = argparse.ArgumentParser(description="Remus API client")
    parser.add_argument("--genome", choices=['hg19', 'hg38'], default='hg19', help="(Required) Genome build to use. Default [hg19]")
    parser.add_argument("--genes", required=True, metavar="gene1,gene2,...", help="(Required) Comma-separated list of target gene symbols")
    parser.add_argument("--tissue_keywords", metavar="kidney,liver,...",
                        help="Comma-separated list of keywords for tissue selection. \
                              Tissues matching any of the keywords will be used.")
    parser.add_argument("--tissue_file", type=argparse.FileType('r'), help="Filename with tissues to use.")
    parser.add_argument("--endpoint", default="http://0.0.0.0:5000", help="Remus endpoint URL. Default [http://0.0.0.0:5000] ")
    parser.add_argument("--verbose", action="store_true", help='Print control messages to stderr')

    args = parser.parse_args()

    if not args.tissue_keywords and not args.tissue_file:
        raise Exception("No tissues selected. Either --tissue_keywords or --tissue_file must be provided. \
                         To select all tissues, use keyword \'ALL\'.")

    genes = args.genes.split(",")
    tissues = \
        get_all_tissues(args.endpoint, args.genome) if args.tissue_keywords == 'ALL' else \
        get_matching_tissues(args.endpoint, args.genome, args.tissue_keywords.split(',')) if args.tissue_keywords else \
        get_tissues_from_file(args.tissue_file)

    if args.tissue_keywords and tissues == []:
        sys.stderr.write('\n'.join(['Could not find any tissues matching: '] + args.tissue_keywords.split(',')) + '\n')
        sys.exit(1)
    elif args.verbose:
        sys.stderr.write('\n'.join(['Tissues used: '] + sorted(tissues)) + '\n')

    print(
        query_default(args.endpoint,
                  args.genome,
                  genes,
                  tissues).decode()
    )