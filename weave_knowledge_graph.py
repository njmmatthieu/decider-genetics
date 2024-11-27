#!/usr/bin/env python3
import argparse
import logging
import pandas as pd

from biocypher import BioCypher
import ontoweaver

if __name__ == "__main__":

    usage = f"Extract nodes and edges from CSV tables from synthetic clinical, copy number variations and single nucleotide variants and prepare a knowledge graph import script."
    parser = argparse.ArgumentParser(
        description=usage)

    parser.add_argument("-cli", "--synthetic_clinical", metavar="CSV", nargs="+",
                        help="Extract from the synthetic clinical CSV file.")

    parser.add_argument("-cna", "--synthetic_cns", metavar="CSV", nargs="+",
                        help="Extract from the synthetic copy number alterations CSV file.")

    parser.add_argument("-snv", "--synthetic_variants", metavar="CSV", nargs="+",
                        help="Extract from the synthetic single nucleotide variants CSV file.")
    
    parser.add_argument("-o", "--oncokb", metavar="CSV", nargs="+",
                        help="Extract from an OncoKB CSV file.")
    
    levels = {
        "DEBUG": logging.DEBUG,
        "INFO": logging.INFO,
        "WARNING": logging.WARNING,
        "ERROR": logging.ERROR,
        "CRITICAL": logging.CRITICAL
    }

    parser.add_argument("-v", "--verbose", choices = levels.keys(), default = "WARNING",
                        help="Set the verbose level (default: %(default)s).")
    
    asked = parser.parse_args()
    bc = BioCypher(
        biocypher_config_path="config/biocypher_config.yaml",
        schema_config_path="config/schema_config.yaml"
    )

    # Actually extract data.
    nodes = []
    edges = []

    data_mappings = {}

    # Extract from databases not requiring preprocessing.
    if asked.synthetic_clinical:
        logging.info(f"Weave synthetic clinical data...")
        for file_path in asked.synthetic_clinical:
            data_mappings[file_path] =  "./decider_genetics/adapters/clinical.yaml"

    if asked.synthetic_cns:
        logging.info(f"Weave synthetic copy number alterations data")
        for file_path in asked.synthetic_cns:
            data_mappings[file_path] =  "./decider_genetics/adapters/cn_genes.yaml"

    if asked.synthetic_variants:
        logging.info(f"Weave synthetic single nucleotide variants data")
        for file_path in asked.synthetic_variants:
            data_mappings[file_path] = "./decider_genetics/adapters/all_variants.yaml"
    
    if asked.oncokb:
        logging.info(f"Weave OncoKB data...")
        for file_path in asked.oncokb:
            data_mappings[file_path] =  "./decider_genetics/adapters/oncokb.yaml"

    # Write everything.
    n, e = ontoweaver.extract(data_mappings)
    nodes += n
    edges += e

    import_file = ontoweaver.reconciliate_write(nodes, edges, "config/biocypher_config.yaml", "config/schema_config.yaml", separator=", ")

    print(import_file)