#!/usr/bin/env python3
import argparse
import logging
import pandas as pd
import yaml

from biocypher import BioCypher
import decider_genetics.adapters.preprocessing_dataframes as preprocess
import ontoweaver

if __name__ == "__main__":

    usage = f"Extract nodes and edges from CSV tables from synthetic clinical, copy number variations and single nucleotide variants and prepare a knowledge graph import script."
    parser = argparse.ArgumentParser(
        description=usage)

    parser.add_argument("-cli", "--synthetic_clinical", metavar="CSV", nargs="+",
                        help="Extract from the synthetic clinical CSV file.")

    parser.add_argument("-cna", "--synthetic_cna", metavar="CSV", nargs="+",
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
        
        clinical_df = pd.read_csv(asked.synthetic_clinical[0], sep=';')
        preprocessed_clinical_df = preprocess.preprocess_clinical(clinical_df)

        mapping_file = "./decider_genetics/adapters/clinical.yaml"
        with open(mapping_file) as fd:
            conf = yaml.full_load(fd)

        adapter = ontoweaver.tabular.extract_all(df=preprocessed_clinical_df, config=conf,separator = None, affix= "none")

        nodes += adapter.nodes
        edges += adapter.edges

        logging.info(f"Wove Clinical: {len(nodes)} nodes, {len(edges)} edges.")

    if asked.synthetic_cna:
        logging.info(f"Weave synthetic copy number alterations data")
        
        cna_df = pd.read_csv(asked.synthetic_cna[0], sep='\t')
        preprocessed_cna_df = preprocess.preprocess_cna(cna_df)

        mapping_file = "./decider_genetics/adapters/cna_genes.yaml"
        with open(mapping_file) as fd:
            conf = yaml.full_load(fd)

        adapter = ontoweaver.tabular.extract_all(df=preprocessed_cna_df, config=conf,separator = None, affix= "none")

        nodes += adapter.nodes
        edges += adapter.edges

        logging.info(f"Wove Clinical: {len(nodes)} nodes, {len(edges)} edges.")

    if asked.synthetic_variants:
        logging.info(f"Weave synthetic single nucleotide variants data")
        for file_path in asked.synthetic_variants:
            data_mappings[file_path] = "./decider_genetics/adapters/all_variants.yaml"
    
    if asked.oncokb:
        logging.info(f"Weave OncoKB data...")
        for file_path in asked.oncokb:
            data_mappings[file_path] =  "./decider_genetics/adapters/oncokb.yaml"

    # # Write everything.
    # n, e = ontoweaver.extract(data_mappings)
    # nodes += n
    # edges += e

    import_file = ontoweaver.reconciliate_write(nodes, edges, "config/biocypher_config.yaml", "config/schema_config.yaml", separator=", ")
    # bc.write_schema_info(as_node=True)

    print(import_file)