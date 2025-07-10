#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Date: 2023-07-10
# Author: panyq
# Description: Make GO, PO, TO annotation long table from oryzabase downloadable file.

import datetime
from pathlib import Path
import re

import pandas as pd

from onto_wrapper import Onto

config = {
    "url": {
        "oryzabase_gene_list": "https://shigen.nig.ac.jp/rice/oryzabase/gene/download?classtag=GENE_EN_LIST",
        "go_owl": "http://purl.obolibrary.org/obo/go.owl",
        "po_owl": "http://purl.obolibrary.org/obo/po.owl",
        "to_owl": "http://purl.obolibrary.org/obo/to.owl"
    },
    "gene_id_regex": {
        "RAP ID": r"Os\d{2}g\d{7}",
        "MSU ID": r"LOC_Os\d{2}g\d{5}"
    },
    "onto_id_regex": {
        "Gene Ontology": r"GO:\d{7}",
        "Plant Ontology": r"PO:\d{7}",
        "Trait Ontology": r"TO:\d{7}"
    },
    "table": {
        # "RAP ID", "Gene Ontology"... are the column names of OryzabaseGeneList file.
        "RAP_GO": {"gene": "RAP ID", "onto": "Gene Ontology"},
        "RAP_PO": {"gene": "RAP ID", "onto": "Plant Ontology"},
        "RAP_TO": {"gene": "RAP ID", "onto": "Trait Ontology"},
        "MSU_GO": {"gene": "MSU ID", "onto": "Gene Ontology"},
        "MSU_PO": {"gene": "MSU ID", "onto": "Plant Ontology"},
        "MSU_TO": {"gene": "MSU ID", "onto": "Trait Ontology"}
    },
    "out_dir": "results"
}

def main():

    print("Downloading OryzabaseGeneList ...")
    oryzabase_gene_list = pd.read_table(config["url"]["oryzabase_gene_list"])
    print("Done")

    onto_obj = {
        "Gene Ontology": Onto(config["url"]["go_owl"]),
        "Plant Ontology": Onto(config["url"]["po_owl"]),
        "Trait Ontology": Onto(config["url"]["to_owl"])
    }

    df_dict = {}

    for table_name in config["table"]:

        gene = config["table"][table_name]["gene"]
        onto = config["table"][table_name]["onto"]

        id_to_onto = oryzabase_gene_list[[gene, onto]]

        print(f"Processing {table_name} ...")
        term_to_gene = onto_obj[onto].get_clean_onto_table(
            id_to_onto,
            id_regex=config["gene_id_regex"][gene],
            onto_regex=config["onto_id_regex"][onto],
        )[["OntoID", "GeneID"]]

        term_to_description = pd.DataFrame({"OntoID": term_to_gene["OntoID"].unique()})
        term_to_description["Description"] = term_to_description["OntoID"].map(onto_obj[onto].get_onto_label)

        if "GO" in table_name:
            term_to_description["Description"] = term_to_description["Description"] + " (" + term_to_description["OntoID"].map(onto_obj[onto].get_go_category) + ")"

        term_to_description = term_to_description.sort_values("OntoID")

        df_dict[table_name] = {
            "term_to_gene": term_to_gene,
            "term_to_description": term_to_description
        }
        print("Done")


    out_dir = Path(config['out_dir'])
    print(f"Writing to {str(out_dir)} ...")

    if not out_dir.exists():
        out_dir.mkdir(parents=True)

    for k1, v1 in df_dict.items():
        for k2, v2 in v1.items():
            v2.to_csv(out_dir / f"oryzabase.{k1}.{k2}.tsv.gz", sep="\t", header = False, index = False)

    print("Done")

if __name__ == "__main__":
    main()

