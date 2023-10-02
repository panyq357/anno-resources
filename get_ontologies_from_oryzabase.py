#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Date: 2023-07-10
# Author: panyq
# Description: Make GO, PO, TO annotation long table from oryzabase downloadable file.

import datetime
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
        "MSU_TO": {"gene": "MSU ID", "onto": "Trait Ontology"},
    },
    "out_path": f"oryzabase-ontologies-{str(datetime.date.today())}.xlsx"
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
        df = onto_obj[onto].get_clean_onto_table(
            id_to_onto,
            id_regex=config["gene_id_regex"][gene],
            onto_regex=config["onto_id_regex"][onto],
        )

        if "GO" in table_name:
            df["Category"] = df["OntoID"].map(onto_obj[onto].get_go_category)

        df_dict[table_name] = df
        print("Done")

    print(f"Writing to {config['out_path']} ...")
    with pd.ExcelWriter(config["out_path"]) as writer:
        for df_name, df in df_dict.items():
            df.to_excel(writer, df_name, index=False)
    print("Done")

if __name__ == "__main__":
    main()

