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

        print(f"Processing {table_name} ...")
        df = extract_oryzabase_onto_table(
            oryzabase_gene_list=oryzabase_gene_list,
            id_column_name=gene,
            id_regex=config["gene_id_regex"][gene],
            onto_column_name=onto,
            onto_regex=config["onto_id_regex"][onto],
            onto=onto_obj[onto]
        )

        if "GO" in table_name:
            df["Category"] = df["OntoID"].map(lambda x: get_go_category(x, onto_obj[onto]))

        df_dict[table_name] = df
        print("Done")

    print(f"Writing to {config['out_path']} ...")
    with pd.ExcelWriter(config["out_path"]) as writer:
        for df_name, df in df_dict.items():
            df.to_excel(writer, df_name, index=False)
    print("Done")


def id_extractor(id_regex):
    id_regex = re.compile(id_regex)
    def _id_ex(str_containing_id):
        id_list = id_regex.findall(str_containing_id)
        return list(set(id_list))
    return _id_ex


def extract_oryzabase_onto_table(
    oryzabase_gene_list,
    id_column_name,
    id_regex,
    onto_column_name,
    onto_regex,
    onto
) -> "Long DataFrame":
    '''
    By specifying the following:
        - name of gene ID containing column
        - gene ID regex
        - name of ontology id containing column
        - ontology id regex
        - Onto instence from onto_wrapper.py

    Extract a clean long table from oryzabase_gene_list that containing three columns:
        - gene ID (duplicated)
        - ontology ID of corresponding gene
        - description of ontology ID

    This DataFrame is suitable for enrichment analysis using R package: clusterProfiler.
    '''

    # Get a two column DataFrame: 1. gene ID ; 2. ontology info.
    row_mask = (oryzabase_gene_list[id_column_name].str.contains(id_regex, na=False) &
        oryzabase_gene_list[onto_column_name].str.contains(onto_regex, na=False))
    id_onto = oryzabase_gene_list.loc[row_mask , [id_column_name, onto_column_name]]

    # For each row, extract gene IDs and onto IDs into lists.
    # e.g.   "geneA,geneB" | "GO:1;GO:2"  ->  ["geneA", "geneB"]| ["GO:1", "GO:2"]
    id_onto[onto_column_name] = id_onto[onto_column_name].map(id_extractor(onto_regex))
    id_onto[id_column_name] = id_onto[id_column_name].map(id_extractor(id_regex))

    # Extend onto ID lists.
    id_onto[onto_column_name] = id_onto[onto_column_name].map(onto.extend_onto_id_list)

    # Explode to long table.
    id_onto = id_onto.explode(onto_column_name).explode(id_column_name)

    # Fetch ontology labels and sort by gene IDs.
    id_onto["Description"] = id_onto[onto_column_name].map(onto.get_onto_label)
    id_onto = id_onto.sort_values(by=id_column_name, ascending=True)

    # Unifiy column names.
    id_onto.columns = ["GeneID", "OntoID", "Description"]

    return id_onto


def get_go_category(go_id, onto):
    '''
    Determine which category a GO ID belongs to.
    '''

    if onto.has_ancestor(go_id, "GO:0008150"):
        return "BP"
    elif onto.has_ancestor(go_id, "GO:0003674"):
        return "MF"
    elif onto.has_ancestor(go_id, "GO:0005575"):
        return "MF"
    else:
        return None


if __name__ == "__main__":
    main()

