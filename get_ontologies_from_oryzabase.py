#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Date: 2023-07-10
# Author: panyq
# Description: Make GO, PO, TO annotation long table from oryzabase downloadable file.

import datetime
import re

import pandas as pd
from owlready2 import get_ontology
from retrying import retry

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
        "Gene Ontology": try_to_load_owl(config["url"]["go_owl"]),
        "Plant Ontology": try_to_load_owl(config["url"]["po_owl"]),
        "Trait Ontology": try_to_load_owl(config["url"]["to_owl"])
    }

    df_dict = {}

    for table_name in config["table"]:

        gene = config["table"][table_name]["gene"]
        onto = config["table"][table_name]["onto"]

        print(f"Processing {table_name} ...")
        df_dict[table_name] = extract_oryzabase_onto_table(
            oryzabase_gene_list=oryzabase_gene_list,
            id_column_name=gene,
            id_regex=config["gene_id_regex"][gene],
            onto_column_name=onto,
            onto_regex=config["onto_id_regex"][onto],
            onto=onto_obj[onto]
        )
        print("Done")

    print(f"Writing to {config['out_path']} ...")
    with pd.ExcelWriter(config["out_path"]) as writer:
        for df_name, df in df_dict.items():
            df.to_excel(writer, df_name, index=False)
    print("Done")


def update_onto_id(onto_id, onto):
    '''
    Check if a ontology ID is obsolete,
    if it is, replace it with the new one.
    '''

    search_result_list = onto.search(id=onto_id)
    if len(search_result_list) != 1:
        search_result_list = onto.search(hasAlternativeId=onto_id)
        if len(search_result_list) != 1:
            raise Exception(f"error on searching {onto_id}")
    elif len(search_result_list[0].id) != 1:
        raise Exception(f"error on getting ID of {ancestor}")

    return search_result_list[0].id[0]


def get_ancestor_id_list(onto_id, onto):
    '''
    Given a ontology ID (e.g. GO:0080050),
    extract a list of ontology ID of its ancestors from owlready2 ontology object.
    '''
    search_result_list = onto.search(id=onto_id)
    if len(search_result_list) != 1:
        raise Exception(f"error on searching {onto_id}")
    ancestor_list = search_result_list[0].ancestors()
    ancestor_id_list = list()
    for ancestor in ancestor_list:
        if hasattr(ancestor, "id") and len(ancestor.id) == 1:
            ancestor_id_list.append(ancestor.id[0])
        else:
            continue
    return ancestor_id_list


def extend_onto_id_list(onto_id_list, onto):
    '''
    Given a ontology ID list (e.g. ["GO:0080050", "GO:0010468"])
    return a list appended with ancestor ontology ID.
    '''
    extended_list = [update_onto_id(onto_id, onto) for onto_id in onto_id_list]
    for onto_id in extended_list.copy():
        ancestor_id_list = get_ancestor_id_list(onto_id, onto)
        extended_list.extend(ancestor_id_list)
    return sorted(list(set(extended_list)))


def get_onto_label(onto_id, onto):
    '''
    Given a ontology ID,
    return it's label.
    '''
    search_result_list = onto.search(id=onto_id)
    if len(search_result_list) != 1:
        raise Exception(f"error on searching {onto_id}")
    label_list = list(set(search_result_list[0].label))
    if len(label_list) != 1:
        raise Exception(f"error on getting the label of {onto_id}")
    return label_list[0]


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
        - owlready2 ontology object

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
    onto_ex = id_extractor(onto_regex)
    id_onto[onto_column_name] = id_onto[onto_column_name].map(onto_ex)
    id_ex = id_extractor(id_regex)
    id_onto[id_column_name] = id_onto[id_column_name].map(id_ex)

    # Extend onto ID lists.
    id_onto[onto_column_name] = id_onto[onto_column_name].map(lambda x: extend_onto_id_list(x, onto))

    # Explode to long table.
    id_onto = id_onto.explode(onto_column_name).explode(id_column_name)

    # Fetch ontology labels and sort by gene IDs.
    id_onto["Description"] = id_onto[onto_column_name].map(lambda x: get_onto_label(x, onto))
    id_onto = id_onto.sort_values(by=id_column_name, ascending=True)

    # Unifiy column names.
    id_onto.columns = ["GeneID", "OntoID", "Description"]

    return id_onto


@retry
def try_to_load_owl(owl_path_str):
    print(f"Trying to load {owl_path_str} ...")
    ontology = get_ontology(owl_path_str).load()
    print("Done")
    return ontology


if __name__ == "__main__":
    main()

