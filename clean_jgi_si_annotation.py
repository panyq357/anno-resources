#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Date: 2023-07-10
# Author: panyq
# Description: Make GO, PO, TO annotation long table from oryzabase downloadable file.

import pandas as pd

from onto_wrapper import Onto

config = {
    "jgi_si_annotation": "~/Downloads/Sitalica_312_v2.2.annotation_info.txt",
    "go_owl": "http://purl.obolibrary.org/obo/go.owl",
    "gene_id_regex": r"Seita.[1-9]G\d{6}",
    "go_id_regex": r"GO:\d{7}",
    "out_path": "jgi-si-go-v2.2.xlsx"
}

def main():

    si_anno = pd.read_table(config["jgi_si_annotation"])

    go_onto = Onto(config["go_owl"])

    id_to_onto = si_anno[["locusName", "GO"]]

    df = go_onto.get_clean_onto_table(
        id_to_onto=id_to_onto,
        id_regex=config["gene_id_regex"],
        onto_regex=config["go_id_regex"]
    )

    df["Category"] = df["OntoID"].map(go_onto.get_go_category)

    print(f"Writing to {config['out_path']} ...")
    df.to_excel(config['out_path'], index=False)
    print("Done")


if __name__ == "__main__":
    main()

