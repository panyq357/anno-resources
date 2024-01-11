#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Date: 2024-01-11
# Author: panyq357
# Description: Make GO annotation long table from JGI annotation info file.

from pathlib import Path

import pandas as pd

from onto_wrapper import Onto

config = {
    "jgi_si_annotation": "/mnt/d/PublicData/JGI/Setaria_italica_v2.2/rawdata/Sitalica_312_v2.2.annotation_info.txt",
    "go_owl": "http://purl.obolibrary.org/obo/go.owl",
    "gene_id_regex": r"Seita.[1-9]G\d{6}",
    "go_id_regex": r"GO:\d{7}",
    "out_path": "results/jgi-si-go-v2.2.xlsx"
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

    out_path = Path(config['out_path'])

    print(f"Writing to {str(out_path)} ...")

    if not out_path.parent.exists():
        out_path.parent.mkdir(parents=True)

    df.to_excel(out_path, index=False)
    print("Done")


if __name__ == "__main__":
    main()

