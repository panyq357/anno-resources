#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Date: 2024-01-11
# Author: panyq357
# Description: Make GO annotation long table from JGI annotation info file.

from pathlib import Path

import pandas as pd

from onto_wrapper import Onto

config = {
    "ath_go": "ATH_GO_GOSLIM.txt",
    "go_owl": "http://purl.obolibrary.org/obo/go.owl",
    "gene_id_regex": r"AT[1-5CM]G\d{5}",
    "go_id_regex": r"GO:\d{7}",
    "out_path": "results/ath_go.csv.gz"
}

def main():

    ath_go = pd.read_table(config["ath_go"], comment="!", header=None)

    go_onto = Onto(config["go_owl"])

    id_to_onto = ath_go[[0, 5]]

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

    df.to_csv(out_path, index=False)
    print("Done")


if __name__ == "__main__":
    main()

