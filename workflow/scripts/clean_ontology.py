from pathlib import Path

import pandas as pd

from onto_wrapper import Onto

raw_df = eval(snakemake.params[0]["reader"])(snakemake.input[0])

onto = Onto(snakemake.params[0]["owl"])

id_to_onto = raw_df[[snakemake.params[0]["gene_id_col"], snakemake.params[0]["onto_id_col"]]]

term_to_gene = onto.get_clean_onto_table(
    id_to_onto = id_to_onto,
    id_regex = snakemake.params[0]["gene_id_regex"],
    onto_regex = snakemake.params[0]["onto_id_regex"]
)[["OntoID", "GeneID"]]

term_to_description = pd.DataFrame({"OntoID": term_to_gene["OntoID"].unique()})

term_to_description["Description"] = term_to_description["OntoID"].map(onto.get_onto_label)

if "GO" in snakemake.params[0]["onto_id_regex"]:
    term_to_description["Description"] = term_to_description["Description"] + " (" + term_to_description["OntoID"].map(onto.get_go_category) + ")"

term_to_description = term_to_description.sort_values("OntoID")

term_to_gene.to_csv(snakemake.output["term_to_gene"], sep="\t", header = False, index = False)
term_to_description.to_csv(snakemake.output["term_to_description"], sep="\t", header = False, index = False)
