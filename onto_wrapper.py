# description: A wrapper class for common operation of ontologies
# date: 2023-09-30
# author: panyq357

import re

from owlready2 import get_ontology
from retrying import retry

class Onto():

    def __init__(self, owl_path_str):
        self.onto = self.__try_to_load_owl(owl_path_str)

    @retry
    def __try_to_load_owl(self, owl_path_str):
        print(f"Trying to load {owl_path_str} ...")
        ontology = get_ontology(owl_path_str).load()
        print("Done")
        return ontology

    def update_onto_id(self, onto_id):
        '''
        Check if a ontology ID is obsolete (e.g. GO:0030722),
        if it is, replace it with the new one.
        '''

        search_result_list = self.onto.search(id=onto_id)

        # The correct search result is a list with only one element.
        if len(search_result_list) < 1:  # Not Found, maybe a obsolete ID.
            search_result_list = self.onto.search(hasAlternativeId=onto_id)  # Search for up to date ID.
            if len(search_result_list) < 1:
                raise Exception(f"Error: {onto_id} not found!")
        if len(search_result_list) > 1:  # Found multiple ID (may never happen).
            raise Exception(f"Error: search {onto_id} found multiple result {search_result_list}")

        id_list = search_result_list[0].id

        # The correct ID list is a list with only one string element.
        if len(id_list) != 1:
            raise Exception(f"Error on getting ID of {onto_id}")

        return id_list[0]

    def get_ancestor_id_list(self, onto_id):
        '''
        Given a ontology ID (e.g. GO:0080050),
        extract a list of ontology ID of its ancestors.
        '''
        search_result_list = self.onto.search(id=onto_id)

        # The correct search result is a list with only one element.
        if len(search_result_list) != 1:
            raise Exception(f"Error on searching {onto_id}, maybe it is obsolete?")

        ancestor_list = search_result_list[0].ancestors()
        ancestor_id_list = list()
        for ancestor in ancestor_list:
            if hasattr(ancestor, "id") and len(ancestor.id) == 1:
                ancestor_id_list.append(ancestor.id[0])
            else:
                # Thing, BFO:0000001, BFO:0000003, BFO:0000015 have no ID, but they are everyone's ancestor.
                continue

        return ancestor_id_list


    def extend_onto_id_list(self, onto_id_list):
        '''
        Given a ontology ID list (e.g. ["GO:0080050", "GO:0010468"]),
        update it if necessary, and return a list appended with ancestor ontology ID.
        '''

        extended_list = [self.update_onto_id(onto_id) for onto_id in onto_id_list]

        for onto_id in extended_list.copy():
            ancestor_id_list = self.get_ancestor_id_list(onto_id)
            extended_list.extend(ancestor_id_list)
        return sorted(list(set(extended_list)))

    def get_onto_label(self, onto_id):
        '''
        Given a ontology ID,
        return it's label.
        '''

        search_result_list = self.onto.search(id=onto_id)

        if len(search_result_list) != 1:
            raise Exception(f"error on searching {onto_id}")

        label_list = list(set(search_result_list[0].label))

        if len(label_list) != 1:
            raise Exception(f"error on getting the label of {onto_id}")

        return label_list[0]

    def has_ancestor(self, onto_id, ancestor_id):
        '''
        Given a ontology ID, check if it has a specific ancestor.
        '''
        ancestor_id_list = self.get_ancestor_id_list(onto_id)

        if ancestor_id in ancestor_id_list:
            return True
        else:
            return False

    def get_clean_onto_table(self, id_to_onto, id_regex, onto_regex):
        '''
        By specifying the following:
            - A two column DataFrame, 1st column contains gene ID, 2nd column contains Onto ID
            - Gene ID regex
            - Ontology id regex
            - Onto instence from onto_wrapper.py

        Extract a clean long table from id_to_onto that containing three columns:
            - Gene ID (duplicated)
            - Ontology ID of corresponding gene
            - Description of ontology ID

        This DataFrame is suitable for enrichment analysis using R package: clusterProfiler.
        '''

        def id_extractor(id_regex):
            '''
            Generate a function, which do something like this:
                "GeneA,GeneB|GeneC" -> ["GeneA", "GeneB", "GeneC"]
            '''
            id_regex = re.compile(id_regex)
            def _id_ex(str_containing_id):
                id_list = id_regex.findall(str_containing_id)
                return list(set(id_list))
            return _id_ex

        id_to_onto.columns = ["GeneID", "OntoID"]

        # Pre-filter rows of two column DataFrame.
        row_mask = (id_to_onto["GeneID"].str.contains(id_regex, na=False) &
            id_to_onto["OntoID"].str.contains(onto_regex, na=False))
        id_to_onto = id_to_onto.loc[row_mask,:]

        # For each row, extract gene IDs and onto IDs into lists.
        # e.g.   "geneA,geneB" | "GO:1;GO:2"  ->  ["geneA", "geneB"]| ["GO:1", "GO:2"]
        id_to_onto.loc[:,"GeneID"] = id_to_onto["GeneID"].map(id_extractor(id_regex))
        id_to_onto.loc[:,"OntoID"] = id_to_onto["OntoID"].map(id_extractor(onto_regex))

        # Extend onto ID lists.
        id_to_onto.loc[:,"OntoID"] = id_to_onto["OntoID"].map(self.extend_onto_id_list)

        # Explode to long table.
        id_to_onto = id_to_onto.explode("GeneID").explode("OntoID")

        # Fetch ontology labels and sort by gene IDs.
        id_to_onto["Description"] = id_to_onto["OntoID"].map(self.get_onto_label)
        id_to_onto = id_to_onto.sort_values(by="GeneID", ascending=True)

        return id_to_onto

    def get_go_category(self, go_id):
        '''
        Determine which category a GO ID belongs to.
        '''

        if self.has_ancestor(go_id, "GO:0008150"):
            return "BP"
        elif self.has_ancestor(go_id, "GO:0003674"):
            return "MF"
        elif self.has_ancestor(go_id, "GO:0005575"):
            return "CC"
        else:
            return None

