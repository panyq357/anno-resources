# Oryzabase Ontologies

## Requirements

```bash
python3 -m pip install pandas owlready2 retrying openpyxl
```

## Fetch Ontologies from Oryzabase

```bash
python3 get_ontologies_from_oryzabase.py
```

## Enrichment Analysis using clusterProfiler

```r
onto = readxl::read_excel("oryzabase-ontologies-2023-05-27.xlsx", sheet="RAP_GO")

gene_list = c("Os01g0118100", "Os01g0549700", "Os02g0710800", "Os03g0108600", "Os03g0158200", "Os03g0746500")
universe = unique(onto[["GeneID"]])

enrich_result = clusterProfiler::enricher(
    gene=gene_list,
    universe=universe,
    TERM2GENE=onto[c("OntoID", "GeneID")],
    TERM2NAME=onto[c("OntoID", "Description")]
)

svg("demo_dotplot.svg")
clusterProfiler::dotplot(enrich_result)
dev.off()
```

![demo dotplot](demo_dotplot.svg)

