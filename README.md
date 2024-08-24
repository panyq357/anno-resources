Some scripts to prepare ontology long table for enrichment and GSEA analysis.

## Requirements

```bash
python3 -m pip install pandas owlready2 retrying openpyxl
```

## Fetch ontologies from Oryzabase

Just run this, annotation will be downloaded from [Oryzabase annotation](https://shigen.nig.ac.jp/rice/oryzabase/download/gene).

```bash
python3 scripts/get_ontologies_from_oryzabase.py
```

## Clean JGI Sitalica annotation info

First, download `Sitalica_312_v2.2.annotation_info.txt` from [JGI Sitalica annotation](https://data.jgi.doe.gov/refine-download/phytozome?organism=Sitalica)

Then replace the `jgi_si_annotation` path in `scripts/clean_jgi_si_annotation.py`.

Finally, run this.

```bash
python3 scripts/clean_jgi_si_annotation.py
```

## Clean GO annotation from TAIR

First download GO annotation from TAIR

```bash
aria2c "https://www.arabidopsis.org/api/download-files/download?filePath=GO_and_PO_Annotations/Gene_Ontology_Annotations/ATH_GO_GOSLIM.txt.gz"
gzip -d ATH_GO_GOSLIM.txt.gz
```

Then run the script.

```bash
python3 scripts/clean_ath_go.py
```

## Enrichment Analysis using clusterProfiler

```r
install.packages("readxl")
install.packages("BiocManager")
BiocManager::install("clusterProfiler")
```

```r
onto = readxl::read_excel("results/oryzabase-ontologies.xlsx", sheet="RAP_GO")

gene = c("Os01g0118100", "Os01g0549700", "Os02g0710800", "Os03g0108600", "Os03g0158200", "Os03g0746500")
universe = NULL

enrich_res = clusterProfiler::enricher(
    gene=gene,
    universe=universe,
    TERM2GENE=onto[c("OntoID", "GeneID")],
    TERM2NAME=onto[c("OntoID", "Description")]
)

write.csv(as.data.frame(enrich_res), "enrich_res.csv")

svg("demo_dotplot.svg")
clusterProfiler::dotplot(enrich_res)
dev.off()
```

![demo dotplot](demo_dotplot.svg)

