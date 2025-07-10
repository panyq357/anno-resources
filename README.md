Some scripts to prepare ontology long table for enrichment and GSEA analysis.

This pipeline requires these dependencies.

```bash
python3 -m pip install pandas owlready2 retrying openpyxl
```

following data require manually download.

- `Sitalica_312_v2.2.annotation_info.txt`: from [JGI Sitalica annotation](https://data.jgi.doe.gov/refine-download/phytozome?organism=Sitalica)
- `ATH_GO_GOSLIM.txt`: from [TAIR Download GO_and_PO_Annotations](https://www.arabidopsis.org/download/list?dir=GO_and_PO_Annotations%2FGene_Ontology_Annotations)

For TAIR gene search results, in [TAIR Gene Search](https://www.arabidopsis.org/search/genes) page, use different chromosome as filter,
after clicking "Submit Query" bottom and search results displayed, click "Download All" bottom to download result tsv file.


A example of enrichment analysis using clusterProfiler.

```r
install.packages("readxl")
install.packages("BiocManager")
BiocManager::install("clusterProfiler")
```

```r
term_to_gene = readr::read_tsv("results/oryzabase.RAP_GO.term_to_gene.tsv.gz", col_names=F)
term_to_description = readr::read_tsv("results/oryzabase.RAP_GO.term_to_description.tsv.gz", col_names=F)

gene = c("Os01g0118100", "Os01g0549700", "Os02g0710800", "Os03g0108600", "Os03g0158200", "Os03g0746500")
universe = NULL

enrich_res = clusterProfiler::enricher(
    gene=gene,
    universe=universe,
    TERM2GENE=term_to_gene,
    TERM2NAME=term_to_description
)

write.csv(as.data.frame(enrich_res), "enrich_res.csv")

svg("demo_dotplot.svg")
clusterProfiler::dotplot(enrich_res)
dev.off()
```

![demo dotplot](demo_dotplot.svg)

