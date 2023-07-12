library(readxl)
library(clusterProfiler)

config = list(
    oryzabase_ontologies_xlsx_path = "oryzabase-ontologies-2023-07-10.xlsx",
    sheet_list = list("RAP_GO", "RAP_TO", "RAP_PO"),
    p_column_name = "pvalue",
    p_list = list(0.001),
    fc_column_name = "log2FoldChange",
    fc_list = list(log2(1)),
    fc_rev_fun = function(x) { 2 ^ x },
    outdir = "results/oryzabase_anno_pipeline"
)

de_result_df = read.csv("deseq_results.csv", row.names=1)


main = function() {

    if (! dir.exists(config[["outdir"]])) dir.create(config[["outdir"]], recursive=T)

    onto_list = read_oryzabase_onto(config[["oryzabase_ontologies_xlsx_path"]], sheet_list=config[["sheet_list"]])

    for (p in config[["p_list"]]) {
        for (fc in config[["fc_list"]]) {

            gene_list = select_gene_from_de_result(de_result_df, config[["p_column_name"]], p, config[["fc_column_name"]], fc)
            gsea_gene_list = get_gsea_gene_from_de_result(de_result_df, config[["fc_column_name"]])

            enrich_result_list = oryzabase_enricher(gene_list, onto_list)
            gsea_result_list = oryzabase_GSEA(gsea_gene_list, onto_list)
            names(gsea_result_list) = paste0(names(gsea_result_list), "_GSEA")

            all_result_list = c(enrich_result_list, gsea_result_list)

            save_result_list(all_result_list, file.path(config[["outdir"]], sprintf("p%s.fc%s.xlsx", p, config$fc_rev_fun(fc))))
            save_dotplot(all_result_list, file.path(config[["outdir"]], sprintf("p%s.fc%s.pdf", p, config$fc_rev_fun(fc))))
            saveRDS(all_result_list, file.path(config[["outdir"]], sprintf("p%s.fc%s.rds", p, config$fc_rev_fun(fc))))
        }
    }
}


select_gene_from_de_result = function(de_result_df, p_column_name, max_p, fc_column_name, min_fc) {
    gene_list = list()
    class(gene_list) = "gene_list"
    gene_list[["target"]] = de_result_df |>
        subset(de_result_df[[p_column_name]] <= max_p) |>
        subset(de_result_df[[fc_column_name]] > min_fc) |>
        row.names()
    gene_list[["universe"]] = row.names(de_result_df)
    return(gene_list)
}

get_gsea_gene_from_de_result = function(de_result_df, fc_column_name) {
    gsea_gene_list = de_result_df[[fc_column_name]] |>
        setNames(row.names(de_result_df)) |>
        sort(decreasing=T)
    class(gsea_gene_list) = c("numeric", "gsea_gene_list")
    return(gsea_gene_list)
}

read_oryzabase_onto = function(xlsx_path, sheet_list) {
    onto_list = list()
    class(onto_list) = "onto_list"
    for (sheet in sheet_list)
        onto_list[[sheet]] = read_excel(xlsx_path, sheet=sheet)
    return(onto_list)
}

oryzabase_enricher = function(gene_list, onto_list) {
    enrich_result_list = list()
    class(enrich_result_list) = "enrich_result_list"
    for (onto_name in names(onto_list)) {
        res = enricher(
            gene=gene_list[["target"]],
            universe=gene_list[["universe"]],
            TERM2GENE=onto_list[[onto_name]][c("OntoID", "GeneID")],
            TERM2NAME=onto_list[[onto_name]][c("OntoID", "Description")],
            pvalueCutoff = 1,
            qvalueCutoff = 1
        )
        if (is.null(res)) {
            warning(sprintf("%s result is NULL.", onto_name))
        } else {
            enrich_result_list[[onto_name]] = res
        }
    }
    return(enrich_result_list)
}

oryzabase_GSEA = function(gsea_gene_list, onto_list) {
    gsea_result_list = list()
    class(gsea_result_list) = "gsea_result_list"
    for (onto_name in names(onto_list)) {
        res = GSEA(
            geneList=gsea_gene_list,
            TERM2GENE=onto_list[[onto_name]][c("OntoID", "GeneID")],
            TERM2NAME=onto_list[[onto_name]][c("OntoID", "Description")],
            pvalueCutoff = 1
        )
        if (is.null(res)) {
            warning(sprintf("%s result is NULL.", onto_name))
        } else {
            gsea_result_list[[onto_name]] = res
        }
    }
    return(gsea_result_list)
}

save_result_list = function(result_list, out_xlsx) {
    result_list |>
        lapply(as.data.frame) |>
        writexl::write_xlsx(out_xlsx)
}

save_dotplot = function(result_list, out_pdf) {
    pdf(out_pdf)
    for (res_name in names(result_list))
        dotplot(result_list[[res_name]], title=res_name) |> plot()
    dev.off()
}

main()
