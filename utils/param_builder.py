


def build_gdc_params(dataset, data_normalizaton):
    gene_expression_file_name = "TCGA-{}.htseq_{}.tsv".format(dataset, data_normalizaton)
    pval_preprocessing_file_name = "pvals_protein_coding_{}.txt".format(data_normalizaton)
    phenotype_file_name = "TCGA-{}.GDC_phenotype.tsv".format(dataset)

    return gene_expression_file_name, phenotype_file_name, pval_preprocessing_file_name