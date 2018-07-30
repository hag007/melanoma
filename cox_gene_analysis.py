from lifelines.datasets import load_dd
from lifelines import KaplanMeierFitter
from matplotlib import pyplot as plt
data = load_dd()
from lifelines.datasets import load_rossi
from lifelines import CoxPHFitter
from infra import *
import pandas

FILTER_IN = "filter in"
FILTER_OUT = "filter out"
OS_FIELDS = ["_OS", "_OS_IND"]

pheno_start = 1
pheno_limit = 2000


def cox_gene(test_independently, filtered_out, filtered_in, filter_type, filter_na_by_rows,tested_genes_file_name ,gene_expression_file_name, phenotype_file_name, survival_file_name, filter_expression=None, var_th_index=None):
    tested_genes = load_gene_list(tested_genes_file_name)
    survival_dataset = np.array(load_survival_data(survival_file_name, survival_list_path=None))
    gene_expression_dataset = np.array(load_gene_expression_profile_by_genes(tested_genes_file_name, gene_expression_file_name))
    if filter_expression is not None:
        filtered_patients = divided_patient_ids_by_label(phenotype_file_name, groups=filter_expression)[0]
        print "number of filtered patients from phenotypes: {}".format(len(filtered_patients))
    else:
        filtered_patients = np.append(gene_expression_dataset[0,1:] ,survival_dataset[1:,0])
    filtered_gene_expression_bool = np.in1d(gene_expression_dataset[0], filtered_patients)
    filtered_gene_expression_bool[0]=True
    filtered_survival_bool = np.in1d(survival_dataset[:,0], filtered_patients)
    filtered_survival_bool[0] = True
    print "Total n patients in expression before filtering: {}".format(np.shape(gene_expression_dataset)[1] - 1)
    gene_expression_dataset = gene_expression_dataset[:,filtered_gene_expression_bool]
    print "Total n patients in expression after filtering: {}".format(np.shape(gene_expression_dataset)[1]-1)
    dataset_headers_rows, dataset_headers_columns, dataset = separate_headers(gene_expression_dataset)
    gene_expression_top_var, gene_expression_top_var_headers_rows, gene_expression_top_var_headers_columns = filter_top_var_genes(dataset ,dataset_headers_columns, dataset_headers_rows, var_th_index)

    # flipping
    tmp=gene_expression_top_var_headers_rows
    gene_expression_top_var_headers_rows = gene_expression_top_var_headers_columns
    gene_expression_top_var_headers_columns = tmp
    gene_expression_top_var = np.rot90(np.flip(gene_expression_top_var, 1), k=-1, axes=(1, 0))
    #

    expression_survival_integrated = {}
    for i, cur_expression in enumerate(gene_expression_top_var):
        expression_survival_integrated[gene_expression_top_var_headers_rows[i]] = cur_expression[:pheno_limit]
    for cur_survival in survival_dataset[1:]:
        if expression_survival_integrated.has_key(cur_survival[0]):
            expression_survival_integrated[cur_survival[0]] = list(expression_survival_integrated[cur_survival[0]]) + list(cur_survival[4:])
    for k, v in dict(expression_survival_integrated).iteritems():
        if len(v) != len(gene_expression_top_var_headers_columns[:pheno_limit]) + len(survival_dataset[0][4:]):
            expression_survival_integrated.pop(k, None)

    if filter_type==FILTER_IN:
        filtered_list = tested_genes
    else:
        tested_genes = load_gene_list("")
        filtered_list = filtered_out
    filtered_list = [x.split(".")[0] for x in filtered_list]

    if test_independently:
        for i, cur in enumerate(filtered_list):
            results = test_cox(gene_expression_top_var, gene_expression_top_var_headers_columns, survival_dataset, expression_survival_integrated, [cur]+ OS_FIELDS, filter_type, filter_na_by_rows)
            if i == 0:
                return results
            else:
                return results[results.index('\n')+1:]
    else:
        return test_cox(gene_expression_top_var, gene_expression_top_var_headers_columns, survival_dataset, expression_survival_integrated, filtered_list+ OS_FIELDS, filter_type, filter_na_by_rows)

def test_cox(phenotype_dataset, gene_expression_top_var_headers_columns, survival_dataset, pheno_survival_integrated, filtered_list, filter_type, filter_na_by_rows):
    headers = ["ids"] + list(gene_expression_top_var_headers_columns[:pheno_limit]) + list(survival_dataset[0][4:])
    pandas.set_option("mode.use_inf_as_na", True)
    df = pandas.DataFrame(columns=headers, data=[[k] + v for k, v in
                                                 pheno_survival_integrated.iteritems()])  # np.array().astype(np.float32)
    for cur_header in headers[1:]:
        if filter_type == FILTER_IN and cur_header.split(".")[0] not in filtered_list:
            df = df.drop(cur_header, 1)
            print "column {} was dropped as it's not filtered in".format(cur_header)
            continue

        if filter_type == FILTER_OUT and cur_header in filtered_list:
            df = df.drop(cur_header, 1)
            print "column {} was dropped as it has low variance".format(cur_header)
            continue

        if df[[cur_header]].isnull().values.any() and (not filter_na_by_rows or filter_type == FILTER_OUT):
            df = df.drop(cur_header, 1)
            print "column {} was dropped as it has NaN values".format(cur_header)
            continue

        try:
            df[[cur_header]] = df[[cur_header]].apply(pandas.to_numeric)
            print "column {} has numeric values".format(cur_header)

        except ValueError:
            print "{} cannot be converted to numeric. converting to categorical instead".format(cur_header)
            df[cur_header] = df[cur_header].astype('category')
            df[cur_header] = df[cur_header].cat.codes
    if filter_na_by_rows and filter_type == FILTER_IN:
        print "remove NaN values by row"
        df.dropna(inplace=True)
    df = df.drop(headers[0], 1)

    print "shape : {}".format(df.shape)
    cph = CoxPHFitter()
    cph.fit(df, duration_col='_OS', event_col='_OS_IND', show_progress=True, step_size=0.0001)
    return cph.print_summary()  # access the results using cph.summary
    # print re.sub( '  +', '\t',  cph.summary.to_string(float_format=lambda f: '{:4.4f}'.format(f)))
    # print cph
    # print cph._log_likelihood



