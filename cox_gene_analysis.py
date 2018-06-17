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
pheno_limit = 135


def cox_gene(test_independently, filtered_out, filtered_in, filter_type, filter_na_by_rows):
    tested_gene_list = "test.txt"
    survival_dataset = load_survival_data("TCGA-SKCM.survival.tsv", survival_list_path=None)
    gene_expression_dataset = load_gene_expression_profile_by_patients(tested_gene_list, "TCGA-SKCM.htseq_counts.tsv")

    expression_survival_integrated = {}
    for cur_expression in gene_expression_dataset[1:]:
        expression_survival_integrated[cur_expression[0]] = cur_expression[pheno_start:pheno_limit]
    for cur_survival in survival_dataset[1:]:
        if expression_survival_integrated.has_key(cur_survival[0]):
            expression_survival_integrated[cur_survival[0]] = list(expression_survival_integrated[cur_survival[0]]) + cur_survival[4:]
    for k, v in dict(expression_survival_integrated).iteritems():
        if len(v) != len(gene_expression_dataset[0][pheno_start:pheno_limit]) + len(survival_dataset[0][4:]):
            expression_survival_integrated.pop(k, None)

    if filter_type==FILTER_IN:
        filtered_list = filtered_in
    else:
        filtered_list = filtered_out
    if test_independently:
        for i, cur in enumerate(filtered_list):
            results = test_cox(gene_expression_dataset, survival_dataset, expression_survival_integrated, [cur]+ OS_FIELDS, filter_type, filter_na_by_rows)
            if i == 0:
                f.write(results)
            else:
                f.write(results[results.index('\n')+1:])
    else:
        f.write(test_cox(gene_expression_dataset, survival_dataset, expression_survival_integrated, filtered_list+ OS_FIELDS, filter_type, filter_na_by_rows))

def test_cox(phenotype_dataset, survival_dataset, pheno_survival_integrated, filtered_list, filter_type, filter_na_by_rows):
    headers = list(phenotype_dataset[0][0:1]) + list(phenotype_dataset[0][pheno_start:pheno_limit]) + list(survival_dataset[0][4:])
    pandas.set_option("mode.use_inf_as_na", True)
    df = pandas.DataFrame(columns=headers, data=[[k] + v for k, v in
                                                 pheno_survival_integrated.iteritems()])  # np.array().astype(np.float32)
    for cur_header in headers[1:]:
        if filter_type == FILTER_IN and cur_header not in filtered_list:
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
    # df[['_OS','_OS_IND', 'age_at_initial_pathologic_diagnosis']] =  df[['_OS','_OS_IND', 'age_at_initial_pathologic_diagnosis']].apply(pandas.to_numeric)
    # df["anatomic_treatment_site"] = df["anatomic_treatment_site"].astype('category')
    # df['anatomic_treatment_site'] = df['anatomic_treatment_site'].cat.codes
    df = df.drop(headers[0], 1)
    # df = pandas.get_dummies(df)
    # print str(df['tx_on_clinical_trial'])
    print "shape : {}".format(df.shape)
    cph = CoxPHFitter()
    cph.fit(df, duration_col='_OS', event_col='_OS_IND', show_progress=False, step_size=0.001)
    return cph.print_summary()  # access the results using cph.summary
    # print re.sub( '  +', '\t',  cph.summary.to_string(float_format=lambda f: '{:4.4f}'.format(f)))
    # print cph
    # print cph._log_likelihood



