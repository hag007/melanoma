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


def cox_phenotype(test_independently, filtered_out, filtered_in, filter_type, filter_na_by_rows, phenotype_file_name, survival_file_name):
    phenotype_dataset = load_phenotype_data(phenotype_file_name=phenotype_file_name, phenotype_list_path=None)
    survival_dataset = load_survival_data(survival_file_name, survival_list_path=None)
    pheno_survival_integrated = {}
    for cur_pheno in phenotype_dataset[1:]:
        pheno_survival_integrated[cur_pheno[0]] = cur_pheno[pheno_start:pheno_limit]
    for cur_survival in survival_dataset[1:]:
        if pheno_survival_integrated.has_key(cur_survival[0]):
            pheno_survival_integrated[cur_survival[0]] = pheno_survival_integrated[cur_survival[0]] + cur_survival[4:]
    for k, v in dict(pheno_survival_integrated).iteritems():
        if len(v) != len(phenotype_dataset[0][pheno_start:pheno_limit]) + len(survival_dataset[0][4:]):
            pheno_survival_integrated.pop(k, None)

    if filter_type==FILTER_IN:
        filtered_list = filtered_in
    else:
        filtered_list = filtered_out
    if test_independently:
        for i, cur in enumerate(filtered_list):
            results = test_cox(phenotype_dataset, survival_dataset, pheno_survival_integrated, [cur]+ OS_FIELDS, filter_type, filter_na_by_rows)
            if i == 0:
                f.write(results)
            else:
                f.write(results[results.index('\n')+1:])
    else:
        f.write(test_cox(phenotype_dataset, survival_dataset, pheno_survival_integrated, filtered_list+ OS_FIELDS, filter_type, filter_na_by_rows))

def test_cox(phenotype_dataset, survival_dataset, pheno_survival_integrated, filtered_list, filter_type, filter_na_by_rows):
    headers = phenotype_dataset[0][0:1] + phenotype_dataset[0][pheno_start:pheno_limit] + survival_dataset[0][4:]
    pandas.set_option("mode.use_inf_as_na", True)
    df = pandas.DataFrame(columns=headers, data=[[k] + v for k, v in
                                                 pheno_survival_integrated.iteritems()])  # np.array().astype(np.float32)
    decode_categorical_values(df)
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


def decode_categorical_values(df):
    df.replace('', np.nan, inplace=True)
    ### encode ###
    df.replace('stage 0', 0, inplace=True)
    df.replace('stage i', 1, inplace=True)
    df.replace('stage ia', 2, inplace=True)
    df.replace('stage ib', 3, inplace=True)
    df.replace('stage ic', 4, inplace=True)
    df.replace('stage ii', 5, inplace=True)
    df.replace('stage iia', 6, inplace=True)
    df.replace('stage iib', 7, inplace=True)
    df.replace('stage iic', 8, inplace=True)
    df.replace('stage iii', 9, inplace=True)
    df.replace('stage iiia', 10, inplace=True)
    df.replace('stage iiib', 11, inplace=True)
    df.replace('stage iiic', 12, inplace=True)
    df.replace('stage iv', 13, inplace=True)
    df.replace('stage iva', 14, inplace=True)
    df.replace('stage ivb', 15, inplace=True)
    df.replace('stage ivc', 16, inplace=True)
    df.replace('stage v', 17, inplace=True)
    df.replace('stage va', 18, inplace=True)
    df.replace('stage vb', 19, inplace=True)
    df.replace('stage vc', 20, inplace=True)
    df.replace('i/ii nos', np.nan, inplace=True)
    df.replace('not reported', np.nan, inplace=True)
    df.replace('I', 0, inplace=True)
    df.replace('II', 1, inplace=True)
    df.replace('III', 2, inplace=True)
    df.replace('IV', 3, inplace=True)
    df.replace('V', 4, inplace=True)
    df.replace('M0', 0, inplace=True)
    df.replace('M1', 1, inplace=True)
    df.replace('M1a', 2, inplace=True)
    df.replace('M1b', 3, inplace=True)
    df.replace('M1c', 4, inplace=True)
    df.replace('N0', 0, inplace=True)
    df.replace('N1', 1, inplace=True)
    df.replace('N1a', 2, inplace=True)
    df.replace('N1b', 3, inplace=True)
    df.replace('N2', 4, inplace=True)
    df.replace('N2a', 5, inplace=True)
    df.replace('N2b', 6, inplace=True)
    df.replace('N2c', 7, inplace=True)
    df.replace('N3', 8, inplace=True)
    df.replace('NX', np.nan, inplace=True)
    df.replace('T0', 0, inplace=True)
    df.replace('T1', 1, inplace=True)
    df.replace('T1a', 2, inplace=True)
    df.replace('T1b', 3, inplace=True)
    df.replace('T2', 4, inplace=True)
    df.replace('T2a', 5, inplace=True)
    df.replace('T2b', 6, inplace=True)
    df.replace('T3', 7, inplace=True)
    df.replace('T3a', 8, inplace=True)
    df.replace('T3b', 9, inplace=True)
    df.replace('T4', 10, inplace=True)
    df.replace('T4a', 11, inplace=True)
    df.replace('T4b', 12, inplace=True)
    df.replace('Tis', np.nan, inplace=True)
    df.replace('TX', np.nan, inplace=True)
    ### end of encode ###

