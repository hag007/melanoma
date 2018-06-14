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


def main(test_independently, filtered_out, filtered_in, filter_type, filter_na_by_rows):
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



test_independently = False
filtered_out = []
filtered_in = ["ENSG00000004975.10", "ENSG00000010404.16", "ENSG00000016391.9", "ENSG00000016490.14", "ENSG00000042980.11", "ENSG00000049249.7", "ENSG00000052126.13", "ENSG00000057019.14", "ENSG00000064115.9", "ENSG00000065615.12", "ENSG00000068985.4", "ENSG00000074201.7", "ENSG00000075651.14", "ENSG00000083290.18", "ENSG00000084093.14", "ENSG00000086544.2", "ENSG00000087116.12", "ENSG00000088367.19", "ENSG00000088682.12", "ENSG00000100077.13", "ENSG00000100167.18", "ENSG00000101310.13", "ENSG00000104047.13", "ENSG00000104133.13", "ENSG00000105137.11", "ENSG00000105383.13", "ENSG00000106538.8", "ENSG00000108242.11", "ENSG00000109099.12", "ENSG00000109181.10", "ENSG00000113407.12", "ENSG00000115216.12", "ENSG00000118518.14", "ENSG00000120458.8", "ENSG00000121741.15", "ENSG00000124098.9", "ENSG00000124104.17", "ENSG00000124140.11", "ENSG00000124422.10", "ENSG00000124875.8", "ENSG00000125779.20", "ENSG00000125995.14", "ENSG00000126231.12", "ENSG00000129128.11", "ENSG00000129353.13", "ENSG00000129450.7", "ENSG00000129451.10", "ENSG00000131738.8", "ENSG00000132155.10", "ENSG00000132694.17", "ENSG00000133475.15", "ENSG00000134755.13", "ENSG00000134759.12", "ENSG00000134760.5", "ENSG00000134762.15", "ENSG00000136153.18", "ENSG00000136682.13", "ENSG00000136694.8", "ENSG00000137434.10", "ENSG00000137642.11", "ENSG00000137693.12", "ENSG00000137707.12", "ENSG00000137968.15", "ENSG00000138246.14", "ENSG00000140199.10", "ENSG00000140274.12", "ENSG00000141577.12", "ENSG00000142252.9", "ENSG00000142621.18", "ENSG00000143515.15", "ENSG00000143520.6", "ENSG00000143554.12", "ENSG00000144451.17", "ENSG00000145920.13", "ENSG00000147669.9", "ENSG00000147687.15", "ENSG00000148344.10", "ENSG00000153292.14", "ENSG00000153790.10", "ENSG00000154222.13", "ENSG00000156453.12", "ENSG00000158050.4", "ENSG00000158122.10", "ENSG00000158773.13", "ENSG00000159335.14", "ENSG00000159450.11", "ENSG00000159496.13", "ENSG00000159516.8", "ENSG00000161791.12", "ENSG00000162068.1", "ENSG00000162365.10", "ENSG00000163191.5", "ENSG00000163206.5", "ENSG00000163207.6", "ENSG00000163214.19", "ENSG00000163217.1", "ENSG00000163431.12", "ENSG00000164366.3", "ENSG00000164694.15", "ENSG00000164823.8", "ENSG00000164853.8", "ENSG00000165792.16", "ENSG00000165795.19", "ENSG00000165949.11", "ENSG00000166181.11", "ENSG00000166394.13", "ENSG00000166532.14", "ENSG00000166669.12", "ENSG00000166866.11", "ENSG00000167642.11", "ENSG00000167654.16", "ENSG00000167740.8", "ENSG00000167751.11", "ENSG00000167754.11", "ENSG00000167755.12", "ENSG00000167757.12", "ENSG00000167767.12", "ENSG00000167768.4", "ENSG00000168140.4", "ENSG00000168702.15", "ENSG00000168906.11", "ENSG00000169032.8", "ENSG00000169446.5", "ENSG00000169469.8", "ENSG00000169508.6", "ENSG00000169548.3", "ENSG00000169592.13", "ENSG00000170421.10", "ENSG00000170425.3", "ENSG00000170448.10", "ENSG00000170464.8", "ENSG00000170476.14", "ENSG00000170782.3", "ENSG00000171121.15", "ENSG00000171345.12", "ENSG00000171396.11", "ENSG00000171402.13", "ENSG00000171953.14", "ENSG00000172476.3", "ENSG00000172548.13", "ENSG00000172817.3", "ENSG00000172845.12", "ENSG00000174948.5", "ENSG00000175115.10", "ENSG00000175311.6", "ENSG00000175701.9", "ENSG00000175792.10", "ENSG00000176087.13", "ENSG00000176182.5", "ENSG00000178171.9", "ENSG00000178358.4", "ENSG00000178363.4", "ENSG00000178591.6", "ENSG00000178917.13", "ENSG00000178928.7", "ENSG00000179144.4", "ENSG00000179172.8", "ENSG00000179476.6", "ENSG00000181323.7", "ENSG00000182035.10", "ENSG00000182584.4", "ENSG00000183346.6", "ENSG00000183753.8", "ENSG00000184144.8", "ENSG00000184321.1", "ENSG00000184454.6", "ENSG00000184730.9", "ENSG00000185477.4", "ENSG00000185634.10", "ENSG00000185963.12", "ENSG00000186472.18", "ENSG00000186827.9", "ENSG00000186838.12", "ENSG00000186844.5", "ENSG00000187170.4", "ENSG00000187175.5", "ENSG00000187210.11", "ENSG00000188060.6", "ENSG00000188086.11", "ENSG00000188095.4", "ENSG00000188290.9", "ENSG00000188372.13", "ENSG00000188505.4", "ENSG00000188909.4", "ENSG00000188997.6", "ENSG00000189181.4", "ENSG00000189269.11", "ENSG00000189326.4", "ENSG00000189430.11", "ENSG00000196539.3", "ENSG00000196730.11", "ENSG00000196800.5", "ENSG00000197081.11", "ENSG00000197635.8", "ENSG00000198064.12", "ENSG00000198853.10", "ENSG00000203780.9", "ENSG00000203782.5", "ENSG00000203784.2", "ENSG00000203785.7", "ENSG00000203818.6", "ENSG00000204420.7", "ENSG00000204536.12", "ENSG00000204538.3", "ENSG00000204859.10", "ENSG00000205060.9", "ENSG00000205413.6", "ENSG00000206072.11", "ENSG00000206073.9", "ENSG00000215845.9", "ENSG00000223569.4", "ENSG00000235863.3", "ENSG00000237254.2", "ENSG00000240382.3", "ENSG00000241755.1", "ENSG00000244067.2", "ENSG00000256806.3", "ENSG00000261147.1", "ENSG00000264187.1", "ENSG00000265118.4", "ENSG00000268089.2", "ENSG00000269720.1", "ENSG00000281899.1"]
#["melanoma_clark_level_value", "_OS", "_OS_IND"]
filter_type = FILTER_IN
filter_na_by_rows = True
f = file(os.path.join(constants.OUTPUT_DIR,"output_expression", 'w+'))
results = main(test_independently, filtered_out, filtered_in, filter_type, filter_na_by_rows)
f.close()
