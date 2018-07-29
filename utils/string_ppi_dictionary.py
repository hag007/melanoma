from infra import *
import constants




string_ppi_dict = {}

def load_dict():
	u_d = np.array(load_phenotype_data("9606.protein.links.v10.5.txt", phenotype_list_path=constants.DICTIONARIES_DIR))
	u_d_header = u_d[0]
	u_d_content = u_d[1:][u_d[1:, 0].argsort(), :]

	for i, cur_entry in enumerate(u_d_content):
		edge = "{}={}".format(cur_entry[0][5:], cur_entry[1][5:])
		string_ppi_dict[edge] = cur_entry[2]

def get_string_ppi_dict():
	ensg_dict = load_dict()
	return string_ppi_dict

