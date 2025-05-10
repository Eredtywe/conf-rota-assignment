import glob
import pandas as pd
from pathlib import Path

from assignment.rota_assign.final_code.get_rotastate import *

aa_chis_dict = {'ARG': 5, 'LYS': 4, 'GLU': 3, 'GLN': 3, 'MET': 3, 'ASP': 2, 'ASN': 2, 'HIS': 2, 'LEU': 2, 'ILE': 2, 'PHE': 2, 'PRO': 2, 'TRP': 2, 'TYR': 2, 'CYS': 1, 'SER': 1, 'THR': 1, 'VAL': 1, 'ALA': 0, 'GLY': 0}


def assign_rota_states(file_path):
	
    in_file_arr = []

    with open(file_path, 'r') as in_file:
        for line in in_file:
            each_line_arr = line.strip().split()
            in_file_arr.append(each_line_arr)

    rota_state_arr = []


    for each_line in in_file_arr:

        if len(each_line[2]) > 1:
            temp_chain = each_line[2][0]
            temp_pos = each_line[2][1:]
            each_line.remove(each_line[2])
            each_line.insert(2, temp_pos)
            each_line.insert(2, temp_chain)
        else:
            pass

        aa = each_line[1]

        ss_b = each_line[12]
        # ss_ub = each_line[21]

        chis_b = [None if float(x) == 0.00 else float(x) for x in each_line[7: 12]]
        # chis_ub = [None if float(x) == 0.00 else float(x) for x in each_line[16: 21]]

        each_line.insert(13, get_aa_rotamer(aa, ss_b, 'bound', each_line[-1], chis_b[0], chis_b[1], chis_b[2], chis_b[3], chis_b[4]))

        
        # if len(each_line) > 20:
        #     each_line.insert(23, get_aa_rotamer(aa, ss_ub, 'unbound', each_line[-1], chis_ub[0], chis_ub[1], chis_ub[2], chis_ub[3], chis_ub[4]))

        rota_state_arr.append(each_line)
    
    return rota_state_arr


if __name__ == '__main__':



    base_path = Path('/Users/savan/Downloads/rota_assign/raw_data')

    pro_rota_files = list((base_path / 'pro_data_b').glob('*.txt'))

    extract_id = lambda f: f.name[-8:-4]

    pro_rota_complexes = [extract_id(f) for f in pro_rota_files]


    for each_id in pro_rota_complexes:
        
        each_file = f'/Users/savan/Downloads/rota_assign/raw_data/pro_data_b/{each_id}.txt'

        rota_assign_arr = assign_rota_states(each_file)

        rota_assign_df_columns = ['pdb_id', 'AA', 'Chain', 'position', 'b_phi', 'b_psi', 'b_omega', 'b_chi1', 'b_chi2', 'b_chi3', 'b_chi4', 'b_chi5', 'b_ss', 'b_rotastate', 'Interface']

        # ['ub_phi', 'ub_psi', 'ub_omega', 'ub_chi1', 'ub_chi2', 'ub_chi3', 'ub_chi4', 'ub_chi5', 'ub_ss', 'b_rotastate']

        rota_assign_df = pd.DataFrame(rota_assign_arr, columns=rota_assign_df_columns)

        rota_assign_df.to_csv(f'/Users/savan/Downloads/rota_assign/final_code/pro_rota_files/{each_id}.csv', index=False)
