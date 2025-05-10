
import copy
import logging
import pandas as pd
from pathlib import Path


logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(message)s')




def get_pro_bound_info(file_path):

    temp_df = pd.read_csv(file_path)
    temp_up_df = copy.deepcopy(temp_df)

    temp_up_df_cols = ['pdb_id', 'AA', 'Chain', 'position', 'b_phi', 'b_psi', 'b_ss', 'b_rotastate']
    temp_up_df = temp_up_df[temp_up_df_cols]
    temp_up_df.rename(columns={'pdb_id': 'PDB_ID', 'Chain': 'AA_chain', 'position': 'AA_pos'}, inplace=True)
    temp_up_df = temp_up_df.astype(str)

    pro_df_to_merge = copy.deepcopy(temp_up_df)

    return pro_df_to_merge



def get_rna_conf_info(file_path):
    temp_df = pd.read_csv(file_path, compression="gzip", header=0)
    temp_df_columns = ['PDB_ID', 'Chain', 'Nuc', 'Position']

    temp_df[temp_df_columns] = temp_df['step_ID'].str.split('_').apply(lambda x: pd.Series(x[:4]))

    temp_up_df = copy.deepcopy(temp_df)

    temp_up_df = temp_up_df[['PDB_ID', 'Chain', 'Nuc', 'Position', 'NtC', 'CANA', 'confalA']]
    temp_up_df.rename(columns={'Position': 'Nuc_pos', 'Chain': 'Nuc_chain'}, inplace=True)

    temp_up_df['PDB_ID'] = temp_up_df['PDB_ID'].str.upper()
    temp_up_df = temp_up_df.astype(str)

    rna_df_to_merge = copy.deepcopy(temp_up_df)

    return rna_df_to_merge



def get_hb_info(file_path):
    hb_arr = []
    pdb_id = file_path[-8:-4]
    
    for line in open(file_path, "r").readlines():
        if line[1:5].isdigit() and (not ((len(line[6:9].strip()) == 3) and (len(line[20:23].strip()) == 3))):
            d_resNum = int(line[1:5])
            a_resNum = int(line[15:19])
            hb_arr.append([line[6:9].strip(), line[10:13].strip(), line[0], d_resNum, line[20:23].strip(), line[24:27].strip(), line[14], a_resNum, line[28:32].strip()])
        else:
            pass
    
    for each_arr in hb_arr:
        if len(each_arr[0]) != 3:
            each_arr[0], each_arr[4] = each_arr[4], each_arr[0]
            each_arr[1], each_arr[5] = each_arr[5], each_arr[1]
            each_arr[2], each_arr[6] = each_arr[6], each_arr[2]
            each_arr[3], each_arr[7] = each_arr[7], each_arr[3]
        else:
            pass
	
    if not hb_arr:
        logging.info(f'No lines in hb file {pdb_id}')
    else:
        pass
        

    col_names = ['AA', 'AA_atom', 'AA_chain', 'AA_pos', 'Nuc', 'Nuc_atom', 'Nuc_chain', 'Nuc_pos', 'hb_dis']

    hb_df = pd.DataFrame(hb_arr, columns=col_names)
    hb_df.insert(0, 'PDB_ID', pdb_id)
    hb_df['Nuc'] = hb_df['Nuc'].str.strip()
    hb_df = hb_df.astype(str)

    return hb_df



def merge_rota_conf_hb(hb_file_path, pro_file_path, rna_file_path):
    hb_df = get_hb_info(hb_file_path)
    pro_rota_df = get_pro_bound_info(pro_file_path)
    rna_conf_df = get_rna_conf_info(rna_file_path)

    rna_hb_merge_df = hb_df.merge(rna_conf_df, on=['PDB_ID', 'Nuc_chain', 'Nuc_pos', 'Nuc'], how='left')
    final_merge_df = rna_hb_merge_df.merge(pro_rota_df, on=['PDB_ID', 'AA_chain', 'AA_pos', 'AA'], how='left')

    return hb_df, pro_rota_df, rna_conf_df, final_merge_df




if __name__ == '__main__':
	

    base_path = Path('/Users/savan/Downloads/rota_assign/final_code')

    hb_files = list((base_path / 'hb_files').glob('*.hb2'))
    rna_conf_files = list((base_path / 'rna_conf_files').glob('*.csv'))
    pro_rota_files = list((base_path / 'pro_rota_files').glob('*.csv'))

    extract_id = lambda f: f.name[-8:-4]

    hb_complexes = [extract_id(f) for f in hb_files]
    rna_conf_complexes = [extract_id(f) for f in rna_conf_files]
    pro_rota_complexes = [extract_id(f) for f in pro_rota_files]
	
    both_hb_pro_rna = list(set(hb_complexes) & set(pro_rota_complexes) & set(rna_conf_complexes))
	
    not_in_hb_or_pro_or_rna = list((set(hb_complexes) | set(pro_rota_complexes) | set(rna_conf_complexes)) - (set(hb_complexes) & set(pro_rota_complexes) & set(rna_conf_complexes)))
    
    if both_hb_pro_rna:
        pass
    else:
        raise ValueError("Not even one complex has all required files present")

    if not_in_hb_or_pro_or_rna:
        logging.info(f'Protein files which are not merged: {len(set(pro_rota_complexes) - set(both_hb_pro_rna)), set(pro_rota_complexes) - set(both_hb_pro_rna)}')
        logging.info(f'RNA files which are not merged: {len(set(rna_conf_complexes) - set(both_hb_pro_rna)), set(rna_conf_complexes) - set(both_hb_pro_rna)}')
        logging.info(f'One file of {not_in_hb_or_pro_or_rna} is missing')
    else:
        pass
	
    for each_complex in both_hb_pro_rna:
		
        hb_pdb_file_path = f'/Users/savan/Downloads/rota_assign/final_code/hb_files/{each_complex}.hb2'
        rna_conf_file_path = f'/Users/savan/Downloads/rota_assign/final_code/rna_conf_files/{each_complex}.csv'  
        pro_rota_file_path = f'/Users/savan/Downloads/rota_assign/final_code/pro_rota_files/{each_complex}.csv'

        hb_df, pro_rota_df, rna_conf_df, merge_df = merge_rota_conf_hb(hb_pdb_file_path, pro_rota_file_path, rna_conf_file_path)
        
        hb_df.to_csv(f'/Users/savan/Downloads/rota_assign/final_code/hb_lines_files/{each_complex}.csv', index=False)
        merge_df.to_csv(f'/Users/savan/Downloads/rota_assign/final_code/hb_pro_rna_merge_files/{each_complex}.csv', index=False)