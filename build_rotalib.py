import csv
import copy
import json
import pandas as pd

aa_chis_dict = {'ARG': 5, 'LYS': 4, 'GLU': 3, 'GLN': 3, 'MET': 3, 'ASP': 2, 'ASN': 2, 'HIS': 2, 'LEU': 2, 'ILE': 2, 'PHE': 2, 'PRO': 2, 'TRP': 2, 'TYR': 2, 'CYS': 1, 'SER': 1, 'THR': 1, 'VAL': 1, 'ALA': 0, 'GLY': 0}

aaList = ["ARG","ASN","ASP","GLU","GLN","HIS","ILE","LEU","LYS","MET","PHE","PRO","TRP","TYR", "CYS", "SER", "THR", "VAL", "GLY", "ALA"] 





def assign_ptm_state(val):

    """
    Assigns p-t-m state based on predeifned intervals (Both on and off rotamers states).
    :param val: Angle for which ptm state to be assiged. 
    :type val: float
    :returns: ptm state
    :rtype: string
    """

    val = float(val)

    if val > -30 and val <= 30:
        return 'T'
    elif val > 30 and val <= 90:
        return 'p'
    elif val > 90 and val <= 150:
        return 'P'
    elif val > 150 or val <= -150:
        return 't'
    elif val > -150 and val <= -90:
        return 'M'
    elif val > -90 and val <= -30:
        return 'm'
    else:
        return str(int(val))


def chi_state(row):

    """
    Assigns rotamer state based on the chi angels of respective amino acid.
    :param row: information of each amino acid.
    :returns: rotastate of the amino acid.
    :rtype: string
    """
    if row.shape[0] < 9:
        return assign_ptm_state(row['Chi1']) + assign_ptm_state(row['Chi2'])

    else:
        if aa_chis_dict[row['AA']] == 1:
            return assign_ptm_state(row['Chi1'])
        elif aa_chis_dict[row['AA']] == 2:
            return assign_ptm_state(row['Chi1']) + assign_ptm_state(row['Chi2'])
        elif aa_chis_dict[row['AA']] == 3:
            return assign_ptm_state(row['Chi1']) + assign_ptm_state(row['Chi2']) + assign_ptm_state(row['Chi3'])
        elif aa_chis_dict[row['AA']] in [4, 5]:
            return assign_ptm_state(row['Chi1']) + assign_ptm_state(row['Chi2']) + assign_ptm_state(row['Chi3']) + assign_ptm_state(row['Chi4'])
        else:
            return None


def get_rota_stats(file_path):
    """
    Builds dataframes of given rotamer library file (full df(both on anf off combined), on rotamer df, off rotamer df, and array which has the total probabilities of on rotamer states).
    :param file_path: path of the library file (txt/csv/tsv)
    :type file_path: string
    :returns: full_df, on_rotamer_df, off_rotamer_df, prob_arr
    :rtype: df, df, df, list
    """

    if file_path[-3:] == 'txt':
        with open(file_path, 'r') as file:
            lines = file.readlines()

        lines = [line.strip() for line in lines]
        data = [line.split() for line in lines]

        df = pd.DataFrame(data)
        df.columns = ['AA', 'SS', 'Chi1', 'Chi2', 'Chi3', 'Chi4', 'Chi1_std', 'Chi2_std', 'Chi3_std', 'Chi4_std', 'prob']

    elif file_path[-3:] == 'tsv':
        df = pd.read_csv(file_path, skiprows=9, sep='\t')
        df.columns = ['AA', 'SS', 'Chi1', 'Chi2', 'Chi1_std', 'Chi2_std', 'prob']
    
    elif file_path[-3:] == 'csv':
        df = pd.read_csv(file_path, skiprows=9)
        df.columns = ['AA', 'SS', 'Chi1', 'Chi2', 'Chi1_std', 'Chi2_std', 'prob']

    df['rota_state'] = None

    df['rota_state'] = df.apply(chi_state, axis=1)


    prob_arr = []

    on_rota_state_df = pd.DataFrame()
    off_rota_state_df = pd.DataFrame()

    for each_aa in list(set(df['AA'].tolist())):
        for each_ss in list(set(df['SS'].tolist())):
            temp_df = df[(df['AA'] == each_aa) & (df['SS'] == each_ss)]
            if aa_chis_dict[each_aa] != 1:
                prob_arr.append([each_aa, each_ss, sum(temp_df[~temp_df['rota_state'].str[:2].str.contains(r'[A-Z]')]['prob'].astype(float))])

                temp_df_filter_on = temp_df[~temp_df['rota_state'].str[:2].str.contains(r'[A-Z]')]
                temp_df_filter_off = temp_df[temp_df['rota_state'].str[:2].str.contains(r'[A-Z]')]

                on_rota_state_df = pd.concat([on_rota_state_df, temp_df_filter_on])
                off_rota_state_df = pd.concat([off_rota_state_df, temp_df_filter_off])

            else:
                prob_arr.append([each_aa, each_ss, sum(temp_df[~temp_df['rota_state'].str.contains(r'[A-Z]')]['prob'].astype(float))])

                temp_df_filter_on = temp_df[~temp_df['rota_state'].str.contains(r'[A-Z]')]
                temp_df_filter_off = temp_df[temp_df['rota_state'].str.contains(r'[A-Z]')]

                on_rota_state_df = pd.concat([on_rota_state_df, temp_df_filter_on])
                off_rota_state_df = pd.concat([off_rota_state_df, temp_df_filter_off])

    return df, on_rota_state_df, off_rota_state_df, prob_arr


def build_lib(df, rota_lib_file_path=None):
    """
    Builds a 2D-dict which contains the rotamer states of each amino acid based on secondary structure (both on and off).
    """

    final_rota_lib_arr = []

    for i in range(len(df)):

        temp_arr = []

        if aa_chis_dict[df.iloc[i]['AA']] != 1:
            temp_arr.append(df.iloc[i,-1][:-1])
        else:
            temp_arr.append(df.iloc[i,-1])

        for j in range(min(aa_chis_dict[df.iloc[i]['AA']], 4) + 2):
            temp_arr.append(df.iloc[i,j])
            
        temp_arr.append(df.iloc[i,-2])
        final_rota_lib_arr.append(temp_arr)


    all_aa = list(set(df['AA'].tolist()))
    all_ss = list(set(df['SS'].tolist()))


    ss_dep_dict = {}
    final_rota_lib = {}

    temp_dict = {}
    for each_ss in all_ss:
        temp_dict[each_ss] = {}

    ss_dep_dict = temp_dict.copy()

    temp_dict = {}
    for each_aa in all_aa:
        temp_dict[each_aa] = copy.deepcopy(ss_dep_dict)

    final_rota_lib = temp_dict.copy()

    for each_arr in final_rota_lib_arr:
        final_rota_lib[each_arr[1]][each_arr[2]][each_arr[0]] = []

    for each_arr in final_rota_lib_arr:
        formatted_arr = [round(float(value), 3) for value in each_arr[3:]]
        final_rota_lib[each_arr[1]][each_arr[2]][each_arr[0]].append(formatted_arr)

    if rota_lib_file_path is not None:
        with open(rota_lib_file_path, 'w') as file:
            json.dump(final_rota_lib, file, indent=4)

    return final_rota_lib


def run_all_files():

    full_df_b_I, on_rota_state_df_b_I, off_rota_state_df_b_I, prob_arr_b_I = get_rota_stats('/Users/savan/Downloads/rota_assign/bse_lib_files_inter/interface_bound.txt')

    full_df_b_N, on_rota_state_df_b_N, off_rota_state_df_b_N, prob_arr_b_N = get_rota_stats('/Users/savan/Downloads/rota_assign/bse_lib_files_inter/nonInterface_bound.txt')

    full_df_ub_I, on_rota_state_df_ub_I, off_rota_state_df_ub_I, prob_arr_ub_I = get_rota_stats('/Users/savan/Downloads/rota_assign/bse_lib_files_inter/interface_unbound.txt')

    full_df_ub_N, on_rota_state_df_ub_N, off_rota_state_df_ub_N, prob_arr_ub_N = get_rota_stats('/Users/savan/Downloads/rota_assign/bse_lib_files_inter/nonInterface_unbound.txt')


    final_on_rota_lib_b_I = build_lib(on_rota_state_df_b_I, '/Users/savan/Downloads/rota_assign/dump/rota_lib_b_I.txt')

    final_on_rota_lib_b_N = build_lib(on_rota_state_df_b_N, '/Users/savan/Downloads/rota_assign/dump/rota_lib_b_N.txt')

    final_on_rota_lib_ub_I = build_lib(on_rota_state_df_ub_I, '/Users/savan/Downloads/rota_assign/dump/rota_lib_ub_I.txt')

    final_on_rota_lib_ub_N = build_lib(on_rota_state_df_ub_N, '/Users/savan/Downloads/rota_assign/dump/rota_lib_ub_N.txt')


    final_full_rota_lib_b_I = build_lib(full_df_b_I, '/Users/savan/Downloads/rota_assign/dump/full_rota_lib_b_I.txt')

    final_full_rota_lib_b_N = build_lib(full_df_b_N, '/Users/savan/Downloads/rota_assign/dump/full_rota_lib_b_N.txt')

    final_full_rota_lib_ub_I = build_lib(full_df_ub_I, '/Users/savan/Downloads/rota_assign/dump/full_rota_lib_ub_I.txt')

    final_full_rota_lib_ub_N = build_lib(full_df_ub_N, '/Users/savan/Downloads/rota_assign/dump/full_rota_lib_ub_N.txt')


    return final_full_rota_lib_b_I, final_full_rota_lib_b_N, final_full_rota_lib_ub_I, final_full_rota_lib_ub_N, final_on_rota_lib_b_I, final_on_rota_lib_b_N, final_on_rota_lib_ub_I, final_on_rota_lib_ub_N




if __name__ == '__main__':
    
    run_all_files()
