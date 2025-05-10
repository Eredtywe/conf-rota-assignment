#!/usr/bin/env python3
import numpy as np
import sys
import json
import importlib

import build_rotalib

importlib.reload(build_rotalib)

aa_chis_dict = build_rotalib.aa_chis_dict
full_rota_lib_dict_b_I, full_rota_lib_dict_b_N, full_rota_lib_dict_ub_I, full_rota_lib_dict_ub_N, rota_lib_dict_b_I, rota_lib_dict_b_N, rota_lib_dict_ub_I, rota_lib_dict_ub_N = build_rotalib.run_all_files()

def wrap_180(angle):
    return angle - np.ceil(angle / 360.0 - 0.5) * 360.0


def angles_in_valid_range(aa, chis, rotamer, rotastate, rotastate_values):

    width = 30
    temp_possible_states_values = []

    for each_rotastate_value in rotastate_values:
        if all(np.abs(wrap_180(float(chis[each_chi_index]) - float(each_rotastate_value[each_chi_index]))) < width for each_chi_index in range(len(each_rotastate_value) - 1)):
            temp_possible_states_values.append(each_rotastate_value)

    if len(temp_possible_states_values) == 0:
        return [rotamer, 0]
    else:
        sorted_temp_possible_states_values = sorted(temp_possible_states_values, key=lambda x: x[-1])
        return [rotastate + str(round(sorted_temp_possible_states_values[0][-2])), sorted_temp_possible_states_values[0][-1]]



def get_aa_rotamer(aa, ss, b_ub, interface, chi_1=None, chi_2=None, chi_3=None, chi_4=None, chi_5=None):
    
    chis_temp = [chi_1, chi_2, chi_3, chi_4, chi_5]
    chis = [val for val in chis_temp if val is not None]

    aa_chis = aa_chis_dict[aa]

    if len(chis) == 0 and aa in ['ALA', 'GLY']:
        return 'no_chi'
    elif len(chis) < aa_chis:
        return 'MA'
    else:
        pass
    

    if (b_ub == 'bound') and (interface == 'I'):
        lib_dict_keys = full_rota_lib_dict_b_I.keys()
    elif (b_ub == 'bound') and (interface == 'N'):
        lib_dict_keys = full_rota_lib_dict_b_N.keys()
    elif (b_ub == 'unbound') and (interface == 'I'):
        lib_dict_keys = full_rota_lib_dict_ub_I.keys()
    elif (b_ub == 'unbound') and (interface == 'N'):
        lib_dict_keys = full_rota_lib_dict_ub_N.keys()

    if aa not in lib_dict_keys:
        return '!' * aa_chis

    if (b_ub == 'bound') and (interface == 'I'):
        aa_rota_library = full_rota_lib_dict_b_I[aa][ss]
    elif (b_ub == 'bound') and (interface == 'N'):
        aa_rota_library = full_rota_lib_dict_b_N[aa][ss]
    elif (b_ub == 'unbound') and (interface == 'I'):
        aa_rota_library = full_rota_lib_dict_ub_I[aa][ss]
    elif (b_ub == 'unbound') and (interface == 'N'):
        aa_rota_library = full_rota_lib_dict_ub_N[aa][ss]


    final_assign_rotastates = []

    if len(aa_rota_library.keys()) == 0:
        return '!' * aa_chis
    
    for rotastate in aa_rota_library.keys():

        rotamer = '!' * aa_chis
        rotastate_value = aa_rota_library[rotastate]

        final_assign_rotastates.append(angles_in_valid_range(aa, chis, rotamer, rotastate, rotastate_value))


    sorted_temp_possible_states_values = sorted(final_assign_rotastates, key=lambda x: x[-1], reverse=True)
    final_state = sorted_temp_possible_states_values[0][0]

    return final_state