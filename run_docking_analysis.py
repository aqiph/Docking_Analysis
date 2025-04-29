#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 25 10:29:00 2023

@author: guohan

"""

import sys

path_list = sys.path
module_path = '/Users/guohan/Documents/Codes/Docking_Analysis'
if module_path not in sys.path:
    sys.path.append(module_path)
    print('Add module path')

from main import filter_by_dockingScore, filter_by_property, get_TopScoringCmpds_id, get_MMGBSA, select_Cmpds


if __name__ == '__main__':
    ### 1. Apply docking score filter ###
    # input_file_SMILES = 'tests/test_SMILES_file.csv'
    # input_file_dockingScore = 'tests/test_dockingScore_filter.csv'
    # id_column_name = 'ID'
    # filter_by_dockingScore(input_file_SMILES, input_file_dockingScore, id_column_name,
    #                        dockingScore_column_name='docking score', dockingScore_cutoff=-6.9)

    ### 2. Apply property filters ###
    # input_file_SMILES = 'tests/test_SMILES_file.csv'
    # input_file_property = 'tests/test_property_filter.csv'
    # id_column_name = 'ID'
    # property_column_names = ['Docking_Score', 'MW', 'logP', 'HBD', 'HBA', 'TPSA']
    # property_filters = {'MW': lambda x: x <= 650, 'logP': lambda x: x <= 5.5}
    # filter_by_property(input_file_SMILES, input_file_property, id_column_name,
    #                    property_column_names=property_column_names, property_filters=property_filters))

    ### Get unique top-ranking compounds ###
    input_file = 'tests/test_DockingScore/test_get_DockingScore_DockingScore_drugLike.csv'
    # Using top percentage
    # get_TopScoringCmpds_id(input_file, method='percentage', dockingScore_percentage=0.2, total_num_compounds=100)
    # Using docking score cutoff
    get_TopScoringCmpds_id(input_file, method='cutoff', dockingScore_cutoff=-5.0)

    ### Get MM-GBSA dG binding energy ###
    input_file = 'tests/test_MMGBSA/test_get_MMGBSA.csv'
    get_MMGBSA(input_file)

    ### Select compounds based on docking score and MM-GBSA ###
    input_file_dockingScore = 'tests/test_MMGBSA/test_select_Cmpds_DockingScore.csv'
    input_file_mmgbsa = 'tests/test_MMGBSA/test_get_MMGBSA_MMGBSA.csv'
    # select_Cmpds(input_file_dockingScore, input_file_mmgbsa, dockingScore_method='cutoff', mmgbsa_method='cutoff',
    #              merge_method='outer', dockingScore_cutoff = -9.0, mmgbsa_cutoff=-90.0)
    select_Cmpds(input_file_dockingScore, input_file_mmgbsa, dockingScore_method='percentage',
                 mmgbsa_method='percentage', total_num_compounds=100,
                 merge_method='outer', dockingScore_percentage=0.1, mmgbsa_percentage=0.1)