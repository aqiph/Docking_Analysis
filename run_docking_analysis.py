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

from main import get_DockingScore, get_DrugLikeCmpds, get_TopScoringCmpds_id, get_MMGBSA, select_Cmpds


if __name__ == '__main__':
    ### Get docking score ###
    input_file = 'tests/test_DockingScore/test_get_DockingScore.csv'
    get_DockingScore(input_file)

    ### Get drug-like compounds ###
    input_file_dockingScore = 'tests/test_DockingScore/test_get_DockingScore_DockingScore.csv'
    input_file_property = 'tests/test_DockingScore/test_get_DrugLikeCmpds_properties.csv'
    get_DrugLikeCmpds(input_file_dockingScore, input_file_property, cutoff_MW=700.0, cutoff_logP=(0.0, 5.0))

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