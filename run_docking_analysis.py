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

from util import add_docking_score, plot_activity_dockingScore


if __name__ == '__main__':
    input_file = 'tests/example.csv'
    input_file_docking_score = 'tests/example_docking.csv'
    id_column_name = 'ID'
    id_column_name_docking = 'Title'
    add_docking_score(input_file, input_file_docking_score, id_column_name, id_column_name_docking, output_file=None)

    input_file = 'tests/example_docked.csv'
    activity_column_name = 'IC50 Value'
    docking_score_column_list = ['Minimum Docking Score', 'Docking Score Pocket1', 'Docking Score Pocket2']
    plot_activity_dockingScore(input_file, activity_column_name, docking_score_column_list)