#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 24 15:36:00 2023

@author: guohan
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
font = FontProperties()
font.set_size(12)

from tools import remove_unnamed_columns


def add_docking_score(input_file, input_file_docking_score, id_column_name = 'ID', id_column_name_docking = 'ID',  output_file = None):
    """
    add docking scores from input_file_docking_score to input_file as a new column
    :param input_file: str, file path of the input file
    :param input_file_docking_score: file path of the input docking score file
    :param id_column_name: str, name of the ID column in the main file
    :param id_column_name_docking: str, name of the ID column in the docking score file
    :param output_file: str or None, pre-defined output file name
    """
    # output name
    folder, basename = os.path.split(os.path.abspath(input_file))
    if output_file is None:
        output_file = os.path.splitext(basename)[0] + '_docked.csv'
    output_file = os.path.join(folder, output_file)

    # read files
    df = pd.read_csv(input_file)
    df_docking = pd.read_csv(input_file_docking_score)
    # add docking scores
    df['Docking_Score'] = df[id_column_name].apply(lambda id: get_docking_score(id, df_docking, id_column_name_docking))

    # write to file
    df = df.reset_index(drop=True)
    print('Number of rows:', df.shape[0])
    df = remove_unnamed_columns(df)
    df.to_csv(output_file)


def get_docking_score(id, df_docking, id_column_name_docking):
    """
    helper function for add_docking_score, get docking score for each compound based on id
    :param id: ID of the compound
    :param df_docking: pd.DataFrame, DataFrame object containing docking score
    :param id_column_name_docking: str, name of the ID column in the docking score file
    :return: float, the best docking score for this compound
    """
    df_docking_subset = df_docking[df_docking[id_column_name_docking] == id]
    best_score = df_docking_subset['docking score'].min()
    print(id, best_score)
    return best_score


def plot_activity_dockingScore(input_file, activity_column_name, docking_score_column_list, llabel=False, id_column_name = 'ID'):
    """
    plot
    :param input_file: str, file path of the input file
    :param activity_column_name:
    :param docking_score_column_list:
    """
    COLORS = ['red', 'blue', 'orange']

    # output name
    output_file = os.path.splitext(os.path.abspath(input_file))[0]
    output_file = f'{output_file}_correlation.pdf'

    df = pd.read_csv(input_file)
    activity = df[activity_column_name].tolist()
    # activity = 6.0-np.log10(activity)
    print(activity)

    plt.grid(True)
    for i, docking_score_column_name in enumerate(docking_score_column_list):
        docking_score = df[docking_score_column_name]
        plt.scatter(activity, docking_score, c = COLORS[i], label = docking_score_column_name)
        if llabel:
            IDs = df[id_column_name]
            for id, x, y in zip(IDs, activity, docking_score):
                plt.text(x, y, f'{id}', ha='left', va='bottom', fontsize=9)
    plt.xlabel(activity_column_name, fontproperties=font)
    plt.ylabel('Docking Score', fontproperties=font)
    plt.xticks(fontproperties=font)
    # plt.xlim([0.0, 10.0])
    plt.yticks(fontproperties=font)
    plt.legend()

    plt.title('Activity - Docking score Correlation', fontproperties=font)
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.show()
    plt.close()




if __name__ == '__main__':
    # input_file = 'tests/example.csv'
    # input_file_docking_score = 'tests/example_docking.csv'
    # id_column_name = 'ID'
    # id_column_name_docking = 'Title'
    # add_docking_score(input_file, input_file_docking_score, id_column_name, id_column_name_docking, output_file=None)

    input_file = 'tests/example_docked.csv'
    activity_column_name = 'IC50 Value'
    docking_score_column_list = ['Minimum Docking Score', 'Docking Score Pocket1', 'Docking Score Pocket2']
    plot_activity_dockingScore(input_file, activity_column_name, docking_score_column_list)







