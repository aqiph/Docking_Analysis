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
import seaborn as sns
sns.set_palette('hls')

from tools import remove_unnamed_columns


def add_DockingScore(input_file, input_file_docking_score, id_column_name = 'ID', id_column_name_docking = 'ID',  output_file = None):
    """
    Add docking scores from input_file_docking_score to input_file as a new column
    :param input_file: str, path of the input file
    :param input_file_docking_score: path of the input docking score file
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
    df['Docking_Score'] = df[id_column_name].apply(lambda id: get_DockingScore(id, df_docking, id_column_name_docking))

    # write to file
    df = df.reset_index(drop=True)
    print('Number of rows:', df.shape[0])
    df = remove_unnamed_columns(df)
    df.to_csv(output_file)


def get_DockingScore(id, df_docking, id_column_name_docking):
    """
    helper function for add_DockingScore, get docking score for each compound based on id
    :param id: ID of the compound
    :param df_docking: pd.DataFrame, DataFrame object containing docking score
    :param id_column_name_docking: str, name of the ID column in the docking score file
    :return: float, the best docking score for this compound
    """
    df_docking_subset = df_docking[df_docking[id_column_name_docking] == id]
    best_score = df_docking_subset['Docking_Score'].min()
    # print(id, best_score)
    return best_score


def plot_activity_DockingScore(input_file, activity_column_name, docking_score_column_list, llabel=False, id_column_name = 'ID'):
    """
    Plot correlation between docking score and activity
    :param input_file: str, path of the input file
    :param activity_column_name: str, name of the activity column
    :param docking_score_column_list: list of str, list of the score column
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


def plot_DockingScoreDistribution(input_file_list, docking_score_column_name, labels, cutoff=None):
    """
    Plot docking score distribution
    :param input_file_list: list of strs, paths of the input files
    :param docking_score_column_name: str, name of the score column
    :param labels: list of str, label of docking score for legend
    :return:
    """
    COLORS = ['red', 'orange', 'blue']

    # output name
    output_file = os.path.splitext(os.path.abspath(input_file_list[0]))[0]
    output_file = f'{output_file}_distribution.pdf'

    # plot distribution
    Scores = []
    for input_file in input_file_list:
        df = pd.read_csv(input_file)
        score = df[docking_score_column_name].tolist()
        Scores.append(score)

    for i, score in enumerate(Scores):
        sns.histplot(score, bins=50, kde=True, edgecolor='black', stat='density', alpha=0.3, label=labels[i], color=COLORS[i])

    if cutoff is not None:
        plt.axvline(x=cutoff, color='black', linestyle='--', linewidth=2)

    plt.xlim(-10, -3)
    plt.xlabel('Docking score')
    plt.ylabel('Frequency')
    plt.title('Docking Score Distribution')

    plt.legend()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.show()
    plt.close()


def combine_files(input_file_list, columns=None, **kwargs):
    """
    Combine input files
    :param input_file_list: list of strs, paths of the input files
    :param columns: list of strs, list of column names in the output file
    """
    # files
    output_file = os.path.splitext(os.path.abspath(input_file_list[0]))[0] + '_combined.csv'

    # read files
    df_list = []
    for file in input_file_list:
        df_list.append(pd.read_csv(file))

    # concat files
    df = pd.concat(df_list, ignore_index=True, sort=False)
    if columns is None:
        columns = df.columns.tolist()
    df = pd.DataFrame(df, columns = columns)

    flag_deduplication = kwargs.get('flag_deduplication', False)
    if flag_deduplication:
        df.drop_duplicates(columns, keep='first', ignore_index=True, inplace=True)

    # write output file
    df = df.reset_index(drop=True)
    print('Number of rows:', df.shape[0])
    df = remove_unnamed_columns(df)
    df.to_csv(output_file)








if __name__ == '__main__':
    # input_file = 'tests/test_add_DockingScore.csv'
    # input_file_docking_score = 'tests/test_add_DockingScore_dockingFile.csv'
    # id_column_name = 'ID'
    # id_column_name_docking = 'Title'
    # add_DockingScore(input_file, input_file_docking_score, id_column_name, id_column_name_docking, output_file=None)

    # input_file = 'tests/test_plotting/test_plot_activity_DockingScore.csv'
    # activity_column_name = 'IC50 Value'
    # docking_score_column_list = ['Minimum Docking Score', 'Docking Score Pocket1', 'Docking Score Pocket2']
    # plot_activity_DockingScore(input_file, activity_column_name, docking_score_column_list)

    input_file = ['tests/test_plotting/test_plot_DockingScoreDistribution_example0.csv',
                  'tests/test_plotting/test_plot_DockingScoreDistribution_example1.csv',
                  'tests/test_plotting/test_plot_DockingScoreDistribution_example2.csv']
    docking_score_column_name = 'Docking_Score'
    labels = ['R2_S0', 'R2_S1', 'R2_S2']
    cutoff = -6.9
    plot_DockingScoreDistribution(input_file, docking_score_column_name, labels, cutoff)







