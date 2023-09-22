#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 17 15:07:00 2023

@author: guohan
"""

import os
import pandas as pd
import numpy as np

from tools import remove_unnamed_columns


def get_topScoring_from_clusters(input_file, score_column_name, label_column_name, count=1, singularity_method='top', singularity_label=0):
    """
    Get the top-n scoring compounds from each cluster
    :param input_file: str, path of the input file
    :param score_column_name: str, name of the score column
    :param label_column_name: str, name of the cluster label column
    :param count: int, number of top-scoring compounds from each cluster
    :param singularity_method: str, method to address singularity, allowed values include 'all', 'top' and 'none'
    :param singularity_label: int, label of singularity class
    """
    # input file
    df = pd.read_csv(input_file)
    COLUMNS = df.columns.tolist()

    # extract representative compounds
    df_representative = pd.DataFrame(columns=COLUMNS)
    label_set = set(df[label_column_name])

    if singularity_method == 'all' or singularity_method == 'none':
        label_set.remove(singularity_label)

    for label in label_set:
        df_cluster = df[df[label_column_name] == label]
        df_cluster.sort_values(by=[score_column_name])
        n = min(count, df_cluster.shape[0])
        df_new_representative = df_cluster.iloc[:n]
        df_representative = pd.concat([df_representative, df_new_representative], ignore_index=True, sort=False)

    if singularity_method == 'all':
        df_new_representative = df[df[label_column_name] == singularity_label]
        df_representative = pd.concat([df_representative, df_new_representative], ignore_index=True, sort=False)

    df_representative = df_representative.reset_index(drop=True)
    print('Number of rows in representative file:', df_representative.shape[0])
    df_representative = remove_unnamed_columns(df_representative)
    output_file = f'{os.path.splitext(input_file)[0]}_representative_{df_representative.shape[0]}.csv'
    df_representative.to_csv(output_file)






if __name__ == '__main__':
    input_file = 'tests/test_get_topScoring_from_clusters.csv'
    score_column_name = 'Docking_Score'
    label_column_name = 'MCS Cluster'
    get_topScoring_from_clusters(input_file, score_column_name, label_column_name, count=1, singularity_method='all',
                                 singularity_label=0)



