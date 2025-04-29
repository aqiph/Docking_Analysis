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



def filter_by_dockingScore(input_file_SMILES, input_file_dockingScore, id_column_name='ID',
                           dockingScore_column_name='r_i_docking_score', dockingScore_cutoff=0.0):
    """
    Filter compounds based on docking score
    :param input_file_SMILES: str, path of the input SMILES file
    :param input_file_dockingScore: str, path of the input docking score file
    :param id_column_name: str, name of the ID column in input_file_dockingScore
    :param dockingScore_column_name: str, name of the docking score column
    :param dockingScore_cutoff: float, docking score cutoff
    :return: None
    """
    # files
    output_file = os.path.splitext(os.path.abspath(input_file_dockingScore))[0] + '_DockingScore'

    df_SMILES = pd.read_csv(input_file_SMILES)
    # df_SMILES = pd.DataFrame(df_SMILES, columns=['ID', 'SMILES', 'Cleaned_SMILES'])
    print('Number of rows in SMILES file:', df_SMILES.shape[0])
    df = pd.read_csv(input_file_dockingScore)
    df.rename(columns={id_column_name: 'ID', dockingScore_column_name: 'Docking_Score'}, inplace=True)
    df = pd.DataFrame(df, columns=['ID', 'Docking_Score'])
    print('Number of rows in docking score file:', df.shape[0])

    # round
    df['Docking_Score'] = df['Docking_Score'].apply(lambda score: np.round(score, decimals=3))
    # sort and deduplicate
    df.sort_values(by=['Docking_Score'], ascending=True, inplace=True)
    df = df.drop_duplicates(['ID'], keep='first', ignore_index=True)
    # filter
    df_filtered = df[df['Docking_Score'] <= dockingScore_cutoff]
    df_filtered = pd.DataFrame(df_filtered, columns=['ID', 'Docking_Score'])
    # merge
    df_filtered = pd.merge(df_filtered, df_SMILES, how='right', on=['ID'])
    # df_filtered = pd.DataFrame(df_filtered, columns=['ID', 'SMILES', 'Cleaned_SMILES', 'Docking_Score'])

    # write output file
    df_filtered = df_filtered.reset_index(drop=True)
    print('Number of rows in filtered docking score file:', df_filtered.shape[0])
    df_filtered = remove_unnamed_columns(df_filtered)
    df_filtered.to_csv(f'{output_file}_{df_filtered.shape[0]}.csv')
    print('Applying docking score filter is done.')


def filter_by_property(input_file_SMILES, input_file_property, id_column_name='ID',
                       property_column_names=None, property_filters=None):
    """
    Filter compounds based on property cutoff
    :param input_file_SMILES: str, path of the input SMILES file
    :param input_file_property: str, path of the input property file
    :param id_column_name: str, name of the ID column in input_file_property
    :param property_column_names: list of strs or None, names of the property columns
    :param property_filters: dict or None, dict of functions for property filters
    :return: None
    """
    if property_column_names is None:
        property_column_names = []
    if property_filters is None:
        property_filters = {'MW':lambda x: x <= 650, 'logP':lambda x: x <= 5.5}

    # files
    output_file = os.path.splitext(os.path.abspath(input_file_property))[0] + '_Property'

    df_SMILES = pd.read_csv(input_file_SMILES)
    # df_SMILES = pd.DataFrame(df_SMILES, columns=['ID', 'SMILES', 'Cleaned_SMILES'])
    print('Number of rows in SMILES file:', df_SMILES.shape[0])
    df = pd.read_csv(input_file_property)
    df.rename(columns={id_column_name:'ID'}, inplace=True)
    df = pd.DataFrame(df, columns=['ID']+property_column_names)
    print('Number of rows in property file:', df.shape[0])

    # filter
    for column, filter in property_filters.items():
        try:
            df = df[df[column].apply(filter)]
        except Exception:
            print(f'Error: Filter for {column} column is not applied.')
            continue
    # merge
    df = pd.merge(df, df_SMILES, how='left', on=['ID'])
    # df = pd.DataFrame(df, columns=['ID', 'SMILES', 'Cleaned_SMILES']+property_column_names)

    # write output file
    df = df.reset_index(drop=True)
    print('Number of rows in filtered property file:', df.shape[0])
    df = remove_unnamed_columns(df)
    df.to_csv(f'{output_file}_{df.shape[0]}.csv')
    print('Applying property filters is done.')


def get_TopScoringCmpds_id(input_file, method='percentage', **kwargs):
    """
    Get a list of unique top-ranking compounds based on the given percentage.
    :param input_file: str, path of the input file
    :param method: str, method using which to get top-ranking molecules, allowed values include 'cutoff' and 'percentage'
    :param total_num_compounds: int, total number of unique compounds. Use the number of unique compounds from input_file, if not specified.
    :param dockingScore_percentage: float, percentage of the top-ranking compounds, if method == 'percentage'
    :param dockingScore_cutoff: float, cutoff for docking score, if method == 'cutoff'
    """
    # files
    output_file = os.path.splitext(os.path.abspath(input_file))[0] + '_topScoringID.csv'
    output_file_txt = os.path.splitext(os.path.abspath(input_file))[0] + '_topScoringID.txt'
    df = pd.read_csv(input_file)
    print('Number of rows in the original file:', df.shape[0])

    # sort values and deduplicate
    df.sort_values(by=['Docking_Score'], ascending=True, inplace=True)
    df_deduplicated = df.drop_duplicates(['Title'], keep='first', ignore_index=True)

    # get top-ranked molecules
    if method == 'percentage':
        num_cmpds = kwargs.get('total_num_compounds', df_deduplicated.shape[0])
        dockingScore_percentage = kwargs.get('dockingScore_percentage', 0.001)
        num_selected_ideal = int(np.round(float(num_cmpds) * dockingScore_percentage, decimals=0))
        dockingScore_cutoff = df_deduplicated.loc[(num_selected_ideal-1), 'Docking_Score']
        print(f'The ideal number of compounds is {num_selected_ideal}, the docking score cutoff is {dockingScore_cutoff}.')
    elif method == 'cutoff':
        dockingScore_cutoff = kwargs.get('dockingScore_cutoff', 0.0)
        print(f'The docking score cutoff is {dockingScore_cutoff}.')
    else:
        raise Exception('Error: Invalid method for getting top compound IDs.')
    df_id = pd.DataFrame(df_deduplicated[df_deduplicated['Docking_Score'].apply(lambda score: score <= dockingScore_cutoff)])

    # write output file
    df_id = df_id.reset_index(drop=True)
    print('Number of rows:', df_id.shape[0])
    df_id = remove_unnamed_columns(df_id)
    df_id.to_csv(output_file)

    # write txt file
    id_list = df_id.Title.values.tolist()
    with open(output_file_txt, 'w') as fh:
        for i in range(len(id_list)):
            fh.write(f'{id_list[i]}\n')
    print('Getting top-scoring compounds is done.')


def get_MMGBSA(input_file, input_file_mmgbsa):
    """
    Get MM-GBSA dG binding energy from raw file
    :param input_file_mmgbsa: str, path of the input MM-GBSA file
    """
    # files
    output_file = os.path.splitext(os.path.abspath(input_file))[0] + '_MMGBSA'

    df = pd.read_csv(input_file)
    df_mmgbsa = pd.read_csv(input_file_mmgbsa)
    df_mmgbsa.rename(columns={'title': 'ID', 'r_psp_MMGBSA_dG_Bind':'MMGBSA_dG_Bind'}, inplace=True)
    df_mmgbsa = pd.DataFrame(df_mmgbsa, columns=['ID', 'MMGBSA_dG_Bind'])
    print('Number of rows in MM-GBSA file:', df_mmgbsa.shape[0])
    # remove duplicates
    df_mmgbsa.sort_values(by=['MMGBSA_dG_Bind'], ascending=True, inplace=True)
    df_mmgbsa.drop_duplicates(['ID'], keep='first', ignore_index=True, inplace=True)
    # round
    df_mmgbsa.loc[:, 'MMGBSA_dG_Bind'] = df_mmgbsa['MMGBSA_dG_Bind'].apply(lambda score: np.round(score, decimals=2))
    # merge
    df = df.merge(df_mmgbsa, how='left', on=['ID'])

    # write output file
    df = df.reset_index(drop=True)
    print('Number of rows:', df.shape[0])
    df = remove_unnamed_columns(df)
    df.to_csv(f'{output_file}_{df.shape[0]}.csv')
    print('Getting MM-GBSA dG is done.')


def select_Cmpds(input_file_dockingScore, input_file_mmgbsa, dockingScore_method, mmgbsa_method, merge_method='outer', **kwargs):
    """
    Select compounds according to docking score and MM-GBSA dG
    :param input_file_dockingScore: str, path of the input docking score file
    :param input_file_mmgbsa: str, path of the input MM-GBSA file
    :param dockingScore_method: str, method using which to get top-ranking molecules based on docking score, allowed values include 'cutoff' and 'percentage'
    :param mmgbsa_method: str, method using which to get top-ranking molecules based on MM-GBSA dG, allowed values include 'cutoff' and 'percentage'
    :param merge_method: str, method to merge results selected by docking score and MM-GBSA dG
    :param total_num_compounds: int, total number of unique compounds. Use the number of unique compounds from input_file, if not specified.
    :param dockingScore_percentage: float, percentage of the top-ranking compounds according to docking score, if dockingScore_method == 'percentage'
    :param dockingScore_cutoff: float, cutoff for docking score, if dockingScore_method == 'cutoff'
    :param mmgbsa_percentage: float, percentage of the top-ranking compounds according to MM-GBSA dG, if mmgbsa_method == 'percentage'
    :param mmgbsa_cutoff: float, cutoff for MM-GBSA dG, if mmgbsa_method == 'cutoff'
    """
    # files
    output_file = os.path.splitext(os.path.abspath(input_file_dockingScore))[0]
    df_dockingScore = pd.read_csv(input_file_dockingScore)
    df_mmgbsa = pd.read_csv(input_file_mmgbsa)

    # sort values and deduplicate
    df_dockingScore.sort_values(by=['Docking_Score_P1'], ascending=True, inplace=True)
    df_dockingScore.drop_duplicates(['ID'], keep='first', ignore_index=True, inplace=True)
    df_mmgbsa.sort_values(by=['MMGBSA_dG_Bind'], ascending=True, inplace=True)
    df_mmgbsa.drop_duplicates(['ID'], keep='first', ignore_index=True, inplace=True)
    print('Number of rows in the docking socre file:', df_dockingScore.shape[0])
    print('Number of rows in the MM-GBSA file:', df_mmgbsa.shape[0])

    # merge docking score and MM-GBSA
    df = df_dockingScore.merge(df_mmgbsa, how='left', on=['ID'])
    # df = pd.DataFrame(df, columns=['Title', 'SMILES', 'Docking_Score', 'MMGBSA_dG_Bind'])
    num_cmpds = kwargs.get('total_num_compounds', df.shape[0])

    # get top-ranked molecules based on docking score
    df.sort_values(by=['Docking_Score'], ascending=True, inplace=True)
    if dockingScore_method == 'percentage':
        dockingScore_percentage = kwargs.get('dockingScore_percentage', 0.001)
        num_selected_ideal = int(np.round(float(num_cmpds) * dockingScore_percentage, decimals=0))
        dockingScore_cutoff = df.iloc[(num_selected_ideal - 1), 2]
        print(f'Docking score: The ideal number of compounds is {num_selected_ideal}, the docking score cutoff is {dockingScore_cutoff}.')
    elif dockingScore_method == 'cutoff':
        dockingScore_cutoff = kwargs.get('dockingScore_cutoff', 0.0)
        print(f'Docking score: The docking score cutoff is {dockingScore_cutoff}.')
    else:
        raise Exception('Error: Invalid method for selecting compounds based on docking score.')
    df_dockingScore_selected = pd.DataFrame(df[df['Docking_Score'].apply(lambda score: score <= dockingScore_cutoff)])
    print(f'Docking score: The real number of compounds is {df_dockingScore_selected.shape[0]}.')

    # get top-ranked molecules based on MM-GBSA
    df.sort_values(by=['MMGBSA_dG_Bind'], ascending=True, inplace=True)
    if mmgbsa_method == 'percentage':
        mmgbsa_percentage = kwargs.get('mmgbsa_percentage', 0.001)
        num_selected_ideal = int(np.round(float(num_cmpds) * mmgbsa_percentage, decimals=0))
        mmgbsa_cutoff = df.iloc[(num_selected_ideal - 1), 3]
        print(f'MM-GBSA: The ideal number of compounds is {num_selected_ideal}, the MM-GBSA cutoff is {mmgbsa_cutoff}.')
    elif mmgbsa_method == 'cutoff':
        mmgbsa_cutoff = kwargs.get('mmgbsa_cutoff', 0.0)
        print(f'MM-GBSA: The docking score cutoff is {mmgbsa_cutoff}.')
    else:
        raise Exception('Error: Invalid method for selecting compounds based on MM-GBSA.')
    df_mmgbsa_selected = pd.DataFrame(df[df['MMGBSA_dG_Bind'].apply(lambda score: score <= mmgbsa_cutoff)])
    print(f'MM-GBSA: The real number of compounds is {df_mmgbsa_selected.shape[0]}.')

    # combine selected compounds
    df_selected = df_dockingScore_selected.merge(df_mmgbsa_selected, how=merge_method, on=['Title', 'SMILES', 'Docking_Score', 'MMGBSA_dG_Bind'])

    # write output file
    df = df.reset_index(drop=True)
    print('Number of unique compounds:', df.shape[0])
    df = remove_unnamed_columns(df)
    df.to_csv(output_file + '_unique.csv')

    df_selected = df_selected.reset_index(drop=True)
    print('Number of selected compounds:', df_selected.shape[0])
    df_selected = remove_unnamed_columns(df_selected)
    df_selected.to_csv(output_file + '_selected.csv')
    print('Compound selection is done.')


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
        df_cluster = df_cluster.sort_values(by=[score_column_name])
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


def cut_file(input_file):
    output_file = os.path.splitext(os.path.abspath(input_file))[0] + '_cutted.csv'
    df = pd.read_csv(input_file)

    # df.sort_values(by=['Title'], inplace=True)
    df = df.iloc[:110]

    # write output file
    df = df.reset_index(drop=True)
    print('Number of rows:', df.shape[0])
    df = remove_unnamed_columns(df)
    df.to_csv(output_file)



if __name__ == '__main__':
    ### 1. Apply docking score filter ###
    # input_file_SMILES = 'tests/test_SMILES_file.csv'
    # input_file_dockingScore = 'tests/test_dockingScore_filter.csv'
    # id_column_name = 'ID'
    # filter_by_dockingScore(input_file_SMILES, input_file_dockingScore, id_column_name,
    #                     dockingScore_column_name='docking score', dockingScore_cutoff=-6.9)

    ### 2. Apply property filters ###
    # input_file_SMILES = 'tests/test_SMILES_file.csv'
    # input_file_property = 'tests/test_property_filter.csv'
    # id_column_name = 'ID'
    # property_column_names = ['Docking_Score', 'MW', 'logP', 'HBD', 'HBA', 'TPSA']
    # property_filters = {'MW':lambda x: x <= 650, 'logP':lambda x: x <= 5.5}
    # filter_by_property(input_file_SMILES, input_file_property, id_column_name,
    #                 property_column_names=property_column_names, property_filters=property_filters)

    ### Get unique top-ranking compounds ###
    # input_file = 'tests/test_DockingScore/test_get_DockingScore_DockingScore_drugLike.csv'
    # Using top percentage
    # get_TopScoringCmpds_id(input_file, method='percentage', dockingScore_percentage=0.2, total_num_compounds=100)
    # Using docking score cutoff
    # get_TopScoringCmpds_id(input_file, method='cutoff', dockingScore_cutoff=-5.0)

    ### Get _TopCmpds.maegz file
    ### Command on cluster: structsubset --title-file *_topScoringID.txt *_pv.maegz *_TopCmpds.maegz

    ### Get MM-GBSA dG binding energy ###
    # input_file = 'tests/test_MMGBSA/test_get_MMGBSA.csv'
    # get_MMGBSA(input_file)

    ### Select compounds based on docking score and MM-GBSA ###
    # input_file_dockingScore = 'tests/test_MMGBSA/test_select_Cmpds_DockingScore.csv'
    # input_file_mmgbsa = 'tests/test_MMGBSA/test_get_MMGBSA_MMGBSA.csv'
    # select_Cmpds(input_file_dockingScore, input_file_mmgbsa, dockingScore_method='cutoff', mmgbsa_method='cutoff',
    #              merge_method='outer', dockingScore_cutoff = -9.0, mmgbsa_cutoff=-90.0)
    # select_Cmpds(input_file_dockingScore, input_file_mmgbsa, dockingScore_method='percentage', mmgbsa_method='percentage', total_num_compounds=100,
    #              merge_method='outer', dockingScore_percentage=0.1, mmgbsa_percentage=0.1)

    ### Select compounds from clusters ###
    input_file = 'tests/test_DockingScore/test_get_topScoring_from_clusters.csv'
    score_column_name = 'Docking_Score'
    label_column_name = 'MCS Cluster'
    get_topScoring_from_clusters(input_file, score_column_name, label_column_name, count=1, singularity_method='all',
                                 singularity_label=0)


    ### Testing ###
    # input_file = 'tests/prime_mmgbsa_SP_6RR7_P1_TopCmps-out.csv'
    # cut_file(input_file)




