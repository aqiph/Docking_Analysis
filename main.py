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



def get_DockingScore(input_file):
    """
    Get docking score from raw file
    :param input_file: str, file path of the input docking score file
    """
    # files
    output_file = os.path.splitext(os.path.abspath(input_file))[0] + '_DockingScore.csv'

    df = pd.read_csv(input_file)
    df = pd.DataFrame(df, columns=['title', 'SMILES', 'r_i_docking_score'])
    df.rename(columns={'title': 'Title'}, inplace=True)
    print('Number of rows in docking score file:', df.shape[0])
    # round
    df['Docking_Score'] = df['r_i_docking_score'].apply(lambda score: np.round(score, decimals=3))
    df = pd.DataFrame(df, columns=['Title', 'SMILES', 'Docking_Score'])

    # write output file
    df = df.reset_index(drop=True)
    print('Number of rows:', df.shape[0])
    df = remove_unnamed_columns(df)
    df.to_csv(output_file)
    print('Getting docking score is done.')


def get_DrugLikeCmpds(input_file_dockingScore, input_file_property, cutoff_MW=700.0, cutoff_logP=(0.0, 5.0)):
    """
    Combine drug-like properties and docking score, filter compounds according to properties.
    :param input_file_dockingScore: str, file path of the input docking score file
    :param input_file_property: str, file path of the input property file
    :param cutoff_MW: float, cutoff value for MW
    :param cutoff_logP: tuple of two floats, cutoff values for logP, i.e., (min, max)
    """
    # files
    output_file= os.path.splitext(os.path.abspath(input_file_dockingScore))[0] + '_drugLike.csv'

    df_docking = pd.read_csv(input_file_dockingScore)
    df_docking = pd.DataFrame(df_docking, columns=['Title', 'SMILES', 'Docking_Score'])
    print('Number of rows in docking score file:', df_docking.shape[0])
    df_property = pd.read_csv(input_file_property)
    df_property = pd.DataFrame(df_property, columns=['Title', 'MW', 'logP'])
    print('Number of rows in property file:', df_property.shape[0])

    # merge docking score and properties, and apply filters
    df = pd.merge(df_docking, df_property, how='left', on=['Title'])
    df = pd.DataFrame(df, columns=['Title', 'SMILES', 'Docking_Score', 'MW', 'logP'])
    df = pd.DataFrame(df[df['MW'].apply(lambda MW: MW_filter(MW, cutoff_MW))])
    df = pd.DataFrame(df[df['logP'].apply(lambda logP: logP_filter(logP, cutoff_logP))])

    # write output file
    df = df.reset_index(drop=True)
    print('Number of rows:', df.shape[0])
    df = remove_unnamed_columns(df)
    df.to_csv(output_file)
    print('Getting drug-like compounds is done.')


def MW_filter(MW, cutoff_MW):
    return MW <= cutoff_MW


def logP_filter(logP, cutoff_logP):
    if (logP >= cutoff_logP[0]) and (logP <= cutoff_logP[1]):
        return True
    else:
        return False


def get_TopScoringCmpds_id(input_file, method='percentage', **kwargs):
    """
    Get a list of unique top-ranking compounds based on the given percentage.
    :param input_file: str, file path of the input file
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


def get_MMGBSA(input_file):
    """
    Get MM-GBSA dG binding energy from raw file
    :param input_file: str, file path of the input MM-GBSA file
    """
    # files
    output_file = os.path.splitext(os.path.abspath(input_file))[0] + '_MMGBSA.csv'

    df = pd.read_csv(input_file)
    df = pd.DataFrame(df, columns=['title', 'r_psp_MMGBSA_dG_Bind'])
    df.rename(columns={'title': 'Title'}, inplace=True)
    print('Number of rows in MM-GBSA file:', df.shape[0])
    # round
    df['MMGBSA_dG_Bind'] = df['r_psp_MMGBSA_dG_Bind'].apply(lambda score: np.round(score, decimals=2))
    df = pd.DataFrame(df, columns=['Title', 'MMGBSA_dG_Bind'])

    # write output file
    df = df.reset_index(drop=True)
    print('Number of rows:', df.shape[0])
    df = remove_unnamed_columns(df)
    df.to_csv(output_file)
    print('Getting MM-GBSA dG is done.')


def select_Cmpds(input_file_dockingScore, input_file_mmgbsa, dockingScore_method, mmgbsa_method, merge_method='outer', **kwargs):
    """
    Select compounds according to docking score and MM-GBSA dG
    :param input_file_dockingScore:
    :param input_file_mmgbsa:
    :param dockingScore_method:
    :param mmgbsa_method:
    :param merge_method:
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
    df_dockingScore.sort_values(by=['Docking_Score'], ascending=True, inplace=True)
    df_dockingScore.drop_duplicates(['Title'], keep='first', ignore_index=True, inplace=True)
    df_mmgbsa.sort_values(by=['MMGBSA_dG_Bind'], ascending=True, inplace=True)
    df_mmgbsa.drop_duplicates(['Title'], keep='first', ignore_index=True, inplace=True)
    print('Number of rows in the docking socre file:', df_dockingScore.shape[0])
    print('Number of rows in the MM-GBSA file:', df_mmgbsa.shape[0])

    # merge docking score and MM-GBSA
    df = df_dockingScore.merge(df_mmgbsa, how='left', on=['Title'])
    df = pd.DataFrame(df, columns=['Title', 'SMILES', 'Docking_Score', 'MMGBSA_dG_Bind'])
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
    ### Get docking score ###
    # input_file = 'tests/test_DockingScore/test_get_DockingScore.csv'
    # get_DockingScore(input_file)

    ### Get drug-like compounds ###
    # input_file_dockingScore = 'tests/test_DockingScore/test_get_DockingScore_DockingScore.csv'
    # input_file_property = 'tests/test_DockingScore/test_get_DrugLikeCmpds_properties.csv'
    # get_DrugLikeCmpds(input_file_dockingScore, input_file_property, cutoff_MW=700.0, cutoff_logP=(0.0, 5.0))

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




