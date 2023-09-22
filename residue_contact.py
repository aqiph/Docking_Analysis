#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 4 14:15:00 2023

@author: guohan
"""

import os, subprocess
import pandas as pd
from Bio import PDB
import argparse

# global variable
BASH_PATH = '/home/hanguo02/Code/Docking_Analysis/residue_contact'


def preprocessing(input_ligand_sdf, working_dir_path, num_molecules):
    """
    Data preprocessing in working_dir_path. Generate .mol2 and .pdb file for each compound. Generate files for complex.
    :param input_ligand_sdf: str, path of the input ligand sdf file.
    :param working_dir_path: str, path of the working directory.
    :param num_molecules: int, 'number of ligands.'
    """
    if not os.path.exists(working_dir_path):
        os.makedirs(working_dir_path)
        print(f'Directory created: {working_dir_path}')
    else:
        print(f'Directory already exists: {working_dir_path}')
    os.chdir(working_dir_path)

    # run bash script to generate com.pdb for complex
    print('Starting preprocessing ...')
    subprocess.run([f"{BASH_PATH}/preprocessing.sh", input_ligand_sdf, str(num_molecules+1)])
    print('Preprocessing done.')


def pdb_to_dataframe(pdb_path):
    """
    Convert comp.pdb to DataFrame object
    :param pdb_path: str, path of the protein pdb file.
    :return: pd.DataFrame object
    """
    # read a PDB file, get structure
    PDBParser = PDB.PDBParser(QUIET=True)
    structure = PDBParser.get_structure("ATOM", pdb_path)

    # get structure coordinats
    atom_data = []
    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    atom_info = {
                        "Record_Type": atom.get_full_id()[0],
                        "Atom_Serial_Number": atom.get_serial_number(),
                        "Atom_Name": atom.get_name(),
                        "Amino_Acid": residue.resname,
                        "Chain_ID": chain.id,
                        "Residue_Sequence_Number": residue.id[1],
                        "X_Coordinate": atom.coord[0],
                        "Y_Coordinate": atom.coord[1],
                        "Z_Coordinate": atom.coord[2],
                        "Element": atom.element,
                    }
                    atom_data.append(atom_info)

    # create a DataFrame object
    df = pd.DataFrame(atom_data)

    return df


def convert_residue_number(old_pdb_path, new_pdb_path, old_res_num):
    """
    Get the residue number in new_pdb_path based on the residue number in old_pdb_path
    :param old_pdb_path: str, path of the old protein pdb file before updating.
    :param new_pdb_path: str, path of the new protein pdb file after updating.
    :param old_res_num: int, old residue number.
    :return: int, new residue number.
    """
    # read pdb files
    df_old = pdb_to_dataframe(old_pdb_path)
    df_new = pdb_to_dataframe(new_pdb_path)

    # given a residue number in old .pdb file, get the coordinates of one atom in that residue
    query_res_num = f"Residue_Sequence_Number == '{old_res_num}'"
    df_resi_num = df_old.query(query_res_num)
    x, y, z = df_resi_num['X_Coordinate'].values[0], df_resi_num['Y_Coordinate'].values[0], df_resi_num['Z_Coordinate'].values[0]
    amino_acid_old = df_resi_num['Amino_Acid'].values[0]

    # given coordinates of one atom in the residue, get the residue number in new .pdb file
    query_xyz = f"X_Coordinate == {x} and Y_Coordinate == {y} and Z_Coordinate == {z}"
    df_xyz = df_new.query(query_xyz)
    new_res_num = df_xyz['Residue_Sequence_Number'].values[0]
    amino_acid_new = df_xyz['Amino_Acid'].values[0]

    assert amino_acid_old == amino_acid_new
    print(f'Residue {amino_acid_old}{old_res_num} in new .pdp file is {amino_acid_new}{new_res_num}.')

    return new_res_num


def main(args):
    """

    :return:
    """
    # working directory
    input_file = os.path.abspath(args.input_file)
    folder, _ = os.path.split(input_file)
    input_protein_pdb = os.path.abspath(args.input_protein_pdb)
    input_ligand_sdf = os.path.abspath(args.input_ligand_sdf)
    working_dir_path = os.path.join(folder, args.working_dir)
    num_molecules = args.num_molecules

    # run preprocessing process
    if not os.path.exists(working_dir_path):
        preprocessing(input_ligand_sdf, working_dir_path, num_molecules)

    # change directory
    os.chdir(working_dir_path)

    # read lig_atom_list, res_num, res_atom and dist
    df_input = pd.read_csv(input_file, sep='\s+')
    num_contacts = df_input.shape[0]
    lig_atom_list = df_input['Ligand_Atom'].tolist()
    res_num_list_old = df_input['Residue_Number'].tolist()
    res_atom_list = df_input['Residue_Atom'].tolist()
    dist_list = df_input['Distance'].tolist()
    dist_list = [str(d) for d in dist_list]

    # update residue numbers
    new_pdb_path = f'{working_dir_path}/lig2/comp.pdb'
    res_num_list_new = [convert_residue_number(input_protein_pdb, new_pdb_path, old_res_num) for old_res_num in res_num_list_old]
    res_num_list_new = [str(n) for n in res_num_list_new]

    # run bash script to calculate surface
    subprocess.run([f"{BASH_PATH}/cal_surface.sh", str(num_molecules+1)])

    # run bash script to calculate residue contacts
    df_CONT = pd.DataFrame()
    for i in range(num_contacts):
        subprocess.run([f"{BASH_PATH}/cal_residue_contact.sh", str(num_molecules+1), lig_atom_list[i], res_num_list_new[i], res_atom_list[i], dist_list[i], str(i)])
        df_cont = pd.read_csv(f'contact{i}_raw.csv', header=None, index_col=None, sep='\s+')
        df_cont = pd.DataFrame(df_cont.iloc[:, 1])
        df_cont.columns = [f'Contact{i}']
        df_CONT = pd.concat([df_CONT, df_cont], axis=1)
        os.remove(f'contact{i}_raw.csv')

    # concate results
    df_id = pd.read_csv('id_raw.csv', header=None, index_col=None, sep='\s+')
    df_id.columns = ['ID']
    df_name = pd.read_csv('name_raw.csv', header=None, index_col=None, sep='\s+')
    df_name = pd.DataFrame(df_name.iloc[:, 1])
    df_name.columns = ['Name']
    df_surf = pd.read_csv('surf_raw.csv', header=None, index_col=None, sep='\s+')
    df_surf = pd.DataFrame(df_surf.iloc[:, 1])
    df_surf.columns = ['Surface']
    df = pd.concat([df_id, df_name, df_surf, df_CONT], axis=1)

    for file in ['id_raw.csv', 'name_raw.csv', 'surf_raw.csv']:
        os.remove(file)

    # write output file
    df = df.reset_index(drop=True)
    print('Number of rows in file:', df.shape[0])
    df.to_csv(os.path.join(os.path.dirname(working_dir_path), args.output_file))


def get_parser():
    """
    generate parser
    :return:
    """
    argparser = argparse.ArgumentParser()

    argparser.add_argument('--input_file', default='tests/test_residue_contact/residue_contact_input.csv', type=str, help='Path of the file specifying residue contact info.')
    argparser.add_argument('--input_protein_pdb', default='tests/test_residue_contact/8DCS-prepared.pdb', type=str, help='Path of the protein pdb file.')
    argparser.add_argument('--input_ligand_sdf', default='tests/test_residue_contact/top_molecules.sdf', type=str, help='Path of the ligand sdf file.')
    argparser.add_argument('--working_dir', default='contact_count', type=str, help='Name of the working directory.')
    argparser.add_argument('--num_molecules', default=1, type=int, help='Number of ligands.')
    argparser.add_argument('--output_file', default='properties.csv', type=str, help='Path of the output file.')

    args = argparser.parse_args()
    return args



if __name__ == '__main__':
    # input_ligand_sdf = 'top_molecules.sdf'
    # num_molecules = 2
    # preprocessing(input_ligand_sdf, working_dir = 'contact_count', num_molecules=num_molecules)

    args = get_parser()
    main(args)
