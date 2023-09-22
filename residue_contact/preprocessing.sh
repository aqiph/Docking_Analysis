#!/bin/bash

# generate .pdb and .mol2 file for protein and ligands, process protein file prot.pdb
source /home/soft/amber22/amber.sh
obabel $1 -O lig.pdb -m
obabel $1 -O lig.mol2 -m
grep -v 'CONECT' lig1.pdb | grep -v 'MASTER' | sed 's/END/TER/g' > prot.pdb

upper_bound=$2
for num in $(seq 2 $upper_bound);do
    # process ligand, generate LIG_raw2.mol2
    mkdir lig$num
    cd ./lig$num || exit
    cp ../lig$num.mol2 ./LIG_raw.mol2
    grep -v 'CONECT' LIG_raw.mol2 | grep -v 'MASTER' | sed 's/UNL/LIG/g' | sed 's/CL/Cl/g' |sed 's/BR/Br/g' > LIG_raw2.mol2
   
    # creat the top crd for lig
    /home/soft/amber22/bin/antechamber -i LIG_raw2.mol2 -fi mol2 -o LIG.mol2 -fo mol2 -c bcc -s 2 -dr no > antechamber.log
    /home/soft/amber22/bin/parmchk2 -i LIG.mol2 -f mol2 -o LIG.frcmod
    rm leap*.in cpp*.in

    # tleap -s leap.in&
    cat >>leap_lig.in<<EOF
  source leaprc.protein.ff19SB
  source leaprc.gaff2
  LIG= loadmol2 LIG.mol2
  loadamberparams LIG.frcmod
  saveoff LIG LIG.lib
  saveamberparm LIG LIG.prmtop LIG.inpcrd
  savepdb LIG LIG_raw.pdb
  quit
EOF
    /home/soft/amber22/bin/tleap -s -f leap_lig.in > leap_lig.log  

    grep -v 'CONECT' LIG_raw.pdb | sed 's/END/TER/g' | sed 's/CL/Cl/g' |sed 's/BR/Br/g' > LIG.pdb
    cat LIG.pdb ../prot.pdb > comp_raw.pdb
    source /home/soft/amber22/amber.sh
    /home/soft/amber22/bin/reduce -Trim -Quiet comp_raw.pdb > comp_raw2.pdb

    # tleap -s leap.in&
    cat >>leap_comp.in<<EOF
  source leaprc.protein.ff19SB
  source leaprc.gaff2
  loadoff LIG.lib 
  loadamberparams LIG.frcmod
  com = loadpdb comp_raw2.pdb
  saveamberparm com comp.prmtop comp.inpcrd
  savepdb com comp.pdb
  quit
EOF
    /home/soft/amber22/bin/tleap -s -f leap_comp.in >leap_comp.log
    cd ..

done