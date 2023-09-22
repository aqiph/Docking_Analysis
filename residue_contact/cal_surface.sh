#!/bin/bash


source /home/soft/amber22/amber.sh

upper_bound=$1
for num in $(seq 2 $upper_bound);do
    cd ./lig$num || exit

    # calculate surface
    cat >>cpp.in<<EOF
  parm comp.prmtop
  trajin comp.inpcrd
  molsurf   :1 out surf$num.dat probe 1
  go
  quit
EOF

    /home/soft/amber22/bin/cpptraj -i cpp.in >> cpp.log

    # results
    echo $num >> ../id_raw.csv
    grep -i 'Row' ../lig$num.pdb  >> ../name_raw.csv
    tail -n 1 surf$num.dat >>../surf_raw.csv

    rm cpp*.in
    cd ..
    
done

