#!/bin/bash


source /home/soft/amber22/amber.sh

upper_bound=$1
for num in $(seq 2 $upper_bound);do
    cd ./lig$num || exit

    # calculate contacts
    cat >>cpp.in<<EOF
  parm comp.prmtop
  trajin comp.inpcrd
  nativecontacts   :1$2 :$3$4 out contact$num.dat distance $5
  go
  quit
EOF

    /home/soft/amber22/bin/cpptraj -i cpp.in >> cpp.log

    # results
    tail -n 1 contact$num.dat>> ../contact$6_raw.csv

    rm cpp*.in
    cd ..
    
done

